// pushrelabelmain.cpp : Defines the entry point for the console application.

#include "stdafx.h"
#include <ctime>

#define INF (1 << 30)
#define MIN(A,B) ((A < B) ? (A) : (B))

const static float global_update_freq = 0.5;
const static int alpha = 6;
const static int beta = 12;

static int NUM_NODES = 0;
static int NUM_EDGES = 0;

static int max_distance = 0;
static int max_active = 0;
static int min_active = 0;

typedef struct node {
  int excess;
  int height;
  int next_active;
  int next_inactive;
  int prev_inactive;
  // int firstarc
} Node;


typedef struct arc {
  int capacity;
  int flow;
  //int available;
  //int firstNode
  //int secondNode
} Arc;


typedef struct layer {
  int first_active;
  int first_inactive;
} Layer;


void add_to_active_layer(Layer *layers, int layer, Node *nodes, int i) {
  nodes[i].next_active = layers[layer].first_active;
  layers[layer].first_active = i;
  if (min_active > nodes[i].height) { min_active = nodes[i].height; }
  if (max_active < nodes[i].height) { max_active = nodes[i].height; }
}


void remove_from_active_layer(Layer *layers, int layer, Node *nodes, int i) {
  layers[layer].first_active = nodes[i].next_active;
}


//TODO: I'm pretty sure layer is the same as node.height
void add_to_inactive_layer(Layer *layers, int layer, Node *nodes, int i) {
  nodes[i].next_inactive = layers[layer].first_inactive;
  nodes[nodes[i].next_inactive].prev_inactive = i;
  nodes[i].prev_inactive = -layer; // flag indicating layer struct
  layers[layer].first_inactive = i;
}


void remove_from_inactive_layer(Layer *layers, int layer, Node *nodes, int i) {
  if (layers[layer].first_active == i) {
    layers[layer].first_inactive = nodes[i].next_inactive;
    nodes[nodes[i].next_inactive].prev_inactive = -layer;
  } else {
    nodes[nodes[i].prev_inactive].next_inactive = nodes[i].next_inactive;
    nodes[nodes[i].next_inactive].prev_inactive = nodes[i].prev_inactive;
  }
}


void net_from_dimacs(PNEANet G, char *filename, Arc **arcs) {
  FILE *f = fopen(filename, "r");
  if (f == NULL) { printf("Couldn't find file!\n"); return; }
  int a, b, c, e_id, max_edges, max_nodes;
  char z, buf[1024];
  int *capacities;
  while (fgets(buf, 1024, f) != NULL) {
    if (buf[0] == 'p') {
      sscanf(buf, "p max %d %d", &max_nodes, &max_edges);
      NUM_NODES = max_nodes;
      capacities = (int *) calloc(2*max_edges, sizeof(int));
    } else if (buf[0] == 'a') {
      sscanf(buf, "%c %d %d %d", &z, &a, &b, &c);
      if (!G->IsNode(a-1)) { G->AddNode(a-1); }
      if (!G->IsNode(b-1)) { G->AddNode(b-1); }
      if (!G->IsEdge(a-1,b-1)) { G->AddEdge(a-1,b-1); NUM_EDGES++; }
      if (!G->IsEdge(b-1,a-1)) { G->AddEdge(b-1,a-1); NUM_EDGES++; }
      e_id = G->GetEId(a-1,b-1);
      capacities[e_id] = c;
    }
  }
  *arcs = (Arc *) calloc(NUM_EDGES, sizeof(Arc));
  for (int i = 0; i < NUM_EDGES; ++i) {
    (*arcs)[i].capacity = capacities[i];
  }
  free(capacities);
  fclose(f);
}


static inline void push(PNEANet G, int u, int v, Node *nodes, Arc *arcs) {
  int e_id = G->GetEId(u,v);
  int rev_e_id = G->GetEId(v,u);
  int f = arcs[e_id].flow;

  int delta = arcs[e_id].capacity - f;
  int ex = nodes[u].excess;
  int amt = MIN( ex, delta );
  arcs[e_id].flow = f + amt;
  arcs[rev_e_id].flow = -(f + amt);
  nodes[u].excess -= amt;
  nodes[v].excess += amt;
}


static inline void global_relabel(PNEANet G, int t, Node *nodes, Arc *arcs, Layer *layers) {
  //TODO: globals are bad!
  static int *bfs_queue = (int *) malloc(G->GetNodes()*sizeof(int));
  for (int i = 0; i < NUM_NODES; ++i) {
    nodes[i].height = NUM_NODES; //INF;
  }
  nodes[t].height = 0;
  TNEANet::TNodeI NI = G->GetNI(t);


  //TODO: doesn't need to do all layers, can be more efficient later after gaps
  for (int i = 0; i <= max_distance; ++i) {
    layers[i].first_active = -1;
    layers[i].first_inactive = -1;
  }

  max_distance = max_active = 0;
  min_active = NUM_NODES;


  int left = 0;
  int right = 1;
  bfs_queue[0] = t;
  while (left != right) {
    int cur_node = bfs_queue[left++];
    NI = G->GetNI(cur_node);
    for (int i = 0; i < NI.GetInDeg(); ++i) {
      int prev_node = NI.GetInNId(i);
      if (nodes[prev_node].height > nodes[cur_node].height + 1) {
        int e_id = G->GetEId(prev_node, cur_node);
        if (arcs[e_id].capacity - arcs[e_id].flow > 0) {
          nodes[prev_node].height = nodes[cur_node].height + 1;

          if (max_distance < nodes[prev_node].height) {
            max_distance = nodes[prev_node].height;
          }

          if (nodes[prev_node].excess > 0) {
            add_to_active_layer(layers, nodes[prev_node].height, nodes, prev_node);
          } else {
            add_to_inactive_layer(layers, nodes[prev_node].height, nodes, prev_node);
          }

          bfs_queue[right++] = prev_node;
        }
      }
    }
  }
}


static inline void gap_relabel(Layer *layers, Node *nodes, int empty_distance) {
  int cutoff = empty_distance - 1;
  for (int i = empty_distance + 1; i < max_distance; ++i) {
    int node_index = layers[i].first_inactive;
    while (node_index >= 0) {
      nodes[node_index].height = NUM_NODES;
      node_index = nodes[node_index].next_inactive;
    }
    layers[i].first_inactive = -1;
  }
  max_distance = max_active = cutoff;
}


static inline void relabel(PNEANet G, int u, int t, Node *nodes, Arc *arcs, Layer *layers) {
  static int counter = 0;
  counter++;
  if (counter > global_update_freq*NUM_NODES) {
    counter = 0;
    global_relabel(G, t, nodes, arcs, layers);
  }
  /*TNEANet::TNodeI NI = G->GetNI(u);
  int min_neighbor = INF;
  for (int i = 0; i < NI.GetOutDeg(); ++i) {
    int v = NI.GetOutNId(i);
    int e_id = G->GetEId(u,v);
    if (capacities[e_id] - flows[e_id] > 0) {
      if (min_neighbor > h[v]) { min_neighbor = h[v]; }
    }
  }
  h[u] = min_neighbor + 1;
  */
}

void discharge(PNEANet G, Layer *layers, Node *nodes, Arc *arcs) { return; }


int push_relabel(PNEANet G, int s, int t, Arc *arcs) {
  int u, v;

  Node *nodes = (Node *) calloc(NUM_NODES, sizeof(Node));
  Layer *layers = (Layer *) calloc(NUM_NODES + 1, sizeof(Layer));
  for (int i = 0; i <= NUM_NODES; ++i) {
    layers[i].first_active = -1;
    layers[i].first_inactive = -1;
  }
  for (int i = 0; i < NUM_NODES; ++i) {
    nodes[i].prev_inactive = -1;
    nodes[i].next_inactive = -1;
    nodes[i].next_active = -1;
  }

  //bool *in_queue = (bool *) calloc(NUM_NODES, sizeof(bool)); //TODO: Deprecated

  int min_neighbor;
  min_active = NUM_NODES;
  max_active = 0;

  //global_relabel(G, t, nodes, arcs, layers); // CANNOT DO AFTER HEIGHT/EXCESS
  nodes[s].height = NUM_NODES;
  nodes[s].excess = INF;

  TNEANet::TNodeI NI = G->GetNI(s);
  //TSnapQueue<int> node_queue(NUM_NODES);
  for (int i = 0; i < NI.GetOutDeg(); ++i) {
    v = NI.GetOutNId(i);
    push(G, s, v, nodes, arcs);
    //node_queue.Push(v);
    //in_queue[v] = 1;
  }
  max_active = 0;
  min_active = NUM_NODES;
  max_distance = NUM_NODES - 1;

  for (int i = 0; i < NUM_NODES; ++i) {
    if (i != s && i != t) {
      if (nodes[i].excess > 0) {
        add_to_active_layer(layers, nodes[i].height, nodes, i);
      } else if (nodes[i].height < NUM_NODES) {
        add_to_inactive_layer(layers, nodes[i].height, nodes, i);
      }
    }
  }
  printf("\n%d %d\n", max_active, min_active);

  while (max_active >= min_active) {
  printf("\n%d %d\n", max_active, min_active);
    if (layers[max_active].first_active < 0) {
      max_active--;
    } else {
      u = layers[max_active].first_active;
      remove_from_active_layer(layers, nodes[u].height, nodes, u);
      discharge(G, layers, nodes, arcs);

      /*
    //while (node_queue.Empty() == 0) {
      min_neighbor = INF;
      //u = node_queue.Top();
      NI = G->GetNI(u);
      for (int i = 0; i < NI.GetOutDeg(); ++i) {
        v = NI.GetOutNId(i);
        int e_id = G->GetEId(u,v);
        if (arcs[e_id].capacity - arcs[e_id].flow > 0) {
          if (nodes[u].height > nodes[v].height) {
            push(G, u, v, nodes, arcs);
            if (v != s && v != t && nodes[v].excess == 0) {
              remove_from_inactive_layer(layers, nodes[v].height, nodes, v);
              add_to_active_layer(layers, nodes[v].height, nodes, v);
            //if (in_queue[v] == 0 && v != s && v != t) {
              //in_queue[v] = 1;
              //node_queue.Push(v);


            }
          }
          if (min_neighbor > nodes[v].height) { min_neighbor = nodes[v].height; }
        }
      }

      if (nodes[u].excess == 0) {
        //in_queue[u] = 0;
        //node_queue.Pop();
      } else if (nodes[u].excess > 0) {
        int old_height = nodes[u].height;
        nodes[u].height = min_neighbor + 1;
        if (layers[old_height].first_inactive < 0 && layers[old_height].first_active < 0) {
          gap_relabel(layers, nodes, old_height);
        }
        //relabel(G, u, t, nodes, arcs, layers);
      }*/
    }
  }

  return nodes[t].excess;
}


void usage() { printf("USAGE: pushrelabel filename\n"); }


int main(int argc, char* argv[]) {
  if (argc <= 1) { usage(); return 0; }
  char *filename = argv[1];
  PNEANet G = PNEANet::New();

  Arc *arcs;

  //clock_t start_init = clock();
  net_from_dimacs(G, filename, &arcs);
  printf("NODES: %d EDGES: %d\n", NUM_NODES, NUM_EDGES);
  //clock_t end_init = clock();

  //printf("%s\t", filename);
  //printf("Init:%f\t", ((float)end_init - start_init)/CLOCKS_PER_SEC);
  //fflush(stdout);
  clock_t start = clock();
  int flow = push_relabel(G, 0, G->GetNodes()-1, arcs);
  clock_t end = clock();
  printf("%d\t%d\t", NUM_NODES, NUM_EDGES);
  printf("%d\t", flow);
  printf("%f\t", ((float)end - start)/CLOCKS_PER_SEC);

  return 0;
}
