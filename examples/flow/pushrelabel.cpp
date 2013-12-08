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

typedef struct node {
  int excess;
  int height;
  int next_active;
  int next_inactive;
  int prev_inactive;
} Node;


typedef struct layer {
  int first_active;
  int first_inactive;
} Layer;


void add_to_active_layer(Layer *layers, int layer, Node *nodes, int i) {
  nodes[i].next_active = layers[layer].first_active;
  layers[layer].first_active = i;
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


void remove_From_inactive_layer(Layer *layers, int layer, Node *nodes, int i) {
  if (layers[layer].first_active == i) {
    layers[layer].first_inactive = nodes[i].next_inactive;
    nodes[nodes[i].next_inactive].prev_inactive = -layer;
  } else {
    nodes[nodes[i].prev_inactive].next_inactive = nodes[i].next_inactive;
    nodes[nodes[i].next_inactive].prev_inactive = nodes[i].prev_inactive;
  }
}


/*typedef struct arc {
  int capacity;
  int flow;
  //int available;
} Arc;*/


void net_from_dimacs(PNEANet G, char *filename, int **capacities, int **flows) {
  FILE *f = fopen(filename, "r");
  if (f == NULL) { printf("Couldn't find file!\n"); return; }
  int a, b, c, e_id, max_edges, max_nodes;
  char z, buf[1024];
  while (fgets(buf, 1024, f) != NULL) {
    if (buf[0] == 'p') {
      sscanf(buf, "p max %d %d", &max_nodes, &max_edges);
      NUM_NODES = max_nodes;
      *capacities = (int *) calloc(2*max_edges, sizeof(int));
      *flows = (int *) calloc(2*max_edges, sizeof(int));
    } else if (buf[0] == 'a') {
      sscanf(buf, "%c %d %d %d", &z, &a, &b, &c);
      if (!G->IsNode(a-1)) { G->AddNode(a-1); }
      if (!G->IsNode(b-1)) { G->AddNode(b-1); }
      if (!G->IsEdge(a-1,b-1)) { G->AddEdge(a-1,b-1); NUM_EDGES++; }
      if (!G->IsEdge(b-1,a-1)) { G->AddEdge(b-1,a-1); NUM_EDGES++; }
      e_id = G->GetEId(a-1,b-1);
      (*capacities)[e_id] = c;
    }
  }
  fclose(f);
}


static inline void push(PNEANet G, int u, int v, Node *nodes, int *capacities, int *flows) {
  int e_id = G->GetEId(u,v);
  int rev_e_id = G->GetEId(v,u);
  int f = flows[e_id];

  int delta = capacities[e_id] - f;
  int ex = nodes[u].excess;
  int amt = MIN( ex, delta );
  flows[e_id] = f + amt;
  flows[rev_e_id] = -(f + amt);
  nodes[u].excess -= amt;
  nodes[v].excess += amt;
}


static inline void global_relabel(PNEANet G, int t, Node *nodes, int *capacities, int *flows) {
  //TODO: can be done without linear pass, also use a constant buffer
  static int *bfs_queue = (int *) malloc(G->GetNodes()*sizeof(int));
  for (int i = 0; i < NUM_NODES; ++i) {
    nodes[i].height = NUM_NODES; //INF;
  }
  nodes[t].height = 0;
  TNEANet::TNodeI NI = G->GetNI(t);
  int left = 0;
  int right = 1;
  bfs_queue[0] = t;
  while (left != right) {
    int cur_node = bfs_queue[left];
    left++;
    NI = G->GetNI(cur_node);
    for (int i = 0; i < NI.GetInDeg(); ++i) {
      int prev_node = NI.GetInNId(i);
      if (nodes[prev_node].height > nodes[cur_node].height + 1) {
        int e_id = G->GetEId(prev_node, cur_node);
        if (capacities[e_id] - flows[e_id] > 0) {
          nodes[prev_node].height = nodes[cur_node].height + 1;
          bfs_queue[right] = prev_node;
          right++;
        }
      }
    }
  }
}


static inline void relabel(PNEANet G, int u, int t, Node *nodes, int *capacities, int *flows) {
  static int counter = 0;
  counter++;
  if (counter > global_update_freq*NUM_NODES) {
    counter = 0;
    global_relabel(G, t, nodes, capacities, flows);
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


int push_relabel(PNEANet G, int s, int t, int *capacities, int *flows) {
  int u, v;

  Node *nodes = (Node *) calloc(NUM_NODES, sizeof(Node));
  bool *in_queue = (bool *) calloc(NUM_NODES, sizeof(bool)); //TODO: Deprecated

  int min_neighbor;

  global_relabel(G, t, nodes, capacities, flows); // CANNOT DO AFTER HEIGHT/EXCESS
  nodes[s].height = NUM_NODES;
  nodes[s].excess = INF;

  TNEANet::TNodeI NI = G->GetNI(s);
  TSnapQueue<int> node_queue(NUM_NODES);
  for (int i = 0; i < NI.GetOutDeg(); ++i) {
    v = NI.GetOutNId(i);
    push(G, s, v, nodes, capacities, flows);
    node_queue.Push(v);
    in_queue[v] = 1;
  }

  while (node_queue.Empty() == 0) {
    min_neighbor = INF;
    u = node_queue.Top();
    NI = G->GetNI(u);
    for (int i = 0; i < NI.GetOutDeg(); ++i) {
      v = NI.GetOutNId(i);
      int e_id = G->GetEId(u,v);
      if (capacities[e_id] - flows[e_id] > 0) {
        if (nodes[u].height > nodes[v].height) {
          push(G, u, v, nodes, capacities, flows);
          if (in_queue[v] == 0 && v != s && v != t) {
            in_queue[v] = 1;
            node_queue.Push(v);
          }
        }
        if (min_neighbor > nodes[v].height) { min_neighbor = nodes[v].height; }
      }
    }

    if (nodes[u].excess == 0) {
      in_queue[u] = 0;
      node_queue.Pop();
    } else if (nodes[u].excess > 0) {
      nodes[u].height = min_neighbor + 1;
      relabel(G, u, t, nodes, capacities, flows);
    }
  }

  return nodes[t].excess;
}


void usage() { printf("USAGE: pushrelabel filename\n"); }


int main(int argc, char* argv[]) {
  if (argc <= 1) { usage(); return 0; }
  char *filename = argv[1];
  PNEANet G = PNEANet::New();

  int *capacities, *flows;

  //clock_t start_init = clock();
  net_from_dimacs(G, filename, &capacities, &flows);
  //clock_t end_init = clock();

  //printf("%s\t", filename);
  //printf("Init:%f\t", ((float)end_init - start_init)/CLOCKS_PER_SEC);
  //fflush(stdout);
  clock_t start = clock();
  int flow = push_relabel(G, 0, G->GetNodes()-1, capacities, flows);
  clock_t end = clock();
  printf("%d\t%d\t", NUM_NODES, NUM_EDGES);
  printf("%d\t", flow);
  printf("%f\t", ((float)end - start)/CLOCKS_PER_SEC);

  return 0;
}
