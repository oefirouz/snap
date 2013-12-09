// pushrelabelmain.cpp : Defines the entry point for the console application.

#include "stdafx.h"
#include <ctime>

#define INF (1 << 30)

const static float global_update_freq = 0.5;
const static int alpha = 6;
const static int beta = 12;

static int num_nodes = 0;
static int num_edges = 0;

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


static int work_since_update = INF;
static int s;
static int t;
static PNEANet G;
static Arc *arcs;
static Node *nodes;
static Layer *layers;


static inline void add_to_active_layer(int i) {
  int layer = nodes[i].height;
  nodes[i].next_active = layers[layer].first_active;
  layers[layer].first_active = i;
  if (min_active > nodes[i].height) { min_active = nodes[i].height; }
  if (max_active < nodes[i].height) { max_active = nodes[i].height; }
}


static inline void remove_from_active_layer(int i) {
  int layer = nodes[i].height;
  layers[layer].first_active = nodes[i].next_active;
}


/*void add_to_inactive_layer(int i) {
  int layer = nodes[i].height;
  nodes[i].next_inactive = layers[layer].first_inactive;
  if (nodes[i].next_inactive >= 0) {
    nodes[nodes[i].next_inactive].prev_inactive = i;
  }
  layers[layer].first_inactive = i;
}


void remove_from_inactive_layer(int i) {
  int layer = nodes[i].height;
  if (layers[layer].first_inactive == i) {
    layers[layer].first_inactive = nodes[i].next_inactive;
  } else {
    nodes[nodes[i].prev_inactive].next_inactive = nodes[i].next_inactive;
  }

  if (nodes[i].next_inactive >= 0) {
    nodes[nodes[i].next_inactive].prev_inactive = nodes[i].prev_inactive;
  }
}*/


void net_from_dimacs(char *filename) {
  FILE *f = fopen(filename, "r");
  if (f == NULL) { printf("Couldn't find file!\n"); return; }
  int a, b, c, e_id, max_edges, max_nodes;
  char z, buf[1024];
  int *capacities;
  while (fgets(buf, 1024, f) != NULL) {
    if (buf[0] == 'p') {
      sscanf(buf, "p max %d %d", &max_nodes, &max_edges);
      num_nodes = max_nodes;
      capacities = (int *) calloc(2*max_edges, sizeof(int));
    } else if (buf[0] == 'a') {
      sscanf(buf, "%c %d %d %d", &z, &a, &b, &c);
      if (!G->IsNode(a-1)) { G->AddNode(a-1); }
      if (!G->IsNode(b-1)) { G->AddNode(b-1); }
      if (!G->IsEdge(a-1,b-1)) { G->AddEdge(a-1,b-1); num_edges++; }
      if (!G->IsEdge(b-1,a-1)) { G->AddEdge(b-1,a-1); num_edges++; }
      e_id = G->GetEId(a-1,b-1);
      capacities[e_id] = c;
    }
  }
  arcs = (Arc *) calloc(num_edges, sizeof(Arc));
  for (int i = 0; i < num_edges; ++i) {
    arcs[i].capacity = capacities[i];
  }
  free(capacities);
  fclose(f);
}


static inline void push(int u, int v) {
  int e_id = G->GetEId(u,v);
  int rev_e_id = G->GetEId(v,u);
  int f = arcs[e_id].flow;
  int delta = arcs[e_id].capacity - f;
  int ex = nodes[u].excess;
  if (ex > delta) { ex = delta; }
  arcs[e_id].flow = f + ex;
  arcs[rev_e_id].flow = -(f + ex);
  nodes[u].excess -= ex;
  nodes[v].excess += ex;
}


static inline void global_relabel() {
  static int *bfs_queue = (int *) malloc(G->GetNodes()*sizeof(int));
  for (int i = 0; i < num_nodes; ++i) {
    nodes[i].height = num_nodes; //INF;
    layers[i].first_active = -1;
  }

  max_active = 0;
  min_active = num_nodes; //TODO: maybe min active isn't very good

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
        if (arcs[e_id].capacity - arcs[e_id].flow > 0) {
          nodes[prev_node].height = nodes[cur_node].height + 1;
          bfs_queue[right] = prev_node;
          right++;

          if (nodes[prev_node].excess > 0) {
            add_to_active_layer(prev_node);
          } else {

          }
        }
      }
    }
  }
}


static inline void relabel(int u) {
  work_since_update++;
  TNEANet::TNodeI NI = G->GetNI(u);
  int min_neighbor = INF;
  for (int i = 0; i < NI.GetOutDeg(); ++i) {
    int v = NI.GetOutNId(i);
    int e_id = G->GetEId(u,v);
    if (arcs[e_id].capacity - arcs[e_id].flow > 0) {
      if (min_neighbor > nodes[v].height) { min_neighbor = nodes[v].height; }
    }
  }
  nodes[u].height = min_neighbor + 1;
}


static inline void gap(int missing_height) {
  for (int i = missing_height + 1; i < num_nodes; ++i) {
    int cur = layers[i].first_inactive;
    while (cur >= 0) {
      nodes[cur].height = num_nodes;
      cur = nodes[cur].next_inactive;
    }
    layers[i].first_inactive = -1;
  }
  max_active = missing_height - 1;
}


static inline void discharge(int u) {
  TNEANet::TNodeI NI = G->GetNI(u);
  while (1) {
    for (int i = 0; i < NI.GetOutDeg(); ++i) {
      int v = NI.GetOutNId(i);
      int e_id = G->GetEId(u,v);
      if (arcs[e_id].capacity - arcs[e_id].flow > 0) {
        if (nodes[u].height == nodes[v].height + 1) {
          if (v != s && v != t && nodes[v].excess == 0) {
            //TODO remove_from_inactive_layer(v);
            add_to_active_layer(v);
          }
          push(u, v);
          if (nodes[u].excess == 0) {
            break;
          }
        }
      }
    }

    if (nodes[u].excess > 0) {
      int old_height = nodes[u].height;
      relabel(u);
      if (layers[old_height].first_inactive < 0 && layers[old_height].first_active < 0) {
        gap(old_height);
      }
      //TODO: Gap
      if (nodes[u].height == num_nodes) { //IMPORTANT
        break;
      }
    } else {
      //TODO: add_to_inactive_layer(u);
      break;
    }
  }
}


int push_relabel() {
  //push_relabel_init();
  int u, v;

  nodes = (Node *) calloc(num_nodes, sizeof(Node));
  layers = (Layer *) calloc(num_nodes + 1, sizeof(Layer));


  for (int i = 0; i < num_nodes; ++i) {
    nodes[i].height = 1;
    if (i == s) {
      nodes[i].height = num_nodes;
    } else if (i == t) {
      nodes[i].height = 0;
    }
  }

  for (int i = 0; i < num_nodes + 1; ++i) {
    layers[i].first_active = -1;
    //TODO: layers[i].first_inactive = -1;
  }
  for (int i = 0; i < num_nodes; ++i) {
    //TODO nodes[i].prev_inactive = -1;
    //TODO nodes[i].next_inactive = -1;
    nodes[i].next_active = -1;
  }

  min_active = num_nodes;
  max_active = 0;

  nodes[s].height = num_nodes;
  nodes[s].excess = INF;

  TNEANet::TNodeI NI = G->GetNI(s);
  for (int i = 0; i < NI.GetOutDeg(); ++i) {
    v = NI.GetOutNId(i);
    push(s, v);
    add_to_active_layer(v);
  }

  for (int i = 0; i < num_nodes; ++i) {
    if (i != s && i != t) {
      if (nodes[i].excess > 0) {
        add_to_active_layer(i);
      } else {
        //add_to_inactive_layer(i);
      }
    }
  }


  while (max_active >= min_active) {
    if (layers[max_active].first_active < 0) {
      max_active--;
    } else {
      u = layers[max_active].first_active;
      remove_from_active_layer(u);
      discharge(u);
      if (work_since_update > global_update_freq*num_nodes) {
        work_since_update = 0;
        global_relabel();
      }
    }
  }

  return nodes[t].excess;
}


void usage() { printf("USAGE: pushrelabel filename\n"); }


int main(int argc, char* argv[]) {
  if (argc <= 1) { usage(); return 0; }
  char *filename = argv[1];
  G = PNEANet::New();
  net_from_dimacs(filename);
  s = 0;
  t = num_nodes - 1;

  clock_t start = clock();
  int flow = push_relabel();
  clock_t end = clock();
  printf("%d\t%d\t", num_nodes, num_edges);
  printf("%d\t", flow);
  printf("%f\t", ((float)end - start)/CLOCKS_PER_SEC);

  return 0;
}
