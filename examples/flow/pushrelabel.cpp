// pushrelabelmain.cpp : Defines the entry point for the console application.

#include "stdafx.h"
#include <ctime>

#define INF (1 << 30)

static const float global_update_freq = 0.5;
static const int alpha = 6;
static const int beta = 12;

static int num_nodes;
static int num_edges;
static int normalized_size;   //alpha * num_nodes + num_edges/2

static int max_active;
static int min_active;
static int max_distance;

#include <list>

typedef struct node {
  int excess;
  int height;
  std::list<int>::iterator list_ptr;
} Node;


typedef struct arc {
  int capacity;
  int flow;
} Arc;


typedef struct layer {
  std::list<int> active;
  std::list<int> inactive;
} Layer;


static int work_since_update;
static int s;
static int t;
static PNEANet G;
static Arc *arcs;
static Node *nodes;
static Layer *layers;


void net_from_dimacs(char *filename) {
  FILE *f = fopen(filename, "r");
  if (f == NULL) { printf("Couldn't find file!\n"); exit(0); }
  int a, b, c, e_id, max_edges, max_nodes;
  char z, buf[1024];
  int *capacities;
  //TODO: source/sink defined in DIMACS
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


static inline void add_to_active_layer(int i) {
  int layer = nodes[i].height;
  //static bool flip = 0;
  //if (flip == 0) {
    layers[layer].active.push_front(i);
    nodes[i].list_ptr = layers[layer].active.begin();
  //} else {
  /*  layers[layer].active.push_back(i);
    nodes[i].list_ptr = layers[layer].active.end();
    --nodes[i].list_ptr;
  }
  flip = 1 - flip;*/
  if (min_active > layer) { min_active = layer; }
  if (max_active < layer) { max_active = layer; }
}


static inline void remove_from_active_layer(int i) {
  int layer = nodes[i].height;
  layers[layer].active.erase(nodes[i].list_ptr);
}


static inline void add_to_inactive_layer(int i) {
  int layer = nodes[i].height;
  layers[layer].inactive.push_front(i);
  nodes[i].list_ptr = layers[layer].inactive.begin();
}


static inline void remove_from_inactive_layer(int i) {
  int layer = nodes[i].height;
  layers[layer].inactive.erase(nodes[i].list_ptr);
}


static void global_relabel() {
  static int *bfs_queue = (int *) malloc(num_nodes*sizeof(int));

  for (int i = 0; i < num_nodes; ++i) { nodes[i].height = num_nodes;  }
  for (int i = 0; i <= max_distance; ++i) {
    layers[i].active.clear();
    layers[i].inactive.clear();
  }

  max_active = max_distance = 0;
  min_active = num_nodes;
  nodes[t].height = 0;

  TNEANet::TNodeI NI = G->GetNI(t);
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
          bfs_queue[right++] = prev_node;

          if (max_distance < nodes[prev_node].height) {
            max_distance = nodes[prev_node].height;
          }

          if (nodes[prev_node].excess > 0) {
            add_to_active_layer(prev_node);
          } else {
            add_to_inactive_layer(prev_node);
          }
        }
      }
    }
  }
}


static inline void push(int u, int v, int e_id) {
  int rev_e_id = G->GetEId(v,u);
  int delta = arcs[e_id].capacity - arcs[e_id].flow;
  int ex = nodes[u].excess;
  if (ex > delta) { ex = delta; }
  arcs[e_id].flow += ex;
  arcs[rev_e_id].flow -= ex;
  nodes[u].excess -= ex;
  nodes[v].excess += ex;
}


static inline void relabel(int u) {
  work_since_update += beta;
  TNEANet::TNodeI NI = G->GetNI(u);
  nodes[u].height = num_nodes;
  int min_neighbor = num_nodes;
  for (int i = 0; i < NI.GetOutDeg(); ++i) {
    work_since_update++;
    int v = NI.GetOutNId(i);
    if (min_neighbor > nodes[v].height) {
      int e_id = G->GetEId(u,v);
      if (arcs[e_id].capacity - arcs[e_id].flow > 0) {
        min_neighbor = nodes[v].height;
      }
    }
  }
  if (min_neighbor < num_nodes - 1) {
    nodes[u].height = min_neighbor + 1;
    if (max_distance < nodes[u].height) {
      max_distance = nodes[u].height;
    }
  }
}


static inline void gap(int missing_height) {
  for (int h = missing_height + 1; h < max_distance; ++h) {
    std::list<int>::iterator li;
    for (li = layers[h].inactive.begin(); li != layers[h].inactive.end(); ++li) {
      nodes[*li].height = num_nodes;
    }
    layers[h].inactive.clear();
  }
  max_distance = missing_height - 1;
  max_active = max_distance;

}


static inline void discharge(int u) {
  TNEANet::TNodeI NI = G->GetNI(u);
  while (1) {
    for (int i = 0; i < NI.GetOutDeg(); ++i) {
      int v = NI.GetOutNId(i);
      if (nodes[u].height == nodes[v].height + 1) {
        int e_id = G->GetEId(u,v);
        if (arcs[e_id].capacity - arcs[e_id].flow > 0) {
          if (v != t && nodes[v].excess == 0) {
            remove_from_inactive_layer(v);
            add_to_active_layer(v);
          }
          push(u, v, e_id);
          if (nodes[u].excess == 0) {
            break;
          }
        }
      }
    }

    if (nodes[u].excess > 0) {
      int old_height = nodes[u].height;
      relabel(u);
      if (layers[old_height].active.empty() && layers[old_height].inactive.empty()) {
        gap(old_height);
      }
      if (nodes[u].height == num_nodes) { //IMPORTANT
        break;
      }
    } else {
      add_to_inactive_layer(u);
      break;
    }
  }
}


int push_relabel() {
  int u, v;


  for (int i = 0; i < num_nodes; ++i) { nodes[i].height = 1; }
  nodes[s].height = num_nodes;
  nodes[t].height = 0;

  max_active = 0;
  min_active = num_nodes;
  max_distance = num_nodes - 1;
  nodes[s].excess = INF;

  TNEANet::TNodeI NI = G->GetNI(s);
  for (int i = 0; i < NI.GetOutDeg(); ++i) {
    v = NI.GetOutNId(i);
    int e_id = G->GetEId(s, v);
    push(s, v, e_id);
  }
  nodes[s].excess = 0;

  for (int i = 0; i < num_nodes; ++i) {
    if (nodes[i].excess > 0) {
      add_to_active_layer(i);
    } else if (nodes[i].height < num_nodes) {
      add_to_inactive_layer(i);
    }
  }
  global_relabel();

  while (max_active >= min_active) {
    std::list<int>::iterator next = layers[max_active].active.begin();
    if (next == layers[max_active].active.end()) {
      max_active--;
    } else {
      u = *next;
      remove_from_active_layer(u);
      discharge(u);
      if (work_since_update*global_update_freq > normalized_size) {
        work_since_update = 0;
        global_relabel();
      }
    }
  }

  return nodes[t].excess;
}


void usage() { printf("USAGE: pushrelabel filename\n"); }


void reset() {
  for (int i = 0; i < num_nodes; ++i) {
    layers[i].active.clear();
    layers[i].inactive.clear();
    nodes[i].height = 0;
    nodes[i].excess = 0;
  }
  for (int j = 0; j < num_edges; ++j) {
    arcs[j].flow = 0;
  }
  work_since_update = 0;
}


int main(int argc, char* argv[]) {
  if (argc <= 1) { usage(); return 0; }
  char *filename = argv[1];
  G = PNEANet::New();
  net_from_dimacs(filename);
  normalized_size = alpha*num_nodes + num_edges/2;

  nodes = (Node *) calloc(num_nodes, sizeof(Node));
  layers = new Layer[num_nodes];

  s = 0;
  t = num_nodes - 1;

  clock_t start = clock();
  int flow = push_relabel();
  clock_t end = clock();
  printf("%f\t", ((float)end - start)/CLOCKS_PER_SEC);
  /*for (int i = 1; i <= 8; ++i) {
    global_update_freq = (float)i/4;
    printf("%f\t", global_update_freq);
  }
  printf("\n");*/
  /*
  if (num_nodes > 150000) { return 0; }
  float updates[6] = {0.25, 0.5, 1.0, 2.0, 4.0, 8.0};
  for (int i = 0; i < 6; ++i) {
    global_update_freq = updates[5];
    reset();
  }*/
  //printf("%d\t%d\t", num_nodes, num_edges);
  //printf("%d\t", flow);
  //printf("%f\t", ((float)end - start)/CLOCKS_PER_SEC);
  fflush(stdout);

  return 0;
}
