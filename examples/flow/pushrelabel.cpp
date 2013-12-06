// pushrelabelmain.cpp : Defines the entry point for the console application.

#include "stdafx.h"
#include <ctime>

#define INF (1 << 30)
#define MIN(A,B) ((A < B) ? (A) : (B))
#define GLOBAL_UPDATE_FREQ 1.0

void net_from_dimacs(PNEANet G, char *filename, int **capacities, int **flows) {
  FILE *f = fopen(filename, "r");
  if (f == NULL) { printf("Couldn't find file!\n"); return; }
  int a, b, c, e_id, max_edges, max_nodes;
  char z, buf[1024];
  while (fgets(buf, 1024, f) != NULL) {
    if (buf[0] == 'p') {
      sscanf(buf, "p max %d %d", &max_nodes, &max_edges);
      *capacities = (int *) calloc(2*max_edges, sizeof(int));
      *flows = (int *) calloc(2*max_edges, sizeof(int));
    }
    if (buf[0] == 'a') {
      sscanf(buf, "%c %d %d %d", &z, &a, &b, &c);
      if (!G->IsNode(a-1)) { G->AddNode(a-1); }
      if (!G->IsNode(b-1)) { G->AddNode(b-1); }
      if (!G->IsEdge(a-1,b-1)) { G->AddEdge(a-1,b-1); }
      //if (!G->IsEdge(b-1,a-1)) { G->AddEdge(b-1,a-1); }
      e_id = G->GetEId(a-1,b-1);
      G->AddIntAttrDatE(e_id, c, "capacity");
      (*capacities)[e_id] = c;
      //e_id = G->GetEId(b-1,a-1);
      //G->AddIntAttrDatE(e_id, 0, "capacity");
    }
  }
  // Make Undirected
  for (TNEANet::TEdgeI EI = G->BegEI(); EI < G->EndEI(); EI++) {
    int s_id = EI.GetSrcNId();
    int d_id = EI.GetDstNId();
    if (! G->IsEdge(d_id, s_id)) {
      int e_id = G->AddEdge(d_id, s_id);
      G->AddIntAttrDatE(e_id, 0, "capacity");
    }
  }
  fclose(f);
}


void init_flow(PNEANet G) {
  for (TNEANet::TEdgeI EI = G->BegEI(); EI < G->EndEI(); EI++) {
    G->AddIntAttrDatE(EI.GetId(), 0, "flow");
  }
}


static inline void push(PNEANet G, int u, int v, int *e) {
  int e_id = G->GetEId(u,v);
  int rev_e_id = G->GetEId(v,u);
  int f = G->GetIntAttrDatE(e_id, "flow");

  int delta = G->GetIntAttrDatE(e_id, "capacity") - f;
  int amt = MIN( e[u], delta );
  G->AddIntAttrDatE(e_id, f + amt, "flow");
  G->AddIntAttrDatE(rev_e_id, -(f + amt), "flow");
  e[u] -= amt;
  e[v] += amt;
}


static inline void global_relabel(PNEANet G, int t, int *h) {
  //TODO: Change this to avoid the fact that we are using a 1 indexed node
  //TODO: can be done without linear pass, also use a constant buffer
  for (int i = 0; i < G->GetNodes(); ++i) {
    h[i] = 10000000;
  }
  h[t] = 0;
  TNEANet::TNodeI NI = G->GetNI(t);
  TSnapQueue<int> node_queue(G->GetNodes());
  node_queue.Push(t);
  while (node_queue.Empty() == 0) {
    int cur_node = node_queue.Top();
    node_queue.Pop();
    NI = G->GetNI(cur_node);
    for (int i = 0; i < NI.GetInDeg(); ++i) {
      int prev_node = NI.GetInNId(i);
      if (h[prev_node] > h[cur_node] + 1) {
        int e_id = G->GetEId(prev_node, cur_node);
        if (G->GetIntAttrDatE(e_id, "capacity") - G->GetIntAttrDatE(e_id, "flow") > 0) {
          h[prev_node] = h[cur_node] + 1;
          node_queue.Push(prev_node);
        }
      }
    }
  }
}


static inline void relabel(PNEANet G, int u, int t, int *h) {
  static int counter = 0;
  counter++;
  if (counter > GLOBAL_UPDATE_FREQ*G->GetNodes()) {
    counter = 0;
    global_relabel(G, t, h);
  }
  TNEANet::TNodeI NI = G->GetNI(u);
  int min_neighbor = INF;
  for (int i = 0; i < NI.GetOutDeg(); ++i) {
    int v = NI.GetOutNId(i);
    int e_id = G->GetEId(u,v);
    if (G->GetIntAttrDatE(e_id, "capacity") - G->GetIntAttrDatE(e_id, "flow") > 0 ) {
      if (min_neighbor > h[v]) { min_neighbor = h[v]; }
    }
  }
  h[u] = min_neighbor + 1;
}


int push_relabel(PNEANet G, int s, int t, int *capacities, int *flows) {
  // Init
  init_flow(G);
  int min_height, u, v, n = G->GetNodes();
  int *height = (int *) calloc(n+1, sizeof(int));
  int *excess = (int *) calloc(n+1, sizeof(int));
  bool *in_queue = (bool *) calloc(n+1, sizeof(bool));
  global_relabel(G, t, height);
  //height[s] = n;
  excess[s] = INF;
  TNEANet::TNodeI NI = G->GetNI(s);
  TSnapQueue<int> node_queue(n);
  for (int i = 0; i < NI.GetOutDeg(); ++i) {
    v = NI.GetOutNId(i);
    push(G, s, v, excess);
    node_queue.Push(v);
    in_queue[v] = 1;
  }

  while (node_queue.Empty() == 0) {
    u = node_queue.Top();
    min_height = INF;
    NI = G->GetNI(u);
    for (int i = 0; i < NI.GetOutDeg(); ++i) {
      v = NI.GetOutNId(i);
      int e_id = G->GetEId(u,v);
      if (G->GetIntAttrDatE(e_id, "capacity") - G->GetIntAttrDatE(e_id, "flow") > 0 ) {
        if (height[u] > height[v]) {
          push(G, u, v, excess);
          if (in_queue[v] == 0 && v != s && v != t) {
            in_queue[v] = 1;
            node_queue.Push(v);
          }
        }
      }
    }

    if (excess[u] == 0) {
      in_queue[u] = 0;
      node_queue.Pop();
    } else if (excess[u] > 0) {
      relabel(G, u, t, height);
    }
  }

  return excess[t];
}


void usage() { printf("USAGE: pushrelabel filename\n"); }


int main(int argc, char* argv[]) {
  if (argc <= 1) { usage(); return 0; }
  char *filename = argv[1];
  PNEANet G = PNEANet::New();

  int *capacities, *flows;

  clock_t start_init = clock();
  net_from_dimacs(G, filename, &capacities, &flows);
  clock_t end_init = clock();

  printf("%s\t", filename);
  printf("Init:%f\t", ((float)end_init - start_init)/CLOCKS_PER_SEC);
  fflush(stdout);
  clock_t start = clock();
  int flow = push_relabel(G, 0, G->GetNodes()-1, capacities, flows);
  clock_t end = clock();
  printf("Flow:%d\t", flow);
  printf("Time:%f\n", ((float)end - start)/CLOCKS_PER_SEC);

  return 0;
}
