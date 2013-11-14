// pushrelabelmain.cpp : Defines the entry point for the console application.

#include "stdafx.h"

#define INF 1073741824	
#define MIN(A,B) ((A < B) ? (A) : (B))
#define GLOBAL_UPDATE_FREQ 0.5

void net_from_dimacs(PNEANet G, char *filename) {
  FILE *f = fopen(filename, "r");
  if (f == NULL) { printf("Couldn't find file!\n"); return; }
  int a, b, c, e_id;
  char z, buf[1024];
  while (fgets(buf, 1024, f) != NULL) {
    if (buf[0] == 'a') {
      sscanf(buf, "%c %d %d %d", &z, &a, &b, &c);
      if (!G->IsNode(a)) { G->AddNode(a); }
      if (!G->IsNode(b)) { G->AddNode(b); }
      e_id = G->AddEdge(a,b);
      G->AddIntAttrDatE(e_id, c, "capacity");
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


/*int push_relabel(PNEANet G, int s, int t) { return 0; }


int push_dfs_flow(PNEANet G, int s, int t, int *dfs_back) {
  int cur_node = t;
  int min_flow = INF;
  while (cur_node != s) {
    int prev_node = dfs_back[cur_node];

    int e_id = G->GetEId(prev_node, cur_node);
    int rev_e_id = G->GetEId(cur_node, prev_node);
    int exc = G->GetIntAttrDatE(e_id, "capacity") - G->GetIntAttrDatE(e_id, "flow");
    exc += G->GetIntAttrDatE(rev_e_id, "flow");

    if (exc < min_flow) { min_flow = exc; }
    cur_node = prev_node;
  }
  cur_node = t;
  while (cur_node != s) {
    int prev_node = dfs_back[cur_node];
    int e_id = G->GetEId(prev_node, cur_node);
    int rev_e_id = G->GetEId(cur_node, prev_node);

    int flow_dif = min_flow;
    int back_flow = G->GetIntAttrDatE(rev_e_id, "flow");
    int cur_flow = G->GetIntAttrDatE(e_id, "flow");
    if (back_flow <= flow_dif) {
      G->AddIntAttrDatE(rev_e_id, 0, "flow");
      G->AddIntAttrDatE(e_id, cur_flow + (flow_dif - back_flow), "flow");
    } else if (back_flow > flow_dif) {
      G->AddIntAttrDatE(rev_e_id, back_flow - flow_dif, "flow");
    }
    cur_node = prev_node;
  }
  return min_flow;
}


int residual_dfs(PNEANet G, int s, int t, int *dfs_back) {
  static int iter_id = 0;
  iter_id++;

  TSnapQueue<int> node_queue(G->GetNodes());
  node_queue.Push(s);

  while (node_queue.Empty() == 0) {
    int cur_node = node_queue.Top();
    node_queue.Pop();
    TNEANet::TNodeI NI = G->GetNI(cur_node);

    for (int i = 0; i < NI.GetOutDeg(); ++i) {
      int n_id = NI.GetOutNId(i);
      int e_id = G->GetEId(cur_node, n_id);
      int rev_e_id = G->GetEId(n_id, cur_node);
      int exc = G->GetIntAttrDatE(e_id, "capacity") - G->GetIntAttrDatE(e_id, "flow");
      exc += G->GetIntAttrDatE(rev_e_id, "flow");

      if (exc > 0) {
        int visited = G->GetIntAttrDatN(n_id, "visited");

        if (visited != iter_id) {
            dfs_back[n_id] = cur_node;
            if (n_id == t) { return push_dfs_flow(G, s, t, dfs_back); }
            G->AddIntAttrDatN(n_id, iter_id, "visited");
            node_queue.Push(n_id);
        }

      }
    }
  }
  return 0;
}



int ford_fulkerson(PNEANet G, int s, int t) {
  if (s == t) { return INF; }
  // TODO: Change to hash map based on node id
  int *dfs_back = (int *) malloc(sizeof(int) * (G->GetNodes() + 1));
  init_flow(G);

  for (TNEANet::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
    G->AddIntAttrDatN(NI.GetId(), 0, "visited");
  }

  int flow = 0, old_flow = -1;
  while (flow != old_flow) {
    old_flow = flow;
    flow += residual_dfs(G, s, t, dfs_back);
    //printf("ITER_FLOW:%d\n", flow);
  }

  return flow;
}*/


void push(PNEANet G, int u, int v, int *e) {
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


void global_relabel_heuristic(PNEANet G, int t, int *h) {
  //TODO: Change this to avoid the fact that we are using a 1 indexed node
  for (int i = 1; i <= G->GetNodes(); ++i) {
    h[i] = 1;
  }
  TNEANet::TNodeI NI = G->GetNI(t);
  TSnapQueue<int> node_queue(G->GetNodes());
  node_queue.Push(t);
  while (node_queue.Empty() == 0) {
    int cur_node = node_queue.Top();
    node_queue.Pop();
    NI = G->GetNI(cur_node);
    for (int i = 0; i < NI.GetInDeg(); ++i) {
      int prev_node = NI.GetInNId(i);
      if (h[prev_node] == 0) {
        int e_id = G->GetEId(prev_node, cur_node);
        if (G->GetIntAttrDatE(e_id, "capacity") - G->GetIntAttrDatE(e_id, "flow") > 0) {
          h[prev_node] = h[cur_node] + 1;
          node_queue.Push(prev_node);
        }
      }
    }
  }
}


void relabel(PNEANet G, int u, int t, int *h) {
  static int counter = 0;
  counter++;
  if (counter > GLOBAL_UPDATE_FREQ*G->GetNodes()) {
    counter = 0;
    global_relabel_heuristic(G, t, h);
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

int push_relabel_global_heuristic(PNEANet G, int s, int t) {
  // Init
  init_flow(G);
  int min_height, u, v, n = G->GetNodes();
  int *height = (int *) calloc(n+1, sizeof(int));
  int *excess = (int *) calloc(n+1, sizeof(int));
  int *in_queue = (int *) calloc(n+1, sizeof(int));
  height[s] = n;
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


int push_relabel(PNEANet G, int s, int t) {
  // Init
  init_flow(G);
  int min_height, u, v, n = G->GetNodes();
  int *height = (int *) calloc(n+1, sizeof(int));
  int *excess = (int *) calloc(n+1, sizeof(int));
  int *in_queue = (int *) calloc(n+1, sizeof(int));
  height[s] = n;
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
        } else {
          min_height = MIN(min_height, height[v]);
        }
      }
    }

    if (excess[u] == 0) {
      in_queue[u] = 0;
      node_queue.Pop();
    } else if (excess[u] > 0) {
      height[u] = min_height + 1;//relabel(G, u, h);
    }
  }

  return excess[t];
}

void usage() { printf("USAGE: pushrelabel filename\n"); }


int main(int argc, char* argv[]) {
  if (argc <= 1) { usage(); return 0; }
  PNEANet G = PNEANet::New();
  char *f_name = argv[1];
  net_from_dimacs(G, f_name);

  printf("%d\n", push_relabel(G, 1, G->GetNodes())); //TODO: Parse s,t

  return 0;
}
