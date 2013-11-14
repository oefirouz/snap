// pushrelabelmain.cpp : Defines the entry point for the console application.

#include <cassert>
#include "stdafx.h"
//#include "agm.h"
//#include "agmfit.h"
//
#define INF (1 << 30)


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
  fclose(f);
}



void init_flow(PNEANet G) {
  int m = G->GetEdges();
  for (int i = 0; i < m; ++i) {
    assert(G->IsEdge(m));
    G->AddIntAttrDatE(i, 0, "flow");
  }
}


int push_relabel(PNEANet G, int s, int t) { return 0; }


int push_dfs_flow(PNEANet G, int s, int t, int *dfs_back) {
  int cur_node = t;
  int min_flow = INF;
  while (cur_node != s) {
    int prev_node = dfs_back[cur_node];
    int e_id = G->GetEId(prev_node, cur_node);
    int exc = G->GetIntAttrDatE(e_id, "capacity") - G->GetIntAttrDatE(e_id, "flow");
    if (exc < min_flow) { min_flow = exc; }
    cur_node = prev_node;
  }
  cur_node = t;
  while (cur_node != s) {
    int prev_node = dfs_back[cur_node];
    int e_id = G->GetEId(prev_node, cur_node);
    int cur_flow = G->GetIntAttrDatE(e_id, "flow");
    G->AddIntAttrDatE(e_id, cur_flow + min_flow, "flow");
    cur_node = prev_node;
  }
  return min_flow;
}


int residual_dfs(PNEANet G, int s, int t, int *dfs_node, int *dfs_back, int *visited) {
  for (int i = 0; i <= G->GetNodes(); ++i) { visited[i] = 0; }
  int stack_pos = 0;
  dfs_node[0] = s;
  visited[s] = 1;

  while (stack_pos >= 0) {
    int cur_node = dfs_node[stack_pos--];
    TNEANet::TNodeI NI = G->GetNI(cur_node);
    for (int i = 0; i < NI.GetOutDeg(); ++i) {
      int n_id = NI.GetOutNId(i);
      int e_id = G->GetEId(cur_node, n_id);
      int exc = G->GetIntAttrDatE(e_id, "capacity") - G->GetIntAttrDatE(e_id, "flow");
      if (visited[n_id] == 0 && exc > 0) {
        dfs_back[n_id] = cur_node;
        if (n_id == t) { return push_dfs_flow(G, s, t, dfs_back); }
        visited[n_id] = 1;
        dfs_node[++stack_pos] = n_id;
      }
    }
  }
  return 0;
}



int ford_fulkerson(PNEANet G, int s, int t) {
  if (s == t) { return INF; }
  int *dfs_node = (int *) malloc(sizeof(int) * (G->GetNodes() + 1));
  int *dfs_back = (int *) malloc(sizeof(int) * (G->GetNodes() + 1));
  int *visited = (int *) malloc(sizeof(int) * (G->GetNodes() + 1));
  init_flow(G);
  int flow = 0;
  int old_flow = 1;
  while (flow != old_flow) {
    old_flow = flow;
    flow += residual_dfs(G, s, t, dfs_node, dfs_back, visited);
  }

  printf("%d\n", flow);
  free(dfs_node);
  free(dfs_back);
  free(visited);
  return 0;
}




int main(int argc, char* argv[]) {
  PNEANet G = PNEANet::New();
  char f_name[] = "Tests/rmflong_16_4_1024_4608.txt";
  net_from_dimacs(G, f_name);
  ford_fulkerson(G, 1, 512);
  //printf("%d\n", n);
  //g->Dump();




  //Env = TEnv(argc, argv, TNotify::StdNotify);
  //Env.PrepArgs(TStr::Fmt("cpm. build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
  //TExeTm ExeTm;
  //Try
  //TStr OutFPrx  = Env.GetIfArgPrefixStr("-o:", "", "Output file name prefix");
  //const TStr InFNm = Env.GetIfArgPrefixStr("-i:", "football.edgelist", "Input edgelist file name. DEMO: AGM with 2 communities");
  //const TStr LabelFNm = Env.GetIfArgPrefixStr("-l:", "football.labels", "Input file name for node names (Node ID, Node label) ");
  //const TInt RndSeed = Env.GetIfArgPrefixInt("-s:", 0, "Random seed for AGM");
  //const TFlt Epsilon = Env.GetIfArgPrefixFlt("-e:", 0, "Edge probability between the nodes that do not share any community (default (0.0): set it to be 1 / N^2)");
  //const TInt Coms = Env.GetIfArgPrefixInt("-c:", 0, "Number of communities (0: determine it by AGM)");

  //PUNGraph G = TUNGraph::New();
  //TVec<TIntV> CmtyVV;
  //TIntStrH NIDNameH;

  //if (InFNm == "DEMO") {
  //  TVec<TIntV> TrueCmtyVV;
  //  TRnd AGMRnd(RndSeed);

  //  //generate community bipartite affiliation
  //  const int ABegin = 0, AEnd = 70, BBegin = 30, BEnd = 100;
  //  TrueCmtyVV.Add(TIntV());
  //  TrueCmtyVV.Add(TIntV());
  //  for (int u = ABegin; u < AEnd; u++) {
  //    TrueCmtyVV[0].Add(u);
  //  }
  //  for (int u = BBegin; u < BEnd; u++) {
  //    TrueCmtyVV[1].Add(u);
  //  }
  //  G = TAGM::GenAGM(TrueCmtyVV, 0.0, 0.2, AGMRnd);
  //}
  //else if (LabelFNm.Len() > 0) {
  //  G = TSnap::LoadEdgeList<PUNGraph>(InFNm);
  //  TSsParser Ss(LabelFNm, ssfTabSep);
  //  while (Ss.Next()) {
  //    if (Ss.Len() > 0) { NIDNameH.AddDat(Ss.GetInt(0), Ss.GetFld(1)); }
  //  }
  //}
  //else {
  //  G = TAGMUtil::LoadEdgeListStr<PUNGraph>(InFNm, NIDNameH);
  //}
  //printf("Graph: %d Nodes %d Edges\n", G->GetNodes(), G->GetEdges());

  //int MaxIter = 50 * G->GetNodes() * G->GetNodes();
  //if (MaxIter < 0) { MaxIter = TInt::Mx; }
  //int NumComs = Coms;
  //if (NumComs < 2) {
  //  int InitComs;
  //  if (G->GetNodes() > 1000) {
  //    InitComs = G->GetNodes() / 5;
  //    NumComs = TAGMUtil::FindComsByAGM(G, InitComs, MaxIter, RndSeed, 1.5, Epsilon, OutFPrx);
  //  } else {
  //    InitComs = G->GetNodes() / 5;
  //    NumComs = TAGMUtil::FindComsByAGM(G, InitComs, MaxIter, RndSeed, 1.2, Epsilon, OutFPrx);
  //  }
  //}
  //TAGMFit AGMFit(G, NumComs, RndSeed);
  //if (Epsilon > 0) { AGMFit.SetPNoCom(Epsilon);  }
  //AGMFit.RunMCMC(MaxIter, 10);
  //AGMFit.GetCmtyVV(CmtyVV, 0.9999);

  //TAGMUtil::DumpCmtyVV(OutFPrx + "cmtyvv.txt", CmtyVV, NIDNameH);
  //TAGMUtil::SaveGephi(OutFPrx + "graph.gexf", G, CmtyVV, 1.5, 1.5, NIDNameH);
  //AGMFit.PrintSummary();

  //Catch

  //printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
  printf("Hello World!\n");

  return 0;
}
