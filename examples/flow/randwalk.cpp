// pushrelabelmain.cpp : Defines the entry point for the console application.

#include "stdafx.h"
#include <ctime>

#define INF (1 << 30)

static int num_nodes;
static int num_edges;
static int num_matched;
static int n;
static int k;

static int s;
static int t;
static PNEANet G;

static bool *is_matched;
static int *partner;
static int *prev_node;
static int *prev_walk_assigned;


void net_from_dimacs(char *filename) {
  FILE *f = fopen(filename, "r");
  if (f == NULL) { printf("Couldn't find file!\n"); return; }
  int a, b, c, max_edges, max_nodes;
  char buf[1024];
  //TODO: source/sink defined in DIMACS
  while (fgets(buf, 1024, f) != NULL) {
    if (buf[0] == 'p') {
      sscanf(buf, "p max %d %d", &max_nodes, &max_edges);
      num_nodes = max_nodes;
      num_edges = max_edges;
    } else if (buf[0] == 'a') {
      //TODO: assert check all capacities 1
      sscanf(buf, "a %d %d %d", &a, &b, &c);
      if (!G->IsNode(a-1)) { G->AddNode(a-1); }
      if (!G->IsNode(b-1)) { G->AddNode(b-1); }
      if (!G->IsEdge(a-1,b-1)) { G->AddEdge(a-1,b-1); }
    }
  }
  is_matched = (bool *) calloc(num_nodes,sizeof(bool));
  partner = (int *) calloc(num_nodes,sizeof(int));
  prev_node = (int *) calloc(num_nodes,sizeof(int));
  prev_walk_assigned = (int *) calloc(num_nodes,sizeof(int));
  n = (num_nodes-2)/2;
  k = (num_edges-2*n)/n;
  fclose(f);
}


static inline int rand_int(int max) {
  float f = (float)rand()/(float)RAND_MAX;
  return f*max;
}


static inline int sample_rand(int u) {
  TNEANet::TNodeI NI = G->GetNI(u);
  return NI.GetOutNId( rand_int( NI.GetOutDeg() ) );
}


static inline int rand_walk() {
  int b = 2*(2 + n/(n-num_matched));
  static int walk_id = 0;
  walk_id++;
  // pick an unmatched guy
  int guy = sample_rand(s);
  int first_guy = guy;
  prev_node[guy] = s;
  prev_walk_assigned[guy] = walk_id;
  while (1) {
    int girl = sample_rand(guy);
    if (girl != first_guy) {
      if (prev_walk_assigned[girl] != walk_id) {
        prev_node[girl] =  guy;
        prev_walk_assigned[girl] = walk_id;
      }
    }
    if (girl == t) { return 1; }
    if (--b <= 0) { return 0; }

    //This is equivalent to shrinking matches to supernodes
    if (is_matched[girl]) {
      guy = partner[girl];
      if (prev_walk_assigned[guy] != walk_id) {
        prev_walk_assigned[guy] = walk_id;
        prev_node[guy] = girl;
      }
      if (guy == first_guy) {
        printf("ERROR THIS SHOULDN'T HAPPEN\n"); fflush(stdout);
      }
    } else {
      guy = girl;
    }
  }
}


static inline void augment() {
  //printf("Starting augment\n"); fflush(stdout);
  int cur = prev_node[t];
  //printf("%d\n", cur);
  int prev = prev_node[cur];
  //printf("%d\n", prev);
  is_matched[cur] = 1; // the only matching that changes it the outside
  while (prev != s) {
    /*int flag = 0;
    for (int i = 1; i <= n; ++i) {
      if (prev_node[i] == s) {
        flag = 1;
      }
    }
    if (flag == 0) {
        printf("ERROR THIS SHOULDN'T HAPPEN NO ROUTE TO SOURCE\n"); fflush(stdout);

    }*/
    //printf("Cur %d Prev %d\n", cur, prev); fflush(stdout);
    //printf("IsPartner? %d\n", partner[prev] == cur);
    if (partner[prev] == cur) {
      // do nothing, will be fixed in a later step
    } else {
      partner[prev] = cur;
      partner[cur] = prev;
    }
    cur = prev;
    prev = prev_node[prev];
    //printf("Cur %d Prev %d\n", cur, prev); fflush(stdout);
    //printf("IsPartner? %d\n", partner[prev] == cur);
  }
  is_matched[cur] = 1;
  G->DelEdge(s, cur);
  //printf("Ending augment\n"); fflush(stdout);
}


void rand_walk_match() {
  num_matched = 0;
  while (num_matched < n) {
    if (rand_walk() == 1) {
      augment();
      num_matched++;
      //printf("Num Matched %d\n", num_matched);
    }
  }
  for (int i = 1; i <= n; ++i) {
    //printf("%d %d\n", i+1, partner[i]+1);
  }
}


void usage() { printf("USAGE: randwalk filename\n"); }


int main(int argc, char* argv[]) {
  srand(0);
  if (argc <= 1) { usage(); return 0; }
  char *filename = argv[1];
  G = PNEANet::New();
  net_from_dimacs(filename);

  s = 0;
  t = num_nodes - 1;

  clock_t start = clock();
  rand_walk_match();
  clock_t end = clock();
  printf("%f\t", ((float)end - start)/CLOCKS_PER_SEC);

  return 0;
}

