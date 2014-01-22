
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>


int n, k, num_nodes, num_edges, s, t, num_left, *potential_mates, *num_mated; 


void usage() {
  printf("Usage: ./make_bipartite [n] [k] [seed]\n"
          "Makes a k regular bipartite graph with 2*n 'real' nodes, and a source/sink");
}


//The key is at any point in time only 0..num_left are valid entries of potential mates
//This is a variation on a Knuth Shuffle
void make_bipartite() {

  for (int i = 1; i <= n; ++i) {
    for (int j = 0; j < k; ++j) {
      //match to a potential mate
      float f = (float)rand()/(float)RAND_MAX;
      int mate = f*num_left;
      int mate_id = potential_mates[mate];
      printf("a\t%d\t%d\t1\n", i+1, mate_id+1);
      //update potential mates
      num_mated[mate_id] += 1;
      if (num_mated[mate_id] == k) {
        int temp_mate_id = potential_mates[num_left-1];
        potential_mates[mate] = temp_mate_id;
        num_left -= 1;
      }
    }
  }
  return;
}

int main(int argc, char *argv[]) {
  if (argc <= 3) { usage(); return 0; }
  n = atoi(argv[1]);
  k = atoi(argv[2]);
  assert(k <= n); assert(k > 0); assert(n > 0);

  int seed = atoi(argv[3]);
  srand(seed);
  num_left = n;
  potential_mates = (int *) calloc(n, sizeof(int));
  num_mated = (int *) calloc(2*(n+1), sizeof(int));
  for (int i = 0; i < n; ++i) {
    potential_mates[i] = 1 + n + i;
  }

  printf("c This file was generated by make_bipartite\n"
         "c Makes a k-regular bipartite graph with 2*n nodes\n"
         "c And a source/sink node for solving of the matching.\n"
         "c The parameters are: n: %d k: %d seed: %d\n", n, k, seed);

  num_nodes = 2*(n+1);
  num_edges = 2*n + n*k;
  s = 0;
  t = num_nodes - 1;
  printf("p max\t%d\t%d\n", num_nodes, num_edges);
  printf("n\t%d s\n", s+1);
  printf("n\t%d t\n", t+1);
  for (int i = 1; i <= n; ++i) {
    printf("a\t%d\t%d\t1\n", s+1, i+1);
  }

  make_bipartite();

  for (int i = 1 + n; i <= 2*n; ++i) {
    printf("a\t%d\t%d\t1\n", i+1, t+1);
  }

  return 0;
}