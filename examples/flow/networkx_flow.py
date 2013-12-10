

import networkx as nx
import sys

def graph_from_dimacs(filename):
    f = open(filename)
    g = nx.Graph()
    for line in f:
        if line[0] == 'a':
            parsed = line.split()
            a = int(parsed[1])
            b = int(parsed[2])
            c = int(parsed[3])
            g.add_edge(a,b,capacity=c)
    f.close()
    return g

for arg in sys.argv:
    if arg[-4:] == ".txt":
        g = graph_from_dimacs(arg)
        print nx.min_cut(g, 1, len(g.nodes()))

