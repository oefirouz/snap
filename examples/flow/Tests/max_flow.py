

import networkx as nx

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

g = graph_from_dimacs('rmflong_16_4_1024_4608.txt')
print nx.min_cut(g,1,len(g.nodes())/2)
