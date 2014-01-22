

#make NAME=biparite all
g++ -O3 make_bipartite.cpp -o bipartite

./bipartite 100000 50 1 > BiTests/bi_100000_50_1.txt
./bipartite 1000000 50 1 > BiTests/bi_1000000_50_1.txt
