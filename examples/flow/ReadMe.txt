========================================================================
    Max-Flow using Push-Relabel and Heuristics 	 
========================================================================

The example implements overlapping community detection by the 
Community-Affiliation Graph Model (AGM).



Fitting procedure and the community-Affiliation Graph Model are 
described in the following paper:
J. Yang and J. Leskovec, Structure and Overlaps of Communities in 
Networks, SNAKDD '12.
J. Yang and J. Leskovec, Community-Affiliation Graph Model for 
Overlapping Community Detection, ICDM '12.

The code works under Windows with Visual Studio or Cygwin with GCC,
Mac OS X, Linux and other Unix variants with GCC. Make sure that a
C++ compiler is installed on the system. Visual Studio project files
and makefiles are provided. For makefiles, compile the code with
"make all".

////////////////////////////////////////////////////////////////////////
Parameters:
   -o: Output file name prefix (default:'')
   -i: Input edgelist file name. DEMO: AGM with 2 communities
   -l:  Input file name for node names (Node ID, Node label) 
   -s: Random seed for AGM
   -e: Edge probability between the nodes that do not share any 
     community: set it to be 1 / N^2)
   -c: Number of communities (0: determine it by AGM)

////////////////////////////////////////////////////////////////////////
Usage:

Detect 12 communities of universities (which correspond to NCAA 
conferences) from the network of NCAA football teams:

agmfitmain -i:football.edgelist -l:football.labels -c:12 -e:0.1
