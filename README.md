# OvMOSP

OvMOSP
Python implementation of OvMOSP, an exact branch-and-bound algorithm for the
Shortest Path over the Efficient Set of the Multi-Objective Shortest Path problem (SPoMOSP):
among the efficient (Pareto-optimal) paths from s to t, find the one that minimizes the
principal objective d (the distance), without enumerating the whole Pareto front.

Requirements
Python 3.7 or later
IBM ILOG CPLEX with its Python API (commercial; free academic licence). The Benson efficiency
test is solved with CPLEX.
psutil (optional, only used to measure memory): pip install psutil

File	Purpose
OvMOSP_random_lambda.py	OvMOSP, Phase 2 with random weight vectors (fixed seed)
OvMOSP_nv_lamda.py	OvMOSP, Phase 2 with the Das-Dennis weight design
extract_5_variants_10_instances.py	Builds the New York subnetworks and the source-target pairs from the DIMACS files
truncate_criteria.py	Derives the instances with fewer secondary criteria from the instances with four
The two OvMOSP scripts share the same code; they differ only in how the weight vectors of
Phase 2 (initialization by scalarization) are generated. Phase 2 only provides an initial upper
bound: the algorithm is exact for any number of weight vectors, including none.

Instances format
_________________________________________
n m p s t
u v d c_1 ... c_p        (one line per arc)
_________________________________________
n vertices, m arcs, p secondary criteria, source s, target t. d is the principal
objective (strictly positive); `c_1 ... c_p` define efficiency (non-negative). Vertices are
numbered from 0 or from 1 (a numbering starting at 1 is converted automatically).


Running OvMOSP
Random weight vectors:
python OvMOSP_random_lambda.py instance.txt [alpha] [seed] [bet_timelimit]
alpha: number of random weight vectors (default 50; 0 disables Phase 2)
seed: random seed (default 12345)
bet_timelimit: time limit in seconds of each Benson test (default 6000)
Das-Dennis weight vectors:
python OvMOSP_nv_lamda.py instance.txt [H] [bet_timelimit]
H: resolution of the design, giving alpha = C(H-1, p-1) weight vectors (default 10; 0 disables Phase 2)
bet_timelimit: as above
Each criterion is divided by its minimum attainable value before the weighted sum is formed.
If a Benson test is not solved to proven optimality (for instance because of the time limit),
the program stops with an error instead of returning a result that is not guaranteed to be
exact; no RESULT_LINE is printed in that case.


Output 
The program prints a detailed report, and a final line RESULT_LINE: key=value;... with
alpha, UB_scalarisation (upper bound after Phase 2), UB_final (optimal distance),
cpu_scalarisation, cpu_total, profondeur_arbre (depth of the search tree), and the memory
fields peak_mem_mb (sampled resident memory), peak_mem_os_mb (peak reported by the
operating system), base_mem_mb (before the search), max_stock (maximum of the total number of
vectors stored in stock[v]), max_UBC and bet_cache. Times are wall-clock times in seconds.

Building the New York instances
The instances are built from the 9th DIMACS Implementation Challenge road network of New York
(distance and travel time) and from the DIMACS-Extended criteria (elevation difference, average
out-degree of adjacent vertices, hop count). All .gr files must list the same arcs in the same
order; the script stops otherwise.
python extract_5_variants_10_instances.py \
    --distance NY_distance.gr --time NY_travel_time.gr \
    --elevation NY_elevation.gr --avg_degree NY_avg_degree.gr \
    --hop_count NY_hop_count.gr --coords USA-road-d.NY.co \
    --radius_km 8 --nb_instances 10 --out_dir instances --seed 42
The principal objective `d` is the distance. The four secondary criteria are, in this order:
travel time, elevation difference, connectivity, hop count. For each borough (Manhattan,
Brooklyn, Queens, Bronx, Staten Island), all vertices within radius_km of a reference point
are kept, then the largest weakly connected component. Source-target pairs are drawn at random
(seed 42) among vertices separated by at least ten arcs (this bound is lowered if too few
candidates exist). Instances with fewer criteria are obtained by removing the last criteria:
python truncate_criteria.py --input_dir instances --out_dir_p3 instances_p3 --out_dir_p2 instances_p2


Licece
MIT (see LICENSE).
