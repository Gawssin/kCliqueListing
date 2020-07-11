# Info
This version can output the largest core, largest truss, and topk maximum cliques to the corresponding file.

# To compile
`g++ "-std=c++11" -O3 -o run main.cpp core_and_truss_decomposition.cpp`

# To run
`./run filepath alg topsize`
* **run** executable file
* **filepath** input file path
* **alg** algorithm, "1" for core, "2" for truss, "3" for maximum cliques
* **topsize** output topsize maximum cliques, only work for "alg=3"
