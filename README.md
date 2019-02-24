# kCliqueListing
7 **k-cliques** algorithms
> Listing k-cliques in Sparse Real-World Graphs


Finding dense subgraphs is an important research area in **graph mining**, with applications in
  * community detection spamlink farms in web graphs;
  * real-time story identification,motif detection in biological networks;
  * epilepsy prediction,graph compression;
  * distance query indexing;
  * finding correlated genes; 
  * finance and many others.

## 7 **k-cliques** algorithms
- **Arboricity** ：Arboricity algorithm；
- **Degree** ：Degree-order algorithm；
- **Degen** ：Degeneracy-order algorithm；
- **DegCol** ：First descending according to degree, then greedy color order algorithm；
- **DegenCol** ：First according to the degeneracy reverse order, then the greedy color order algorithm；
- **DDegCol** ：First calculate out-neighbors according to degeneracy, then perform the order algorithm of greedy color according to degree descending in out-neighbor of each node；
- **DDegree** ：First calculate out-neighbors according to degeneracy, then follow the degree-order algorithm in out-neighbor of each node.；
