/*
Info:
Feel free to use these lines as you wish.
This program iterates over all k-cliques.
This is an improvement of the 1985 algorithm of Chiba And Nishizeki detailed in "Arboricity and subgraph listing".

To compile:
"gcc kClistNodeParallel.c -O9 -o kClistNodeParallel -fopenmp".

To execute:
"./kClistNodeParallel p k edgelist.txt".
"edgelist.txt" should contain the graph: one edge on each line separated by a space.
k is the size of the k-cliques
p is the number of threads
Will print the number of k-cliques.
*/


#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <omp.h>


#define NLINKS 100000000 //maximum number of edges for memory allocation, will increase if needed

typedef struct {
	unsigned s;
	unsigned t;
} edge;


typedef struct {
	unsigned n;//number of nodes
	unsigned e;//number of edges
	edge *edges;//list of edges

} edgelist;

typedef struct {
	unsigned n;
	unsigned *cd;//cumulative degree: (starts with 0) length=n+1
	unsigned *adj;//truncated list of neighbors
	unsigned core;//core value of the graph
} graph;

typedef struct {
	unsigned *n;//n[l]: number of nodes in G_l
	unsigned **d;//d[l]: degrees of G_l
	unsigned *adj;//truncated list of neighbors
	unsigned char *lab;//lab[i] label of node i
	unsigned **nodes;//sub[l]: nodes in G_l
	unsigned core;
} subgraph;

unsigned *d0;

void free_edgelist(edgelist *el) {
	free(el->edges);
	free(el);
}

void free_graph(graph *g) {
	free(g->cd);
	free(g->adj);
	free(g);
}

void free_subgraph(subgraph *sg, unsigned char k) {
	unsigned char i;
	free(sg->n);
	for (i = 2; i<k; i++) {
		free(sg->d[i]);
		free(sg->nodes[i]);
	}
	free(sg->d);
	free(sg->nodes);
	free(sg->lab);
	free(sg->adj);
	free(sg);
}


//Compute the maximum of three unsigned integers.
inline unsigned int max3(unsigned int a, unsigned int b, unsigned int c) {
	a = (a>b) ? a : b;
	return (a>c) ? a : c;
}

edgelist* readedgelist(char* input) {
	unsigned e1 = NLINKS;
	edgelist *el = malloc(sizeof(edgelist));
	FILE *file;

	el->n = 0;
	el->e = 0;
	file = fopen(input, "r");
	el->edges = malloc(e1 * sizeof(edge));
	while (fscanf(file, "%u %u", &(el->edges[el->e].s), &(el->edges[el->e].t)) == 2) {//Add one edge
		el->n = max3(el->n, el->edges[el->e].s, el->edges[el->e].t);
		el->e++;
		if (el->e == e1) {
			e1 += NLINKS;
			el->edges = realloc(el->edges, e1 * sizeof(edge));
		}
	}
	fclose(file);
	el->n++;

	el->edges = realloc(el->edges, el->e * sizeof(edge));

	return el;
}

void relabel(edgelist *el) {
	unsigned i, dsource, dtarget, tmp;

	for (i = 0; i < el->e; i++) {
		dsource = d0[el->edges[i].s];
		dtarget = d0[el->edges[i].t];
		if ( (dsource > dtarget ) || (dsource == dtarget && el->edges[i].s > el->edges[i].t)) {
			tmp = el->edges[i].s;
			el->edges[i].s = el->edges[i].t;
			el->edges[i].t = tmp;
		}
	}
}


//computing degree ordering
void ord_degree(edgelist* el) {

	unsigned i;
	d0 = calloc(el->n, sizeof(unsigned));
	for (i = 0; i < el->e; i++) {
		d0[el->edges[i].s]++;
		d0[el->edges[i].t]++;
	}

}

//////////////////////////
//Building the special graph
graph* mkgraph(edgelist *el) {
	unsigned i, max;
	unsigned *d;
	graph* g = malloc(sizeof(graph));

	d = calloc(el->n, sizeof(unsigned));

	for (i = 0; i<el->e; i++) {
		d[el->edges[i].s]++;
	}

	g->cd = malloc((el->n + 1) * sizeof(unsigned));
	g->cd[0] = 0;
	max = 0;
	for (i = 1; i<el->n + 1; i++) {
		g->cd[i] = g->cd[i - 1] + d[i - 1];
		max = (max>d[i - 1]) ? max : d[i - 1];
		d[i - 1] = 0;
	}
	printf("core value (max truncated degree) = %u\n", max);

	g->adj = malloc(el->e * sizeof(unsigned));

	for (i = 0; i<el->e; i++) {
		g->adj[g->cd[el->edges[i].s] + d[el->edges[i].s]++] = el->edges[i].t;
	}

	free(d);
	g->core = max;
	g->n = el->n;
	return g;
}


subgraph* allocsub(graph *g, unsigned char k) {
	unsigned i;
	subgraph* sg = malloc(sizeof(subgraph));
	sg->n = calloc(k, sizeof(unsigned));
	sg->d = malloc(k * sizeof(unsigned*));
	sg->nodes = malloc(k * sizeof(unsigned*));
	for (i = 2; i<k; i++) {
		sg->d[i] = malloc(g->core * sizeof(unsigned));
		sg->nodes[i] = malloc(g->core * sizeof(unsigned));
	}
	sg->lab = calloc(g->core, sizeof(unsigned char));
	sg->adj = malloc(g->core*g->core * sizeof(unsigned));
	sg->core = g->core;
	return sg;
}

void mksub(graph* g, unsigned u, subgraph* sg, unsigned char k) {
	unsigned i, j, l, v, w;

	static unsigned *old = NULL, *new = NULL;//to improve
#pragma omp threadprivate(new,old)

	if (old == NULL) {
		new = malloc(g->n * sizeof(unsigned));
		old = malloc(g->core * sizeof(unsigned));
		for (i = 0; i<g->n; i++) {
			new[i] = -1;
		}
	}

	for (i = 0; i<sg->n[k - 1]; i++) {
		sg->lab[i] = 0;
	}

	j = 0;
	for (i = g->cd[u]; i<g->cd[u + 1]; i++) {
		v = g->adj[i];
		new[v] = j;
		old[j] = v;
		sg->lab[j] = k - 1;
		sg->nodes[k - 1][j] = j;
		sg->d[k - 1][j] = 0;//new degrees
		j++;
	}

	sg->n[k - 1] = j;

	for (i = 0; i<sg->n[k - 1]; i++) {//reodering adjacency list and computing new degrees
		v = old[i];
		for (l = g->cd[v]; l<g->cd[v + 1]; l++) {
			w = g->adj[l];
			j = new[w];
			if (j != -1) {
				sg->adj[sg->core*i + sg->d[k - 1][i]++] = j;
			}
		}
	}

	for (i = g->cd[u]; i<g->cd[u + 1]; i++) {
		new[g->adj[i]] = -1;
	}
}

void kclique_thread(unsigned char l, subgraph *sg, unsigned long long *n) {
	unsigned i, j, k, end, u, v, w;

	if (l == 2) {
		for (i = 0; i<sg->n[2]; i++) {//list all edges
			u = sg->nodes[2][i];
			(*n) += sg->d[2][u];
			/*
			end = u*sg->core + sg->d[2][u];
			for (j = u*sg->core; j<end; j++) {
				(*n)++;//listing here!!!  // NOTE THAT WE COULD DO (*n)+=g->d[2][u] to be much faster (for counting only); !!!!!!!!!!!!!!!!!!
			}
			*/
		}
		return;
	}

	for (i = 0; i<sg->n[l]; i++) {
		u = sg->nodes[l][i];
		//printf("%u %u\n",i,u);
		sg->n[l - 1] = 0;
		end = u*sg->core + sg->d[l][u];
		for (j = u*sg->core; j<end; j++) {//relabeling nodes and forming U'.
			v = sg->adj[j];
			if (sg->lab[v] == l) {
				sg->lab[v] = l - 1;
				sg->nodes[l - 1][sg->n[l - 1]++] = v;
				sg->d[l - 1][v] = 0;//new degrees
			}
		}
		for (j = 0; j<sg->n[l - 1]; j++) {//reodering adjacency list and computing new degrees
			v = sg->nodes[l - 1][j];
			end = sg->core*v + sg->d[l][v];
			for (k = sg->core*v; k<end; k++) {
				w = sg->adj[k];
				if (sg->lab[w] == l - 1) {
					sg->d[l - 1][v]++;
				}
				else {
					sg->adj[k--] = sg->adj[--end];
					sg->adj[end] = w;
				}
			}
		}

		kclique_thread(l - 1, sg, n);

		for (j = 0; j<sg->n[l - 1]; j++) {//restoring labels
			v = sg->nodes[l - 1][j];
			sg->lab[v] = l;
		}

	}
}

unsigned long long kclique_main(unsigned char k, graph *g) {
	int u;
	unsigned long long n = 0;
	subgraph *sg;
#pragma omp parallel private(sg,u) reduction(+:n)
	{
		sg = allocsub(g, k);
#pragma omp for schedule(dynamic, 1) nowait
		for (u = 0; u<g->n; u++) {
			mksub(g, u, sg, k);
			kclique_thread(k - 1, sg, &n);
		}
	}
	return n;
}

int main(int argc, char** argv) {
	edgelist* el;
	graph* g;
	unsigned char k = atoi(argv[2]);
	unsigned long long n;

	omp_set_num_threads(atoi(argv[1]));

	time_t t0, t1, t2;
	t1 = time(NULL);
	t0 = t1;

	printf("Reading edgelist from file %s\n", argv[3]);

	el = readedgelist(argv[3]);
	printf("Number of nodes = %u\n", el->n);
	printf("Number of edges = %u\n", el->e);

	t2 = time(NULL);
	printf("- Time = %ldh%ldm%lds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
	t1 = t2;

	printf("Building the graph structure\n");
	ord_degree(el);
	relabel(el);
	g = mkgraph(el);

	printf("Number of nodes (degree > 0) = %u\n", g->n);

	free_edgelist(el);
	
	t2 = time(NULL);
	printf("- Time = %ldh%ldm%lds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
	t1 = t2;

	printf("Iterate over all cliques\n");

	n = kclique_main(k, g);

	printf("Number of %u-cliques: %llu\n", k, n);

	t2 = time(NULL);
	printf("- Time = %ldh%ldm%lds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
	t1 = t2;

	free_graph(g);

	printf("- Overall time = %ldh%ldm%lds\n", (t2 - t0) / 3600, ((t2 - t0) % 3600) / 60, ((t2 - t0) % 60));

	return 0;
}
