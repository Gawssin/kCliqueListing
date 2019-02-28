/*
Info:
Feel free to use these lines as you wish.
This program iterates over all k-cliques.


To compile:
"gcc DegColNodeParallel.c -O9 -o DegColNodeParallel -fopenmp".

To execute:
"./DegColNodeParallel p k edgelist.txt".
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
	unsigned id;
	unsigned degree;
} iddegree;

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
	unsigned e;
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
	unsigned *color;
	unsigned *cd;
} subgraph;


int *color;
unsigned *Index;

int cmp(const void* a, const void* b)
{
	iddegree *x = (iddegree*)a, *y = (iddegree*)b;
	return y->degree - x->degree;
}

int cmpadj(const void* a, const void* b)
{
	int *x = (int*)a, *y = (int*)b;
	return color[Index[*y]] - color[Index[*x]];
}

void free_edgelist(edgelist *el) {
	free(el->edges);
	free(el);
}

void free_graph(graph *g) {
	free(g->cd);
	free(g->adj);
	free(g);
}

void free_subgraph(subgraph *sg, unsigned char k)
{
	free(sg->n);
	free(sg->cd);
	free(sg->lab);
	free(sg->color);
	for (int i = 2; i < k; i++)
	{
		free(sg->d[i]);
		free(sg->nodes[i]);
	}
	free(sg->d);
	free(sg->nodes);
}


//Compute the maximum of three unsigned integers.
unsigned int max3(unsigned int a, unsigned int b, unsigned int c);
inline unsigned int max3(unsigned int a, unsigned int b, unsigned int c) {
	a = (a > b) ? a : b;
	return (a > c) ? a : c;
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


//computing degree ordering
void ord_color_relabel(edgelist* g) {
	unsigned i, r = 0, N = g->n;

	iddegree *ig = malloc(g->n * sizeof(iddegree));
	unsigned *d0 = calloc(g->n, sizeof(unsigned));
	unsigned *cd0 = malloc((g->n + 1) * sizeof(unsigned));
	unsigned *adj0 = malloc(2 * g->e * sizeof(unsigned));
	for (i = 0; i < g->e; i++) {
		d0[g->edges[i].s]++;
		d0[g->edges[i].t]++;
	}

	cd0[0] = 0;
	for (i = 1; i < g->n + 1; i++) {
		cd0[i] = cd0[i - 1] + d0[i - 1];
		ig[i - 1].id = i - 1;
		ig[i - 1].degree = d0[i - 1];
		d0[i - 1] = 0;
	}
	for (i = 0; i < g->e; i++) {
		adj0[cd0[g->edges[i].s] + d0[g->edges[i].s]++] = g->edges[i].t;
		adj0[cd0[g->edges[i].t] + d0[g->edges[i].t]++] = g->edges[i].s;
	}

	qsort(ig, N, sizeof(ig[0]), cmp);

	Index = malloc(N * sizeof(unsigned));
	for (int i = 0; i < N; i++)
		Index[ig[i].id] = i;


	color = malloc(N * sizeof(int));
	memset(color, -1, sizeof(int)*N);


	int *C = malloc((ig[0].degree + 1) * sizeof(int));
	memset(C, 0, sizeof(int)*(ig[0].degree + 1));
	color[0] = 0;
	int colorNum = 1;


	for (int i = 1; i < N; i++)
	{
		int tmpdegree = ig[i].degree, tmpid = ig[i].id;
		for (int j = 0; j < tmpdegree; j++)
		{
			int now = Index[adj0[cd0[tmpid] + j]];
			if (color[now] != -1)
				C[color[now]] = 1;
		}
		for (int j = 0; j < ig[0].degree + 1; j++)
			if (C[j] == 0)
			{
				color[i] = j;
				colorNum = j > colorNum ? j : colorNum;
				break;
			}

		for (int j = 0; j < tmpdegree; j++)
		{
			int now = Index[adj0[cd0[tmpid] + j]];
			if (color[now] != -1)
				C[color[now]] = 0;
		}

	}
	printf("color number = %d\n", colorNum);


	for (int i = 0; i < g->e; i++)
	{
		if (color[Index[g->edges[i].s]] < color[Index[g->edges[i].t]])
		{
			int tmp = g->edges[i].s;
			g->edges[i].s = g->edges[i].t;
			g->edges[i].t = tmp;
		}
		else if (color[Index[g->edges[i].s]] == color[Index[g->edges[i].t]])
		{
			if (ig[Index[g->edges[i].s]].id > ig[Index[g->edges[i].t]].id)
			{
				int tmp = g->edges[i].s;
				g->edges[i].s = g->edges[i].t;
				g->edges[i].t = tmp;
			}
		}

	}

	free(C);
	free(ig);
	free(d0);
	free(cd0);
	free(adj0);
}


//////////////////////////
//Building the special graph
graph* mkgraph(edgelist *el) {
	unsigned i, max;
	unsigned *d;
	graph* g = malloc(sizeof(graph));
	g->e = el->e;

	d = calloc(el->n, sizeof(unsigned));

	for (i = 0; i < el->e; i++) {
		d[el->edges[i].s]++;
	}

	g->cd = malloc((el->n + 1) * sizeof(unsigned));
	g->cd[0] = 0;
	max = 0;
	for (i = 1; i < el->n + 1; i++) {
		g->cd[i] = g->cd[i - 1] + d[i - 1];
		max = (max > d[i - 1]) ? max : d[i - 1];
		d[i - 1] = 0;
	}
	printf("max truncated degree = %u\n", max);

	g->adj = malloc(el->e * sizeof(unsigned));

	for (i = 0; i < el->e; i++) {
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

	sg->cd = malloc((g->n + 1) * sizeof(unsigned));
	sg->adj = malloc((g->core * g->core > g->e ? g->e : g->core * g->core) * sizeof(unsigned));
	for (i = 2; i < k; i++) {
		sg->d[i] = malloc(g->core * sizeof(unsigned));
		sg->nodes[i] = malloc(g->core * sizeof(unsigned));
	}

	sg->lab = calloc(g->core, sizeof(unsigned char));
	sg->color = malloc(g->core * sizeof(unsigned));
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
		for (i = 0; i < g->n; i++) {
			new[i] = -1;
		}
	}

	for (i = 0; i < sg->n[k - 1]; i++) {
		sg->lab[i] = 0;
	}

	j = 0;
	for (i = g->cd[u]; i < g->cd[u + 1]; i++) {
		v = g->adj[i];
		new[v] = j;
		old[j] = v;
		sg->lab[j] = k - 1;
		sg->nodes[k - 1][j] = j;
		sg->d[k - 1][j] = 0;//new degrees
		j++;
	}

	sg->n[k - 1] = j;
	//sg->color = malloc(j * sizeof(unsigned));

	int sub_edges = 0;
	for (i = 0; i < sg->n[k - 1]; i++) {//reodering adjacency list and computing new degrees
		v = old[i];
		sg->color[i] = color[Index[v]];

		/*
		for (l = g->cd[v]; l < g->cd[v + 1]; l++) {
			w = g->adj[l];
			j = new[w];
			if (j != -1) {
				sub_edges++;
			}
		}
		*/
	}
	//for (int i = 2; i < k; i++)
	//	sg->tmpadj[i] = malloc(sub_edges * sizeof(unsigned));

	sg->cd[0] = 0;
	for (i = 0; i < sg->n[k - 1]; i++) {
		v = old[i];
		sg->color[i] = color[Index[v]];

		for (l = g->cd[v]; l < g->cd[v + 1]; l++) {
			w = g->adj[l];
			j = new[w];
			if (j != -1) {
				sg->adj[sg->cd[i] + sg->d[k - 1][i]++] = j;
			}
		}
		sg->cd[i + 1] = sg->cd[i] + sg->d[k - 1][i];
	}

	for (i = g->cd[u]; i < g->cd[u + 1]; i++) {
		new[g->adj[i]] = -1;
	}
}

void kclique_thread(unsigned char l, subgraph *sg, unsigned long long *n) {
	unsigned i, j, k, end, u, v, w;

	if (l == 2) {
		for (i = 0; i < sg->n[2]; i++) {//list all edges
			u = sg->nodes[2][i];
			(*n) += sg->d[2][u];
			/*
			end = sg->cd[u] + sg->d[2][u];
			for (j = sg->cd[u]; j < end; j++) {
				(*n)++;//listing here!!!  // NOTE THAT WE COULD DO (*n)+=g->d[2][u] to be much faster (for counting only); !!!!!!!!!!!!!!!!!!
			}
			*/

		}
		return;
	}

	if (l > sg->n[l])
		return;

	for (i = 0; i < sg->n[l]; i++) {

		u = sg->nodes[l][i];
		if (sg->color[u] < l - 1)
			continue;

		sg->n[l - 1] = 0;
		end = sg->cd[u] + sg->d[l][u];
		for (j = sg->cd[u]; j < end; j++) {//relabeling nodes and forming U'.
			v = sg->adj[j];
			if (sg->lab[v] == l) {
				sg->lab[v] = l - 1;
				sg->nodes[l - 1][sg->n[l - 1]++] = v;
				sg->d[l - 1][v] = 0;//new degrees
			}
		}
		for (j = 0; j < sg->n[l - 1]; j++) {
			v = sg->nodes[l - 1][j];
			end = sg->cd[v] + sg->d[l][v];
			int Index = sg->cd[v];
			for (k = sg->cd[v]; k < end; k++) {
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

		for (j = 0; j < sg->n[l - 1]; j++) {//restoring labels
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
		for (u = 0; u < g->n; u++) {
			mksub(g, u, sg, k);
			kclique_thread(k - 1, sg, &n);
		}
		//free_subgraph(sg, k);

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
	ord_color_relabel(el);
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
