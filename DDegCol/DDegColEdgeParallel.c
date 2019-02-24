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
	unsigned *rank;//ranking of the nodes according to degeneracy ordering
				   //unsigned *map;//oldID newID correspondance NOT USED IN THIS VERSION
} edgelist;

typedef struct {
	unsigned n;
	unsigned e;
	edge *edges;
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
} subgraph;

typedef struct {
	unsigned id;
	unsigned degree;
} iddegree;

typedef struct {
	unsigned id;
	unsigned color;
} idcolor;

int cmp(const void* a, const void* b)
{
	iddegree *x = (iddegree*)a, *y = (iddegree*)b;
	return y->degree - x->degree;
}

int cmpadj(const void* a, const void* b)
{
	idcolor *x = (idcolor*)a, *y = (idcolor*)b;
	return y->color - x->color;
}

void free_edgelist(edgelist *el) {
	free(el->edges);
	free(el->rank);
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
	for (i = 2; i < k; i++) {
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

void relabel(edgelist *el) {
	unsigned i, source, target, tmp;

	for (i = 0; i < el->e; i++) {
		source = el->rank[el->edges[i].s];
		target = el->rank[el->edges[i].t];
		if (source < target) {
			tmp = source;
			source = target;
			target = tmp;
		}
		el->edges[i].s = source;
		el->edges[i].t = target;
	}

}

///// CORE ordering /////////////////////
typedef struct {
	unsigned key;
	unsigned value;
} keyvalue;

typedef struct {
	unsigned n_max;	// max number of nodes.
	unsigned n;	// number of nodes.
	unsigned *pt;	// pointers to nodes.
	keyvalue *kv; // nodes.
} bheap;


bheap *construct(unsigned n_max) {
	unsigned i;
	bheap *heap = malloc(sizeof(bheap));

	heap->n_max = n_max;
	heap->n = 0;
	heap->pt = malloc(n_max * sizeof(unsigned));
	for (i = 0; i < n_max; i++) heap->pt[i] = -1;
	heap->kv = malloc(n_max * sizeof(keyvalue));
	return heap;
}

void swap(bheap *heap, unsigned i, unsigned j) {
	keyvalue kv_tmp = heap->kv[i];
	unsigned pt_tmp = heap->pt[kv_tmp.key];
	heap->pt[heap->kv[i].key] = heap->pt[heap->kv[j].key];
	heap->kv[i] = heap->kv[j];
	heap->pt[heap->kv[j].key] = pt_tmp;
	heap->kv[j] = kv_tmp;
}

void bubble_up(bheap *heap, unsigned i) {
	unsigned j = (i - 1) / 2;
	while (i > 0) {
		if (heap->kv[j].value > heap->kv[i].value) {
			swap(heap, i, j);
			i = j;
			j = (i - 1) / 2;
		}
		else break;
	}
}

void bubble_down(bheap *heap) {
	unsigned i = 0, j1 = 1, j2 = 2, j;
	while (j1 < heap->n) {
		j = ((j2 < heap->n) && (heap->kv[j2].value < heap->kv[j1].value)) ? j2 : j1;
		if (heap->kv[j].value < heap->kv[i].value) {
			swap(heap, i, j);
			i = j;
			j1 = 2 * i + 1;
			j2 = j1 + 1;
			continue;
		}
		break;
	}
}

void insert(bheap *heap, keyvalue kv) {
	heap->pt[kv.key] = (heap->n)++;
	heap->kv[heap->n - 1] = kv;
	bubble_up(heap, heap->n - 1);
}

void update(bheap *heap, unsigned key) {
	unsigned i = heap->pt[key];
	if (i != -1) {
		((heap->kv[i]).value)--;
		bubble_up(heap, i);
	}
}

keyvalue popmin(bheap *heap) {
	keyvalue min = heap->kv[0];
	heap->pt[min.key] = -1;
	heap->kv[0] = heap->kv[--(heap->n)];
	heap->pt[heap->kv[0].key] = 0;
	bubble_down(heap);
	return min;
}

//Building the heap structure with (key,value)=(node,degree) for each node
bheap* mkheap(unsigned n, unsigned *v) {
	unsigned i;
	keyvalue kv;
	bheap* heap = construct(n);
	for (i = 0; i < n; i++) {
		kv.key = i;
		kv.value = v[i];
		insert(heap, kv);
	}
	return heap;
}

void freeheap(bheap *heap) {
	free(heap->pt);
	free(heap->kv);
	free(heap);
}

//computing degeneracy ordering and core value
void ord_core(edgelist* el) {
	unsigned i, j, r = 0, n = el->n, e = el->e;
	keyvalue kv;
	bheap *heap;

	unsigned *d0 = calloc(el->n, sizeof(unsigned));
	unsigned *cd0 = malloc((el->n + 1) * sizeof(unsigned));
	unsigned *adj0 = malloc(2 * el->e * sizeof(unsigned));
	for (i = 0; i < e; i++) {
		d0[el->edges[i].s]++;
		d0[el->edges[i].t]++;
	}
	cd0[0] = 0;
	for (i = 1; i < n + 1; i++) {
		cd0[i] = cd0[i - 1] + d0[i - 1];
		d0[i - 1] = 0;
	}
	for (i = 0; i < e; i++) {
		adj0[cd0[el->edges[i].s] + d0[el->edges[i].s]++] = el->edges[i].t;
		adj0[cd0[el->edges[i].t] + d0[el->edges[i].t]++] = el->edges[i].s;
	}

	heap = mkheap(n, d0);

	el->rank = malloc(n * sizeof(unsigned));
	for (i = 0; i < n; i++) {
		kv = popmin(heap);
		el->rank[kv.key] = n - (++r);
		for (j = cd0[kv.key]; j < cd0[kv.key + 1]; j++) {
			update(heap, adj0[j]);
		}
	}

	freeheap(heap);
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
	printf("core value (max truncated degree) = %u\n", max);

	g->adj = malloc(el->e * sizeof(unsigned));

	for (i = 0; i < el->e; i++) {
		g->adj[g->cd[el->edges[i].s] + d[el->edges[i].s]++] = el->edges[i].t;
	}

	free(d);
	g->core = max;
	g->edges = el->edges;
	g->e = el->e;
	g->n = el->n;
	return g;
}


subgraph* allocsub(graph *g, unsigned char k) {
	unsigned i;
	subgraph* sg = malloc(sizeof(subgraph));
	sg->n = calloc(k, sizeof(unsigned));
	sg->d = malloc(k * sizeof(unsigned*));
	sg->nodes = malloc(k * sizeof(unsigned*));
	sg->adj = malloc(g->core*g->core * sizeof(unsigned));
	for (i = 0; i < k; i++) {
		sg->d[i] = malloc(g->core * sizeof(unsigned));
		sg->nodes[i] = malloc(g->core * sizeof(unsigned));
	}
	sg->lab = calloc(g->core, sizeof(unsigned char));

	sg->core = g->core;
	return sg;
}

void mksub(graph* g, edge ed, subgraph* sg, unsigned char k) {
	unsigned i, j, l, x, y;
	unsigned u = ed.s, v = ed.t;

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

	for (i = g->cd[v]; i<g->cd[v + 1]; i++) {
		new[g->adj[i]] = -2;
	}

	j = 0;
	for (i = g->cd[u]; i<g->cd[u + 1]; i++) {
		x = g->adj[i];
		if (new[x] == -2) {
			new[x] = j;
			old[j] = x;
			sg->lab[j] = k - 2;
			sg->nodes[k - 2][j] = j;
			sg->d[k - 2][j] = 0;//new degrees
			j++;
		}
	}

	sg->n[k - 2] = j;

	unsigned *d0 = calloc(j, sizeof(unsigned));
	for (i = 0; i < sg->n[k - 2]; i++) {//reodering adjacency list and computing new degrees
		v = old[i];
		for (l = g->cd[v]; l < g->cd[v + 1]; l++) {
			y = g->adj[l];
			j = new[y];
			if (j < -2) {
				sg->adj[sg->core*i + sg->d[k - 2][i]++] = j;
				sg->adj[sg->core*j + sg->d[k - 2][j]++] = i;
				d0[i]++;
				d0[j]++;

			}
		}
	}

	//color ordering
	unsigned *C = calloc(sg->n[k - 2], sizeof(unsigned));
	int *color = malloc(sg->n[k - 2] * sizeof(int));
	unsigned *Index = malloc(sg->n[k - 2] * sizeof(unsigned));
	iddegree *ig;
	ig = malloc(sg->n[k - 2] * sizeof(iddegree));
	for (i = 0; i < sg->n[k - 2]; i++)
	{
		color[i] = -1;
		ig[i].id = i;
		ig[i].degree = d0[i];
	}
	qsort(ig, sg->n[k - 2], sizeof(ig[0]), cmp);

	for (i = 0; i < sg->n[k - 2]; i++)
		Index[ig[i].id] = i;

	color[0] = 0;
	int colorNum = 0;

	for (int i = 1; i < sg->n[k - 2]; i++)
	{
		int tmpdegree = ig[i].degree, tmpid = ig[i].id;

		for (int j = 0; j < tmpdegree; j++)
		{
			int now = Index[sg->adj[sg->core*tmpid + j]];
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
			int now = Index[sg->adj[sg->core*tmpid + j]];
			if (color[now] != -1)
				C[color[now]] = 0;
		}

	}

	sg->color = malloc(sg->n[k - 2] * sizeof(unsigned));
	for (int i = 0; i < sg->n[k - 2]; i++)
	{
		sg->d[k - 2][i] = 0;
		sg->color[i] = color[Index[i]];
	}

	//relabel
	for (i = 0; i < sg->n[k - 2]; i++) {
		v = old[i];
		for (l = g->cd[v]; l < g->cd[v + 1]; l++) {
			y = g->adj[l];
			j = new[y];
			if (j < -2) {

				if (color[Index[i]] > color[Index[j]])
				{
					sg->adj[sg->core*i + sg->d[k - 2][i]++] = j;
				}
				else
				{
					sg->adj[sg->core*j + sg->d[k - 2][j]++] = i;
				}

			}
		}
	}

	for (i = g->cd[ed.t]; i < g->cd[ed.t + 1]; i++) {
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
			end = u*sg->core + sg->d[2][u];
			for (j = u*sg->core; j < end; j++) {
				(*n)++;//listing here!!!  // NOTE THAT WE COULD DO (*n)+=g->d[2][u] to be much faster (for counting only); !!!!!!!!!!!!!!!!!!
			}
			*/
		}
		return;
	}
	if (l == 1)
	{
		(*n) += sg->n[l];
		return;
	}

	if (l > sg->n[l])
		return;
	for (i = 0; i < sg->n[l]; i++) {
		u = sg->nodes[l][i];

		if (sg->color[u] < l - 1)
			continue;

		sg->n[l - 1] = 0;
		end = u*sg->core + sg->d[l][u];
		for (j = u*sg->core; j < end; j++) {//relabeling nodes and forming U'.
			v = sg->adj[j];
			if (sg->lab[v] == l) {
				sg->lab[v] = l - 1;
				sg->nodes[l - 1][sg->n[l - 1]++] = v;
				sg->d[l - 1][v] = 0;//new degrees
			}
		}
		for (j = 0; j < sg->n[l - 1]; j++) {
			v = sg->nodes[l - 1][j];
			end = sg->core*v + sg->d[l][v];
			int Index = sg->core*v;
			for (k = sg->core*v; k < end; k++) {
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
		for (u = 0; u < g->e; u++) {
			mksub(g, g->edges[u], sg, k);
			kclique_thread(k - 2, sg, &n);
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
	ord_core(el);
	relabel(el);
	g = mkgraph(el);

	printf("Number of nodes (degree > 0) = %u\n", g->n);

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
