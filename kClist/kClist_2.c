/*

Info:
Feel free to use these lines as you wish.
This program iterates over all k-cliques.
This is an improvement of the 1985 algorithm of Chiba And Nishizeki detailed in "Arboricity and subgraph listing".

To compile:
"gcc kClist.c -O9 -o kClist".

To execute:
"./kClist k edgelist.txt".
"edgelist.txt" should contain the graph: one edge on each line separated by a space.
Will print the number of k-cliques.
*/


#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>

#define NLINKS 100000000 //maximum number of edges for memory allocation, will increase if needed

typedef struct {
	unsigned s;
	unsigned t;
} edge;

typedef struct {
	unsigned node;
	unsigned deg;
} nodedeg;


typedef struct {

	unsigned n;//number of nodes
	unsigned e;//number of edges
	edge *edges;//list of edges

	unsigned *ns;//ns[l]: number of nodes in G_l
	unsigned **d;//d[l]: degrees of G_l
	unsigned *cd;//cumulative degree: (starts with 0) length=n+1
	unsigned *adj;//truncated list of neighbors

	unsigned char *lab;//lab[i] label of node i
	unsigned **sub;//sub[l]: nodes in G_l

} specialsparse;


typedef struct {
	unsigned id;
	unsigned degree;
} iddegree;


int *color;
unsigned *Index, **tmpadj;
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

void freespecialsparse(specialsparse *g, unsigned char k) {
	unsigned char i;
	free(g->ns);
	for (i = 2; i < k + 1; i++) {
		free(g->d[i]);
		free(g->sub[i]);
	}
	free(g->d);
	free(g->sub);
	free(g->edges);
	free(g->lab);
	free(g->cd);
	free(g->adj);
	free(g);
}

//Compute the maximum of three unsigned integers.
unsigned int max3(unsigned int a, unsigned int b, unsigned int c);
inline unsigned int max3(unsigned int a, unsigned int b, unsigned int c) {
	a = (a > b) ? a : b;
	return (a > c) ? a : c;
}

specialsparse* readedgelist(char* edgelist) {
	unsigned e1 = NLINKS;
	specialsparse *g = malloc(sizeof(specialsparse));
	FILE *file;

	g->n = 0;
	g->e = 0;
	file = fopen(edgelist, "r");
	g->edges = malloc(e1 * sizeof(edge));
	while (fscanf(file, "%u %u", &(g->edges[g->e].s), &(g->edges[g->e].t)) == 2) {//Add one edge
		g->n = max3(g->n, g->edges[g->e].s, g->edges[g->e].t);
		g->e++;
		if (g->e == e1) {
			e1 += NLINKS;
			g->edges = realloc(g->edges, e1 * sizeof(edge));
		}
	}
	fclose(file);
	g->n++;

	g->edges = realloc(g->edges, g->e * sizeof(edge));

	return g;
}


//computing degeneracy ordering and core value
void ord_color_relabel(specialsparse* g) {
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
		int tmpdegree = ig[i].degree, tmpid=ig[i].id;
		for (int j = 0; j < tmpdegree; j++)
		{
			int now = Index[ adj0[ cd0[tmpid]+j  ] ];
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
			if( ig[Index[g->edges[i].s]].id > ig[Index[g->edges[i].t]].id)
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
//Building the special graph structure
void mkspecial(specialsparse *g, unsigned char k) {
	unsigned i, ns, max;
	unsigned *d, *sub;
	unsigned char *lab;

	d = calloc(g->n, sizeof(unsigned));

	for (i = 0; i < g->e; i++) {
		d[g->edges[i].s]++;
	}

	g->cd = malloc((g->n + 1) * sizeof(unsigned));
	ns = 0;
	g->cd[0] = 0;
	max = 0;
	sub = malloc(g->n * sizeof(unsigned));
	lab = malloc(g->n * sizeof(unsigned char));
	for (i = 1; i < g->n + 1; i++) {
		g->cd[i] = g->cd[i - 1] + d[i - 1];
		max = (max > d[i - 1]) ? max : d[i - 1];
		sub[ns++] = i - 1;
		d[i - 1] = 0;
		lab[i - 1] = k;
	}
	printf("max degree = %u\n", max);

	g->adj = malloc(g->e * sizeof(unsigned));

	for (i = 0; i < g->e; i++) {
		g->adj[g->cd[g->edges[i].s] + d[g->edges[i].s]++] = g->edges[i].t;
	}

	for (int i = 0; i < g->n; i++)
	{
		qsort( g->adj + g->cd[i],d[i],sizeof(unsigned),cmpadj);
	}

	g->ns = malloc((k + 1) * sizeof(unsigned));
	g->ns[k] = ns;

	g->d = malloc((k + 1) * sizeof(unsigned*));
	g->sub = malloc((k + 1) * sizeof(unsigned*));
	tmpadj = malloc((k + 1) * sizeof(unsigned*));
	for (i = 2; i <= k; i++) {
		g->d[i] = malloc(g->n * sizeof(unsigned));
		g->sub[i] = malloc(max * sizeof(unsigned));
		tmpadj[i] = malloc(g->e * sizeof(unsigned));
	}
	g->d[k] = d;
	qsort(sub, g->n, sizeof(unsigned), cmpadj);
	g->sub[k] = sub;
	tmpadj[k] = g->adj;
	g->lab = lab;
}


void kclique(unsigned l, specialsparse *g, unsigned long long *n) {
	unsigned i, j, k, end, u, v, w;

	if (l == 2) {
		for (i = 0; i < g->ns[2]; i++) {//list all edges
			u = g->sub[2][i];
			(*n)+=g->d[2][u];
			/*
			end = g->cd[u] + g->d[2][u];
			for (j = g->cd[u]; j < end; j++) {
				(*n)++;//listing here!!!  // NOTE THAT WE COULD DO (*n)+=g->d[2][u] to be much faster (for counting only); !!!!!!!!!!!!!!!!!!
			}
			*/
		}
		return;
	}

	if (l > g->ns[l])
		return;
	for (i = 0; i < g->ns[l]; i++) {
		

		u = g->sub[l][i];
		if (color[Index[u]] < l-1)
			break;
		//printf("%u %u\n",i,u);
		g->ns[l - 1] = 0;
		end = g->cd[u] + g->d[l][u];
		for (j = g->cd[u]; j < end; j++) {//relabeling nodes and forming U'.
			v = tmpadj[l][j];
			if (g->lab[v] == l) {
				g->lab[v] = l - 1;
				g->sub[l - 1][g->ns[l - 1]++] = v;
				g->d[l - 1][v] = 0;//new degrees
			}
		}
		for (j = 0; j < g->ns[l - 1]; j++) {//reodering adjacency list and computing new degrees
			v = g->sub[l - 1][j];
			end = g->cd[v] + g->d[l][v];
			int Index = g->cd[v];
			for (k = g->cd[v]; k < end; k++) {
				w = tmpadj[l][k];
				if (g->lab[w] == l - 1) {
					g->d[l - 1][v]++;
					tmpadj[l - 1][Index++] = w;
				}
			}
		}

		kclique(l - 1, g, n);

		for (j = 0; j < g->ns[l - 1]; j++) {//restoring labels
			v = g->sub[l - 1][j];
			g->lab[v] = l;
		}

	}
}


int main(int argc, char** argv) {
	specialsparse* g;
	int k = atoi(argv[1]);
	unsigned long long n;
	time_t t0, t1, t2;
	t1 = time(NULL);
	t0 = t1;

	printf("Reading edgelist from file %s\n", argv[2]);

	g = readedgelist(argv[2]);
	printf("Number of nodes = %u\n", g->n);
	printf("Number of edges = %u\n", g->e);

	t2 = time(NULL);
	printf("- Time = %ldh%ldm%lds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
	t1 = t2;

	printf("Building the graph structure\n");

	ord_color_relabel(g);
	

	mkspecial(g, k);

	printf("Number of nodes = %u\n", g->n);
	printf("Number of edges = %u\n", g->e);

	t2 = time(NULL);
	printf("- Time = %ldh%ldm%lds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
	t1 = t2;

	printf("Iterate over all cliques\n");

	n = 0;
	kclique(k, g, &n);

	printf("Number of %u-cliques: %llu\n", k, n);

	t2 = time(NULL);
	printf("- Time = %ldh%ldm%lds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
	t1 = t2;

	freespecialsparse(g, k);
	free(color);
	free(Index);

	printf("- Overall time = %ldh%ldm%lds\n", (t2 - t0) / 3600, ((t2 - t0) % 3600) / 60, ((t2 - t0) % 60));

	return 0;
}
