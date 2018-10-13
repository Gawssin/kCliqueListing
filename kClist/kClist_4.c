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
	unsigned *cd,*cdsub;//cumulative degree: (starts with 0) length=n+1
	unsigned *adj,*adjsub;//truncated list of neighbors
	unsigned *rank;//ranking of the nodes according to degeneracy ordering
				   //unsigned *map;//oldID newID correspondance

	unsigned char *lab;//lab[i] label of node i
	unsigned **sub;//sub[l]: nodes in G_l

} specialsparse;


void freespecialsparse(specialsparse *g, unsigned char k) {
	unsigned char i;
	free(g->ns);
	for (i = 2; i < k + 1; i++) {
		free(g->d[i]);
		free(g->sub[i]);
	}
	free(g->d);
	free(g->sub);
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

void relabel(specialsparse *g) {
	unsigned i, source, target, tmp;

	for (i = 0; i < g->e; i++) {
		source = g->rank[g->edges[i].s];
		target = g->rank[g->edges[i].t];
		if (source < target) {
			tmp = source;
			source = target;
			target = tmp;
		}
		g->edges[i].s = source;
		g->edges[i].t = target;
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
void ord_core(specialsparse* g) {
	unsigned i, j, r = 0, n = g->n;
	keyvalue kv;
	bheap *heap;

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
		d0[i - 1] = 0;
	}
	for (i = 0; i < g->e; i++) {
		adj0[cd0[g->edges[i].s] + d0[g->edges[i].s]++] = g->edges[i].t;
		adj0[cd0[g->edges[i].t] + d0[g->edges[i].t]++] = g->edges[i].s;
	}

	heap = mkheap(n, d0);

	g->rank = malloc(g->n * sizeof(unsigned));
	for (i = 0; i < g->n; i++) {
		kv = popmin(heap);
		g->rank[kv.key] = n - (++r);
		for (j = cd0[kv.key]; j < cd0[kv.key + 1]; j++) {
			update(heap, adj0[j]);
		}
	}
	for (i = 0; i < g->n; i++)
		printf("i = %d rank = %d\n",i, g->rank[i]);

	freeheap(heap);
	free(d0);
	free(cd0);
	free(adj0);
}

//////////////////////////
//Building the special graph structure
unsigned *dsub;
int *color;
unsigned *index;
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
	dsub = calloc(max, sizeof(unsigned));
	index = malloc(max * sizeof(unsigned));
	color = malloc(max * sizeof(int));


	g->adj = malloc(g->e * sizeof(unsigned));

	for (i = 0; i < g->e; i++) {
		g->adj[g->cd[g->edges[i].s] + d[g->edges[i].s]++] = g->edges[i].t;
	}
	free(g->edges);

	g->ns = malloc((k + 1) * sizeof(unsigned));
	g->ns[k] = ns;

	g->d = malloc((k + 1) * sizeof(unsigned*));
	g->sub = malloc((k + 1) * sizeof(unsigned*));
	for (i = 2; i < k; i++) {
		g->d[i] = malloc(g->n * sizeof(unsigned));
		g->sub[i] = malloc(max * sizeof(unsigned));
	}
	g->d[k] = d;
	g->sub[k] = sub;

	g->lab = lab;
}

int K;

typedef struct {
	unsigned id;
	unsigned degree;
} iddegree;

int cmp(const void* a, const void* b)
{
	// qsort'cmp 可以 return 0和负�?or 正数 
	iddegree *x = (iddegree*)a, *y = (iddegree*)b;

	return y->degree - x->degree;
}


unsigned **tmpadj;

int cmpadj(const void* a, const void* b)
{
	// qsort'cmp 可以 return 0和负�?or 正数 
	int *x = (int*)a, *y = (int*)b;

	return color[index[*y]] - color[index[*x]];
}
void mkspecial_sub(specialsparse *g, unsigned char k) {
	unsigned i, ns, max;
	unsigned *d, *sub;
	unsigned char *lab;
	printf("hello k = %d n=%d!\n",k,g->n);
	d = calloc(g->n, sizeof(unsigned));
	//printf("hello kclist222  %d %d!\n", g->n,g->e);
	for (i = 0; i < g->e; i++) {
		//printf("i = %d   g->edges[i].s = %d!\n", i, g->edges[i].s);
		d[g->edges[i].s]++;
		
	}
	//printf("hello kclist2223!\n");
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
	//printf("sub max degree = %u\n", max);
	

	g->adj = malloc(g->e * sizeof(unsigned));

	for (i = 0; i < g->e; i++) {
		g->adj[g->cd[g->edges[i].s] + d[g->edges[i].s]++] = g->edges[i].t;
	}


	//qsort(g->adj, d[0], sizeof(unsigned), cmpadj);
	for (int i = 0; i < g->n; i++)
	{
		qsort(g->adj + g->cd[i], d[i], sizeof(unsigned), cmpadj);
	}


	free(g->edges);


	//printf("sub max\n");
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
	//printf("sub max 2\n");
	g->d[k] = d;
	qsort(sub, g->n, sizeof(unsigned), cmpadj);
	g->sub[k] = sub;
	tmpadj[k] = g->adj;
	g->lab = lab;
	//printf("sub max 3\n");
}


void kclique(unsigned l, specialsparse *g, unsigned long long *n) {
	unsigned i, j, k, end, u, v, w;
	
	if (l == 2) {
		for (i = 0; i < g->ns[2]; i++) {//list all edges
			u = g->sub[2][i];
			//(*n)+=g->d[2][u];

			end = g->cd[u] + g->d[2][u];
			for (j = g->cd[u]; j < end; j++) {
				(*n)++;//listing here!!!  // NOTE THAT WE COULD DO (*n)+=g->d[2][u] to be much faster (for counting only); !!!!!!!!!!!!!!!!!!
			}

		}
		return;
	}
	//printf("l ====== g->ns[l]=%d %u\n", l, g->ns[l]);
	for (i = 0; i < g->ns[l]; i++) {
		//printf("ll ====== %u\n", l);
		u = g->sub[l][i];
		//printf("%u %u ------------------------------\n",i,u);
		g->ns[l - 1] = 0;
		end = g->cd[u] + g->d[l][u];
		for (j = g->cd[u]; j < end; j++) {//relabeling nodes and forming U'.
			v = g->adj[j];
			if (g->lab[v] == l) {		//equal to if(1)
				g->lab[v] = l - 1;
				g->sub[l - 1][g->ns[l - 1]++] = v;
				g->d[l - 1][v] = 0;//new degrees
			}
		}
		//printf("lll ====== %u\n", l);
		
		if (l == K)
		{
			//printf("%u\n", i);
			//printf("0000000000000000000000000000\n");
			memset(dsub,0, g->ns[l-1] * sizeof(unsigned));
			unsigned *ind = malloc(g->n * sizeof(unsigned));

			unsigned *cd0 = malloc((g->ns[l - 1] + 1) * sizeof(unsigned));
			//unsigned *adj0 = malloc(2 * g->e * sizeof(unsigned));


			memset(ind, -1, g->n * sizeof(unsigned));


			//printf("hello kclist5!\n");
			int cnt = -1,edge_num=0;
			for (j = 0; j < g->ns[l - 1]; j++) 
			{//reodering adjacency list and computing new degrees

				v = g->sub[l - 1][j];
				if (ind[v] == -1)
					ind[v] = ++cnt;
				end = g->cd[v] + g->d[l][v];
				for (k = g->cd[v]; k < end; k++) 
				{
					w = g->adj[k];
					
					
					if (g->lab[w] == l - 1) 
					{
						if (ind[w] == -1)
							ind[w] = ++cnt;
						edge_num++;
						//g->d[l - 1][v]++;
						dsub[ind[v]]++;
						dsub[ind[w]]++;
					}
				}


			}
			//printf("cnt = %d\n", cnt);

			//printf("hello kclist6!\n");
			iddegree *ig = malloc(g->ns[l - 1] * sizeof(iddegree));

			unsigned *adj0 = malloc(2 * edge_num * sizeof(unsigned));

			cd0[0] = 0;
			for (int i = 1; i < g->ns[l - 1] + 1; i++) {
				cd0[i] = cd0[i - 1] + dsub[i - 1];
				ig[i - 1].id = i - 1;
				ig[i - 1].degree = dsub[i - 1];
				dsub[i - 1] = 0;
			}

			specialsparse *subg = malloc(sizeof(specialsparse));
			subg->edges = malloc(edge_num * sizeof(edge));

			//printf("hello g->ns[l - 1] = %d !\n", g->ns[l - 1]);
			
			for (j = 0; j < g->ns[l - 1]; j++)
			{//reodering adjacency list and computing new degrees

				v = g->sub[l - 1][j];
				end = g->cd[v] + g->d[l][v];
				//printf("hello kclist666!\n");
				for (k = g->cd[v]; k < end; k++)
				{
					w = g->adj[k];
					//printf(" w = %d  g->lab[w] = %d!\n",w, g->lab[w]);
					if (g->lab[w] == l - 1)
					{

						//printf("hello kclist 777  %d %d %d %d!\n", ind[v], ind[w], cd0[ind[v]], dsub[ind[v]]);
						adj0[cd0[ind[v]] + dsub[ind[v]]++] = ind[w];
						//printf("hello kclist 888!\n");
						adj0[cd0[ind[w]] + dsub[ind[w]]++] = ind[v];
					}
				}
			}



			qsort(ig, g->ns[l - 1], sizeof(ig[0]), cmp);

			
			for (int i = 0; i < g->ns[l - 1]; i++) 
				index[ ig[i].id ] = i;

			//printf("hello %d!\n", g->ns[l - 1]);
			memset(color, -1, sizeof(int)*g->ns[l - 1]);

			//printf("hello 1!\n");
			int *C = calloc((ig[0].degree + 1) , sizeof(int));
			//printf("hello 11!\n");
			int aa = 1;
			//printf("hello 2!\n");
  			//memset(C, 0, sizeof(int)*(ig[0].degree + 1));
			//printf("hello 3!\n");
			color[0] = 0;
			//printf("aab !\n");
			int colorNum = 0;
			//printf("aaa !\n");
			if (i == 3)
				printf("g->ns[l - 1] = %d\n", g->ns[l - 1]);
			int vv = i;
			if(i == 3)
			for (int i = 0; i < g->ns[l - 1]; i++)
				printf("id = %d dg = %d\n", ig[i].id,ig[i].degree);

			for (int i = 1; i < g->ns[l - 1]; i++)
			{
				//printf("loop!\n");
				int tmpdegree = ig[i].degree, tmpid = ig[i].id;
				if(vv ==3)
				printf("ig[i].degree = %d tmpid = %d!\n", ig[i].degree, tmpid);

				for (int j = 0; j < tmpdegree; j++)
				{
					int now = index[adj0[cd0[tmpid] + j]];
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
					int now = index[adj0[cd0[tmpid] + j]];
					if (color[now] != -1)
						C[color[now]] = 0;
				}

			}
			if(i == 3)
				printf("color number = %d\n", colorNum);

			if (i == 3)
			for (int i = 0; i < g->ns[l - 1]; i++)
				printf("color = %d\n", color[i]);



			int e_num = 0;
			for (j = 0; j < g->ns[l - 1]; j++)
			{//reodering adjacency list and computing new degrees

				v = g->sub[l - 1][j];
				end = g->cd[v] + g->d[l][v];
				for (k = g->cd[v]; k < end; k++)
				{
					w = g->adj[k];
					if (g->lab[w] == l - 1)
					{
						printf("%d\n", 1111111);
						if (color[index[ind[v]]] < color[index[ind[w]]])
						{

							subg->edges[e_num].s = ind[w];
							subg->edges[e_num++].t = ind[v];
						}
						else if (color[index[ind[v]]] == color[index[ind[w]]])
						{
							if (ig[index[ind[v]]].id < ig[index[ind[w]]].id)
							{
								subg->edges[e_num].s = ind[v];
								subg->edges[e_num++].t = ind[w];
							}
							else
							{
								subg->edges[e_num].s = ind[w];
								subg->edges[e_num++].t = ind[v];
							}
						}
						else if (color[index[ind[v]]] > color[index[ind[w]]])
						{
							subg->edges[e_num].s = ind[v];
							subg->edges[e_num++].t = ind[w];
						}
					}
				}
			}


			for (int i = 0; i < e_num; i++)
			{
				printf("s = %d t = %d \n", subg->edges[i].s,subg->edges[i].t);

			}


			//printf("hello kclist!\n");
			subg->n = g->ns[l - 1];
			subg->e = edge_num;
			mkspecial_sub(subg, l - 1);

			//printf("hello kclist2   l = %d!\n",l);
			kclique(l - 1, subg, n);
			//printf("hello kclist3!\n");
			for (j = 0; j < g->ns[l - 1]; j++) {//restoring labels
				v = g->sub[l - 1][j];
				g->lab[v] = l;
			}

			free(ind);
			free(cd0);
			free(ig);
			//printf("hello kclist5!\n");
			//printf("hello kclist a!\n");
			free(subg);
			//printf("hello kclist q!\n");
			free(adj0);
			//printf("hello kclist w!\n");
			//free(color);
			//printf("hello kclist e!\n");
			//free(index);
			//printf("hello kclist r!\n");
			free(C);

		}

		else
		{
			//printf("subsub %d \n",l);
			printf("l = %d ns = %d\n",l,g->ns[l]);
			if (l > g->ns[l])
				return;
			//printf("aaaaa ns = %d\n", g->ns[l]);
			for (int i = 0; i < g->ns[l]; i++) {

				//printf("aaaaa i=%d ns = %d\n",i, g->ns[l]);
				u = g->sub[l][i];
				printf("color = %d\n", color[index[u]]);
				if (color[index[u]] < l - 1)
					break;
				//printf("%u %u\n",i,u);
				g->ns[l - 1] = 0;
				end = g->cd[u] + g->d[l][u];
				printf("d = %d\n", g->d[l][u]);
				for (j = g->cd[u]; j < end; j++) {//relabeling nodes and forming U'.
					v = tmpadj[l][j];
					printf("tmpadj[l][j] =  %d g->lab[v] = %d \n", tmpadj[l][j], g->lab[v]);
					if (g->lab[v] == l) {
						printf("label\\\\\\\\\\\\\\\n");
						g->lab[v] = l - 1;
						g->sub[l - 1][g->ns[l - 1]++] = v;
						g->d[l - 1][v] = 0;//new degrees
					}
				}
				for (j = 0; j < g->ns[l - 1]; j++) {//reodering adjacency list and computing new degrees
					v = g->sub[l - 1][j];
					end = g->cd[v] + g->d[l][v];
					int index = g->cd[v];
					for (k = g->cd[v]; k < end; k++) {
						w = tmpadj[l][k];
						if (g->lab[w] == l - 1) {
							g->d[l - 1][v]++;
							tmpadj[l - 1][index++] = w;
						}
						/*
						else {
						g->adj[k--] = g->adj[--end];
						g->adj[end] = w;
						}
						*/
					}
				}
				//printf("wwww l = %d\n",l);
				kclique(l - 1, g, n);

				for (j = 0; j < g->ns[l - 1]; j++) {//restoring labels
					v = g->sub[l - 1][j];
					g->lab[v] = l;
				}

			}
		}

		

	}
}


int main(int argc, char** argv) {
	//sym unweighted
	//2766607 1965206
	//freopen("out.log", "w", stdout);
	specialsparse* g;
	unsigned char k = atoi(argv[1]);
	K = k;
	unsigned long long n;
	time_t t0, t1, t2;
	t1 = time(NULL);
	t0 = t1;

	printf("Reading edgelist from file %s\n", argv[2]);

	g = readedgelist(argv[2]);
	printf("Number of nodes = %u\n", g->n);
	printf("Number of edges = %u\n", g->e);

	t2 = time(NULL);
	printf("- Time = %lldh%lldm%llds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
	t1 = t2;

	printf("Building the graph structure\n");

	ord_core(g);
	relabel(g);

	mkspecial(g, k);

	printf("Number of nodes = %u\n", g->n);
	printf("Number of edges = %u\n", g->e);

	t2 = time(NULL);
	printf("- Time = %lldh%lldm%llds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
	t1 = t2;

	printf("Iterate over all cliques\n");

	n = 0;
	kclique(k, g, &n);

	printf("Number of %u-cliques: %llu\n", k, n);

	t2 = time(NULL);
	printf("- Time = %lldh%lldm%llds\n", (t2 - t1) / 3600, ((t2 - t1) % 3600) / 60, ((t2 - t1) % 60));
	t1 = t2;

	freespecialsparse(g, k);

	printf("- Overall time = %lldh%lldm%llds\n", (t2 - t0) / 3600, ((t2 - t0) % 3600) / 60, ((t2 - t0) % 60));
	free(dsub);
	return 0;
}
