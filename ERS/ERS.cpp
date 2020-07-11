#include <stdlib.h>
#include <cstdio>
#include <stdbool.h>
#include <string.h>
#include <string>
#include <vector>
#include <utility>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <time.h>
#include <cmath>
#include <random>
#include <queue>
#include <cstdarg>
#include <chrono>
#include <cassert>

using namespace std;
using namespace chrono;

class edge
{
public:
	edge(int s = 0, int t = 0) :s(s), t(t) {}

	int s;
	int t;
};

class iddeg
{
public:
	int id;
	int degree;
};

class Graph
{
public:
	Graph();
	Graph(const Graph& obj);
	~Graph();
	void readedgelist(string edgelist);
	void coreDecomposition();
	void mkGraph();
	int outLargeClique();
	bool isEdge(int, int);
	Graph* mksub(int, ...);
	//Graph* mksubMark(int*, int, int*);
	int color(int*);
	bool isCliqueV(int *, int);
	void orderByDegree();
	//bool isEdge(int, int);

	int n;
	int m;
	int maxDeg, maxDegv;
	edge* edges;

	int* deg, *degv;
	int* cd, *cdv;
	int* adj, *adjv;
	int* coreRank;			//increasing core number order
	int* coreNum;			//coreNum[i] is the core number of node i.
	int* bin;
};

double nCr[1001][101];

inline int max3(int a, int b, int c) {
	a = (a > b) ? a : b;
	return (a > c) ? a : c;
}

Graph::Graph(void) {}
Graph::~Graph(void)
{
	if (deg != NULL) delete[] deg;
	if (cd != NULL) delete[] cd;
	if (adj != NULL) delete[] adj;
	if (degv != NULL) delete[] degv;
	if (cdv != NULL) delete[] cdv;
	if (adjv != NULL) delete[] adjv;
	if (coreRank != NULL) delete[] coreRank;
	if (coreNum != NULL) delete[] coreNum;
	if (bin != NULL) delete[] bin;
}
Graph::Graph(const Graph& obj)
{
	n = obj.n, m = obj.m, maxDeg = obj.maxDeg, edges = obj.edges;
	edges = obj.edges;

	if (deg != NULL) delete[] deg;
	if (obj.deg != NULL) deg = new int[n], memcpy(deg, obj.deg, n * sizeof(int));

	if (cd != NULL) delete[] cd;
	if (obj.cd != NULL) cd = new int[n + 1], memcpy(cd, obj.cd, (n + 1) * sizeof(int));

	if (adj != NULL) delete[] adj;
	if (obj.adj != NULL) adj = new int[2 * m], memcpy(adj, obj.adj, 2 * m * sizeof(int));

	if (coreRank != NULL) delete[] coreRank;
	if (obj.coreRank != NULL) coreRank = new int[n], memcpy(coreRank, obj.coreRank, n * sizeof(int));

	if (coreNum != NULL) delete[] coreNum;
	if (obj.coreNum != NULL) coreNum = new int[n], memcpy(coreNum, obj.coreNum, n * sizeof(int));

	if (bin != NULL) delete[] bin;
	if (obj.bin != NULL) bin = new int[maxDeg + 2], memcpy(bin, obj.bin, (maxDeg + 2) * sizeof(int));
}

void Graph::readedgelist(string edgelist)
{
	ifstream file;
	file.open(edgelist);
	file >> n >> m;
	edges = new edge[m];
	m = 0;
	while (file >> edges[m].s >> edges[m].t) m++;
	file.close();
}

void Graph::mkGraph()
{
	deg = new int[n]();
	cd = new int[n + 1];
	adj = new int[2 * m];
	maxDeg = 0;
	for (int i = 0; i < m; i++)
	{
		deg[edges[i].s]++;
		deg[edges[i].t]++;
		maxDeg = max3(maxDeg, deg[edges[i].s], deg[edges[i].t]);
	}
	cd[0] = 0;
	for (int i = 1; i < n + 1; i++)
	{
		cd[i] = cd[i - 1] + deg[i - 1];
		deg[i - 1] = 0;
	}

	for (int i = 0; i < m; i++)
	{
		adj[cd[edges[i].s] + deg[edges[i].s]++] = edges[i].t;
		adj[cd[edges[i].t] + deg[edges[i].t]++] = edges[i].s;
	}

	for (int i = 0; i < n; i++) sort(adj + cd[i], adj + cd[i] + deg[i]);
}

bool cmp(const pair<int, int>& a, const pair<int, int>& b)
{
	return a.second > b.second;
}
bool IGCmp(const iddeg& a, const iddeg& b)
{
	return a.degree == b.degree ? (a.id < b.id) : (a.degree > b.degree);
}

bool Graph::isEdge(int a, int b)
{
	if (deg[a] > deg[b]) a = a ^ b, b = a ^ b, a = a ^ b;
	for (int i = cd[a]; i < cd[a] + deg[a]; i++)
		if (adj[i] == b) return true;
	return false;
}


void Graph::coreDecomposition()
{
	bin = new int[maxDeg + 2]();

	for (int i = 0; i < n; i++)
		bin[deg[i]]++;

	int lastBin = bin[0], nowBin;
	bin[0] = 0;
	for (int i = 1; i <= maxDeg; i++)
	{
		nowBin = lastBin + bin[i - 1];
		lastBin = bin[i];
		bin[i] = nowBin;
	}
	int* vert = new int[n](), *pos = new int[n](), *tmpDeg = new int[n]();
	for (int i = 0; i < n; i++)
	{
		pos[i] = bin[deg[i]];

		vert[bin[deg[i]]++] = i;
		tmpDeg[i] = deg[i];
	}

	bin[0] = 0;
	for (int i = maxDeg; i >= 1; i--)
	{
		bin[i] = bin[i - 1];
	}

	//int core = 0;
	int* cNum = new int[n];
	//int *cNum = (int *)malloc(g->n * sizeof(int));
	for (int i = 0; i < n; i++)
	{
		int id = vert[i], nbr, binFrontId;
		//if (i == bin[core + 1]) ++core;
		cNum[id] = tmpDeg[id];
		for (int i = cd[id]; i < cd[id] + deg[id]; i++)
		{
			nbr = adj[i];

			if (tmpDeg[nbr] > tmpDeg[id])
			{
				binFrontId = vert[bin[tmpDeg[nbr]]];
				if (binFrontId != nbr)
				{

					pos[binFrontId] = pos[nbr];
					pos[nbr] = bin[tmpDeg[nbr]];
					vert[bin[tmpDeg[nbr]]] = nbr;
					vert[pos[binFrontId]] = binFrontId;

				}
				bin[tmpDeg[nbr]]++;
				tmpDeg[nbr]--;

			}

		}

	}

	coreNum = cNum;

	coreRank = vert;

	delete[] tmpDeg;
	delete[] pos;
}


void Graph::orderByDegree()
{
	degv = new int[n]();
	cdv = new int[n];
	adjv = new int[m];
	//cout << "max 5 : " << maxDegv << endl;
	for (int i = 0; i < n; i++)
		for (int j = cd[i]; j < cd[i] + deg[i]; j++)
		{
			int v = adj[j];
			if ((deg[i] < deg[v]) || ((deg[i] == deg[v]) && (i < v))) degv[i]++;
				
		}
	//cout << "max 6 : " << maxDegv << endl;
	maxDegv = max(degv[n - 1], 0);
	cdv[0] = 0, degv[n - 1] = 0;
	for (int i = 1; i < n; i++)
	{
		maxDegv = max(degv[i - 1], maxDegv);
		cdv[i] = cdv[i - 1] + degv[i - 1];
		degv[i - 1] = 0;
	}
	cout << "max v : " << maxDegv << endl;

	for (int i = 0; i < n; i++)
		for (int j = cd[i]; j < cd[i] + deg[i]; j++)
		{
			int v = adj[j];
			if ((deg[i] < deg[v]) || ((deg[i] == deg[v]) && (i < v))) adjv[cdv[i] + degv[i]++] = v;
		}
}

bool Graph::isCliqueV(int *nbr, int k)
{
	if (k == 1) return true;
	int minV = n + 1, minDeg = maxDeg + 1, u;
	for (int i = 0; i < k; i++) minDeg = min(minDeg, deg[nbr[i]]);

	for (int i = 0; i < k; i++) if ((deg[nbr[i]] == minDeg) && (nbr[i] < minV))  u = i, minV = nbr[i];

	int nodeId = nbr[u], flagi, flagj;
	for (int i = 0; i < k; i++)
	{
		if (i == u) continue;
		flagj = 0;
		for (int j = cdv[nodeId]; j < cdv[nodeId] + degv[nodeId]; j++)
		{
			if (adjv[j] == nbr[i])
			{
				flagj = 1;
				break;
			}
		}
		if (flagj == 0) return false;
	}

	swap(nbr[u], nbr[k-1]);
	bool cFlag = isCliqueV(nbr, k - 1);
	swap(nbr[u], nbr[k - 1]);

	return cFlag;
}

long long factorial(int k)
{
	long long ans = 1;
	for (int i = 1; i <= k; i++) ans *= i;
	return ans;
}

void getnCr()
{
	FILE *fp;
	if ((fp = fopen("nCr.txt", "r")) == NULL)
	{
		puts("Fail to open file!");
		exit(0);
	}
	for (int i = 0; i < 1001; i++)
		for (int j = 0; j < 101; j++)
			fscanf(fp, "%lf", &nCr[i][j]);

	fclose(fp);
}

double EstimateClique(Graph *g, int k)
{
	getnCr();
	int S = 0, r = 1000000;
	double *pv = new double[g->n], W = 0;
	for (int i = 0; i < g->n; i++) pv[i] = nCr[g->degv[i]][k-1], W += pv[i];
	cout << "W: " << W << endl;
	vector<double> wts(pv, pv + g->n);

	random_device rd;  //Will be used to obtain a seed for the random number engine
	mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	discrete_distribution<size_t> d{ begin(wts), end(wts) };

	int *nbr = new int[k - 1], *nbrNum = new int[10000];
	for (int i = 0; i < 10000; i++) nbrNum[i] = i;

	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	for (int i = 0; i < r; i++)
	{
		int u = d(gen);
		if (g->degv[u] < k - 1) continue;
		vector<int> uadjv(nbrNum, nbrNum + g->degv[u]);
		random_shuffle(uadjv.begin(), uadjv.end());
		for (int j = 0; j < k-1; j++) nbr[j] = g->adjv[g->cdv[u] + uadjv[j]];

		if (g->isCliqueV(nbr, k - 1) == true) S++;
	}

	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	auto t12 = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

	cout << "a constant-factor estimate time: " << t12 / 1e6  << endl;
	
	cout << "a constant-factor estimate of " << k << "-clique: " << 1.0*S/r*W << endl;

	return 1.0*S / r * W;
}


edge D(Graph* g, int *S, int size, discrete_distribution<int> &d)
{
	random_device rd;  //Will be used to obtain a seed for the random number engine
	mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

	int u = S[d(gen)];

	uniform_int_distribution<int> dNbr(0, g->degv[u] - 1);

	int v = g->adjv[g->cdv[u] + dNbr(gen)];
	return edge(u, v);
}


bool Sample_Degrees_Typical(Graph *g, int k, int m, double Cke, double _eps, double _del, int *T, int &sizeT, int &mT, vector<int> &degT)
{
	double Gamma = min(_eps, 1.0/(4 * pow(m, k/2.0)));
	int end = log(2.0 / Gamma) / log(2) + 1;
	for (int i = 0; i < g->n; i++) T[i] = i;
	for (int i = 0; i < end; i++)
	{
		mT = 0;
		double tv = (10.0*(g->n)*log(2.0*(g->n) / Gamma / Gamma)) / (pow(_eps / k, 2) * 2 * sqrt(1.0*m)) + 1;
		int t = min(tv, g->n*1.0);
		random_shuffle(T, T + g->n);
		
		for (int j = 0; j < t; j++) mT += g->degv[T[j]];
		if (mT <= 4.0*m*t / g->n)
		{
			degT.resize(t);
			for (int j = 0; j < t; j++) degT[i] = g->degv[T[j]];
			sizeT = t;
			return true;
		}

	}
	return false;
}


bool Sample_a_Clique(Graph *g, int *S, int *T, int sizeS, int sizeT, int mS, int mT, int k, int *w, discrete_distribution<int> &dS, discrete_distribution<int> &dT)
{
	random_device rd;  //Will be used to obtain a seed for the random number engine
	mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	
	edge es = D(g, S, sizeS, dS);
	int u = es.s, v = es.t;

	w[0] = u, w[1] = v;
	for (int i = 0; i < k - 2; i++)
	{
		if (g->degv[v] * 1.0 <= 2 * sqrt(g->m))
		{
			double pw = 1.0 * g->degv[v] / (2 * sqrt(g->m)) ;
			if (pw >= 1) pw = 0.9;
			discrete_distribution<size_t> d{ pw, 1 - pw };
			if (d(gen) == 1) return false;
			uniform_int_distribution<int> dNbr(0, g->degv[v] - 1);
			int tm = dNbr(gen); 
			w[2 + i] = g->adjv[g->cdv[v] + tm];
		}
		else
		{
			edge et = D(g, T, sizeT, dT);
			int y = et.t;
			if (g->degv[y] * 1.0 <= 2 * sqrt(g->m)) return false;
			double pw = 1.0 *  mT/(g->degv[y] * sizeT * 2.0 * sqrt(g->m) / g->n);
			discrete_distribution<size_t> d{ pw, 1 - pw };
			if (d(gen) == 1) return false;
			w[2 + i] = y;
			cout << "here 2" << endl;
		}
	}

	for (int i = 0; i < k - 2; i++)
		if (g->deg[v] > g->deg[w[2 + i]] || ((g->deg[v] == g->deg[w[2 + i]]) && (v > w[2 + i]))) return false;

	bool tm = g->isCliqueV(w, k);
	return tm;
}


bool Is_sociable(Graph *g, int u, int *T, int sizeT, int mT, double _eps, double _del, double Cke, int k, discrete_distribution<int> &dS, discrete_distribution<int> &dT)
{
	if (g->degv[u] > (4.0*g->m / pow(_eps *  Cke, 1.0 / k))) return  false;
	int r = (g->degv[u] * pow(2 * sqrt(g->m), k - 2) * 15 * log(g->n / _del)) /
		(factorial(k - 2) * (k*pow(50 * Cke, 1 - 1.0 / k) / pow(_eps, 1.0 / k)) * pow(_eps, 2)) + 1;
	if (r < 0) cout << "Is_sociable error!" << endl;
	int *w = new int[k];
	discrete_distribution<int> du = { g->degv[u] *1.0 };

	int sum = 0;
	r /= 100;
	for (int i = 0; i < r; i++)
		sum += Sample_a_Clique(g, &u, T, 1, sizeT, g->degv[u], mT, k, w, du, dT);
	delete[] w;
	double Cku = (g->degv[u] * 2.0 * sqrt(g->m) * sum) / r;
	return  (Cku < 0.5 * k*pow(50 * Cke, 1 - 1.0 / k) / pow(_eps, 1.0 / k));
}

double Approximate_Cliques(Graph *g, int k, double Cke, double Epsilon, double Delta)
{
	double _eps = Epsilon / 5, _del = Delta / 4;

	int *S, mS = 0;
	double sizeSv = 1.0 * 700 * k*g->n*log(1 / _del) /
		(pow(_eps, 2.0 + 1.0 / k) * pow(Cke, 1.0 / k)) + 1;
	int sizeS = min(sizeSv, g->n*1.0);
	//sizeS /= 100;
	cout << "sizeS: " << sizeS << endl;

	S = new int[g->n];
	vector<int> degS, degT;
	degS.resize(sizeS);
	for (int i = 0; i < g->n; i++) S[i] = i;
	random_shuffle(S, S + g->n);
	for (int i = 0; i < sizeS; i++) mS += g->degv[S[i]], degS[i] = g->degv[S[i]];
	cout << "mS / sizeS : " << 1.0 * mS / sizeS << endl;
	discrete_distribution<int> dS{ begin(degS), end(degS) }, dT;

	int *T = new int[g->n], sizeT, mT;
	bool SDT = Sample_Degrees_Typical(g, k, g->m, Cke, _eps, _del, T, sizeT, mT, degT);
	dT = {begin(degT), end(degT)};
	if (!SDT)
	{
		cout << "Sample Degrees Typical Fail!" << endl;
		return 0;
	}
	double qv = (mS * pow(2 * sqrt(g->m), k - 2) * 10 * log(1 / _del)) /
		(pow(1 - _eps, 3) * factorial(k - 2) * Cke * (1.0 * sizeS / g->n) * pow(_eps, 2)) + 1;

	cout << "qv: " << qv << endl;
	int q = min(qv, g->m*1.0);
	//int q = g->m*5;
	q = 1e8;

	cout << "q: " << q << endl;

	double ESCkPerX = (mS * pow(2 * sqrt(g->m), k - 2) * 1 / q) / (factorial(k - 2) * 1.0 * sizeS / g->n) ;
	cout << "ES per X: " << ESCkPerX << endl;

	//bool *X = new bool[q];
	int *w = new int[k], sumX = 0;
	for (int i = 0; i < q; i++)
	{
		bool SCFlag = Sample_a_Clique(g, S, T, sizeS, sizeT, mS, mT, k, w, dS, dT);
		if (SCFlag == false) sumX += 0;
		else sumX += Is_sociable(g, w[0], T, sizeT, mT, _eps, _del, Cke, k, dS, dT);
	}
	cout << "sumX: " << sumX << endl;
	double ESCk = (mS * pow(2 * sqrt(g->m), k - 2) * sumX / q) / (factorial(k - 2) * 1.0 * sizeS / g->n) ;
	
	return ESCk;
}


int main(int argc, char** argv)
{
	Graph *g = new Graph;
	int k = atoi(argv[1]);
	cout << "Reading edgelist from file " << argv[2] << " k = " << k << endl;
	g->readedgelist(argv[2]);
	cout << "Reading edgelist finished!" << endl;
	cout << "N: " << g->n << endl << "M: " << g->m << endl;
	double Epsilon = atof(argv[3]), Delta = atof(argv[4]);

	g->mkGraph();
	//g->mv = g->m / 2;
	cout << "mkGraph finished!" << endl;
	g->orderByDegree();
	high_resolution_clock::time_point ES_t1 = high_resolution_clock::now();

	double Cke = EstimateClique(g, k);
	double ESC = Approximate_Cliques(g, k, Cke, Epsilon, Delta);
	cout << "Estimate number of " << k << "-clique: " << ESC << endl;
	high_resolution_clock::time_point ES_t2 = high_resolution_clock::now();
	auto ES_t1_2 = std::chrono::duration_cast<std::chrono::microseconds>(ES_t2 - ES_t1).count();
	cout << "toltal time: " << ES_t1_2 / 1e6 << endl;
	return 0;
}