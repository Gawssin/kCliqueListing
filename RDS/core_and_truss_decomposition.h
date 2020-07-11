#pragma once
#ifndef CORE_AND_TRUSS_DECOMPOSITION_H

#include <iostream>
#include <vector>
#include <queue>

#include <assert.h>
#include <string.h>

using namespace std;

#define Long long int

class Decom
{
private:
	Long E, V;
	Long maxcore, maxtruss, max_nbr, global_color;
	vector<Long> vertice, edges, reverse_pos;

	vector<Long> coreness, degrees, colors;
	vector<Long> trussness;

	Long Top_Size, Max_Size, top_num, minsize_of_topnum;
	vector< vector<Long> > Max_C;

public:
	Decom();
	~Decom();

	void read_graph(const char *str);

	void core();
	int com_nbr(Long u, Long v);
	int find_nbr(Long u, Long v);
	void truss();

	void coloring();

	//void RDS();
	//Long Get_Can(Long u, Long lm, vector<Long> & C);

	void SetTopSize(Long Tsize);
	void BK();
	void Search_Clique(vector<Long> &Can, Long c_size, vector<Long> &R, Long r_size, vector<Long> &X, Long x_size);
	void get_difference(Long u, vector<Long> &Can, Long c_size, vector<Long> &Q, Long &q_size);
	void get_intersectionX(vector<Long>& X, Long x_size, Long u, vector<Long>& X_n, Long & xn_size);
	void get_intersectionP(vector<Long> &Can, Long c_size, Long u, vector<Long> &Can_n, Long &cn_size);

	void print(const char *str, Long alg);
};


#endif // !CORE_AND_TRUSS_DECOMPOSITION_H