#include "core_and_truss_decomposition.h"
#include <ctype.h>
#include <sys/time.h>

Decom::Decom()
{
	E = 0; V = 0;
	maxcore = 0;
	maxtruss = 0;
	max_nbr = 0;

	Top_Size = 0;
	Max_Size = 0;
	top_num = 0;
	minsize_of_topnum = 0;

	global_color = 0;
}

Decom::~Decom()
{
}

void Decom::read_graph(const char * str)
{
	char BUFF[1024];
	FILE *in = fopen(str, "r");
	if (in == NULL) {
		cout << "No such the ‘in’ file " << endl;
		exit(1);
	}

	while (fgets(BUFF, 1024, in) != NULL) {
		if (isdigit(*BUFF))
			break;
	}
	sscanf(BUFF, "%ld%ld", &V, &E);

	vector<pair <int,int> > temp_edges;
	vertice.resize(V + 2);
	degrees.resize(V + 1);
	temp_edges.reserve(E);
	Long x, y;

	//for (int i = 0; i < V + 2; ++i)
	//	cout <<i <<" "<< ver[i] << endl;
	E = 0;
	while (fgets(BUFF, 1024, in) != NULL)
	{
		if (isdigit(*BUFF)) {
			sscanf(BUFF, "%ld%ld", &x, &y);
			if (x > y)
				continue;
			assert(x <= V && x >= 0);
			if (y > V || y < 0){
				cout << "y " << y << " E " << E << endl;
				assert(y <= V && y >= 0);
			}
			++vertice[x], ++vertice[y];
			temp_edges.emplace_back(x, y);
			++E;
		}
		if (feof(in))
			break;
	}
	fclose(in);
	//exit(1);
	cout << "E " << E << endl;
	E *= 2;
	edges.resize(E);
	reverse_pos.resize(E);

	Long counts = 0;

	for (Long i = 0, d; i <= V; ++i)
	{
		d = vertice[i];
		degrees[i] = d;
		if (max_nbr < d)
			max_nbr = d;
		vertice[i] = counts;
		counts += d;
		/*cout << i << " " << ver[i] << endl;*/
	}
	vertice[V + 1] = counts;

	for (size_t i = 0; i < temp_edges.size(); ++i)
	{
		x = temp_edges[i].first;
		y = temp_edges[i].second;

		edges[vertice[x]] = y;
		edges[vertice[y]] = x;

		reverse_pos[vertice[x]] = vertice[y];
		reverse_pos[vertice[y]] = vertice[x];

		++vertice[x];
		++vertice[y];
	}

	for (Long i = V + 1; i > 0; --i)
		vertice[i] = vertice[i - 1];
	vertice[0] = 0;

	printf("|V|=%ld, |E|=%ld\n", V, E);
	printf("max nbr: %ld\n", max_nbr);
}

void Decom::core()
{
	coreness.resize(V + 1);
	for (Long i = 0; i <= V; ++i)
		coreness[i] = degrees[i];

	Long v, u, t, s, d, i;
	vector<Long> pos, bin, core_ordering;
	pos.resize(V + 1);
	bin.resize(max_nbr + 1);
	core_ordering.resize(V + 1);


	for (i = 0; i <= V; ++i)
		++bin[coreness[i]];

	for (i = 0, t = 0; i <= max_nbr; ++i) {
		s = bin[i];
		bin[i] = t;
		t += s;
		if (t > V + 1) {
			cout << s << endl;
			exit(1);
		}
	}

	for (i = 0, s = 0; i <= V; ++i) {
		s = bin[coreness[i]]++;
		core_ordering[s] = i;
		pos[i] = s;
	}

	for (i = max_nbr; i > 0; --i)
		bin[i] = bin[i - 1];
	bin[0] = 0;

	int old_pos, new_pos;
	for (i = 0; i <= V; ++i) {
		u = core_ordering[i];
		d = coreness[u];
		s = vertice[u]; t = vertice[u + 1];
		while (s < t) {
			v = edges[s++];
			if (coreness[v] > d) {
				old_pos = pos[v];
				new_pos = bin[coreness[v]];

				if (old_pos != new_pos) {
					core_ordering[old_pos] = core_ordering[new_pos];
					core_ordering[new_pos] = v;

					pos[v] = new_pos;
					pos[core_ordering[old_pos]] = old_pos;
				}
				++bin[coreness[v]]; --coreness[v];
			}
		}
	}

	maxcore = coreness[core_ordering[V]];

	cout << "max core: " << maxcore << endl;
}

int Decom::com_nbr(Long u, Long v)
{
	Long us, ut, vs, vt, nu, nv, len = 0;
	us = vertice[u];
	ut = vertice[u + 1];
	vs = vertice[v];
	vt = vertice[v + 1];
	while (us < ut && vs < vt)
	{
		nu = edges[us];
		nv = edges[vs];
		if (nu > nv)
			++vs;
		else if (nu < nv)
			++us;
		else {
			++len;
			++us; ++vs;
		}
	}
	return len;
}

int Decom::find_nbr(Long u, Long v)
{
	Long ns, nt, w, mid;

	ns = vertice[u];
	nt = vertice[u + 1] - 1;
	while (ns <= nt)
	{
		mid = (ns + nt) / 2;
		w = edges[mid];
		if (w > v)
			nt = mid - 1;
		else if (w < v)
			ns = mid + 1;
		else
			return mid;
	}
	return -1;
}

void Decom::truss()
{
	queue< pair<Long, Long> > Q;
	Long i, u, v, ns, nt, k, countE = 0, esize;
	Long us, ut, vs, vt, nu, nv, pos;
	trussness.resize(E);

	for (i = 0; i <= V; ++i) {
		ns = vertice[i];
		nt = vertice[i + 1];
		while (ns < nt)
		{
			v = edges[ns];
			if (i < v) {
				Long tt = com_nbr(i, v) + 2;
				trussness[ns] = tt;
				//cout <<"("<< i <<","<< v <<"): "<< tt <<endl;

				trussness[reverse_pos[ns]] = tt;
			}
			++ns;
		}
	}
	//exit(1);
	k = 2;
	esize = E / 2;
	while (countE < esize)
	{
		for (i = 0; i <= V; ++i) {
			ns = vertice[i];
			nt = vertice[i + 1];
			while (ns < nt)
			{
				v = edges[ns];
				if (i < v && trussness[ns] == k)
					Q.emplace(i, v);
				++ns;
			}
		}
		//cout << "Q " << Q.size() << endl;
		while (!Q.empty())
		{
			u = Q.front().first;
			v = Q.front().second;
			Q.pop();
			++countE;
			/*if (u > v){
			cout << "u > v" <<endl;
			}*/
			pos = find_nbr(u, v);
			trussness[pos] = 0 - trussness[pos];
			trussness[reverse_pos[pos]] = trussness[pos];

			us = vertice[u];
			ut = vertice[u + 1];
			vs = vertice[v];
			vt = vertice[v + 1];
			while (us < ut && vs < vt)
			{
				if (us < ut && vs < vt) {
					nu = edges[us];
					nv = edges[vs];
					if (nu > nv)
						++vs;
					else if (nu < nv)
						++us;
					else {

						if (trussness[us] >= k && trussness[vs] >= k) {

							if (trussness[us] > k) {
								--trussness[reverse_pos[us]];
								if (--trussness[us] == k)
									Q.emplace(u, nu);
							}

							if (trussness[vs] > k) {
								--trussness[reverse_pos[vs]];
								if (--trussness[vs] == k)
									Q.emplace(v, nv);
							}
						}

						++us; ++vs;
					}
				}
			}

		}
		if (countE >= esize)
			maxtruss = k;
		++k;
	}
	cout << "max turss: " << maxtruss << endl;
	maxtruss = 0;
	for (i = 0; i <= V; ++i) {
		ns = vertice[i];
		nt = vertice[i + 1];
		while (ns < nt)
		{
			v = edges[ns];
			trussness[ns] = 0 - trussness[ns];
			if (maxtruss < trussness[ns])
				maxtruss = trussness[ns];
			++ns;
		}
	}
	//cout << "max turss: " << maxtruss << endl;
}

void Decom::coloring()
{
	if (max_nbr > V)
	{
		printf("coloring(): max_nbr > V\n");
		exit(1);
	}
	Long i, j, counts, dnums;
	Long u, v, us, ut, c;
	vector<Long> ver_sorted, maxd, color_temp;
	colors.resize(V + 1);
	ver_sorted.resize(V + 1);
	maxd.resize(max_nbr + 1);

	for (i = 0; i <= V; ++i) {
		++maxd[degrees[i]];
		colors[i] = -1;
	}

	for (i = 0, counts = 0, dnums; i <= max_nbr; ++i)
	{
		dnums = maxd[i];
		maxd[i] = counts;
		counts += dnums;
	}

	//sort
	for (i = 0; i <= V; ++i)
	{
		dnums = degrees[i];
		ver_sorted[maxd[dnums]++] = i;
	}
	/*for (i = 0; i < V; ++i) {
		int d1 = degrees[ver_sorted[i]];
		int d2 = degrees[ver_sorted[i + 1]];
		if (d1 > d2)
		{
			cout << "not sorted" << endl;
			exit(1);
		}
	}*/

	global_color = 0;
	color_temp.resize(max_nbr);
	for (i = V; i >= 0; --i)
	{
		u = ver_sorted[i];
		us = vertice[u];
		ut = vertice[u + 1];

		for (j = 0; j < global_color; ++j)
			color_temp[j] = 0;
		while (us < ut)
		{
			v = edges[us++];
			c = colors[v];
			if (c >= 0)
				color_temp[c] = 1;
		}
		for (j = 0; j < global_color; ++j)
			if (color_temp[j] == 0) {
				break;
			}
		if (j >= global_color)
			colors[u] = global_color++;
		else
			colors[u] = j;
	}
	cout << "colors: " << global_color << endl;
}

//void Decom::RDS()
//{
//	Long i, c_size;
//	bool found = false;
//	vector<Long> R, omega, Can;
//	R.resize(max_nbr + 1);
//	omega.resize(V + 1);
//
//	for (i = 0; i <= V; ++i)
//	{
//		found = false;
//		c_size = Get_Can(i, i, Can);
//		/*cout << i << " nbr size " << c_size << endl;
//		for (int i = 0; i < c_size; ++i) {
//		cout << Can[i].first << " " << Can[i].second << endl;
//		}
//		cout << endl;*/
//		R[0] = i;
//		//search(Can, c_size, R, 1, found);
//		//omega[i] = MAX_SIZE;
//	}
//}
//
//Long Decom::Get_Can(Long u, Long lm, vector<Long>& C)
//{
//	Long us, ut, w, len = 0;
//	us = vertice[u]; ut = vertice[u + 1];
//	while (us < ut)
//	{
//		w = edges[us++];
//		if (w <= lm)
//			C[len++] = w;
//	}
//	return len;
//}

void Decom::SetTopSize(Long Tsize){
	Top_Size = Tsize;
	//cout << "Top_Size:\t" << Top_Size <<endl;
}

void Decom::BK()
{
	if (Top_Size == 0)
		Top_Size = 1;
	Max_C.resize(Top_Size);

	Long i;
	vector<Long> Can, X, R;
	Can.reserve(V + 1);
	X.resize(V + 1);
	R.resize(max_nbr + 1);
	for (i = 0; i <= V; ++i)
		Can.emplace_back(i);

	Search_Clique(Can, V + 1, R, 0, X, 0);

	cout << "max clique: " << Max_Size << endl;
}

void Decom::Search_Clique(vector<Long>& Can, Long c_size, vector<Long>& R, Long r_size, vector<Long>& X, Long x_size)
{
	if (c_size == 0 && x_size == 0) {

		if (Top_Size > top_num) {
			
			Max_C[top_num].resize(r_size);
			for (int i = 0; i < r_size; ++i)
				Max_C[top_num][i] = R[i];
			++top_num;
			if (minsize_of_topnum == 0)
				minsize_of_topnum = r_size;
			else if (minsize_of_topnum > r_size)
				minsize_of_topnum = r_size;
			if (top_num == Top_Size)
				Max_Size = minsize_of_topnum;
		}
		else if (r_size > Max_Size) {
			bool flag = true;
			minsize_of_topnum = V;
			for (Long m_size, i = 0; i < Top_Size; ++i) {
				m_size = Max_C[i].size();
				if (m_size == Max_Size && flag) {
					flag = false;
					Max_C[i].resize(r_size);
					for (int j = 0; j < r_size; ++j)
						Max_C[i][j] = R[j];
					m_size = r_size;
				}
				if (minsize_of_topnum > m_size)
					minsize_of_topnum = m_size;
			}
			Max_Size = minsize_of_topnum;
		}
		return;
	}
	if (c_size == 0)
		return;

	Long i, j, u, v, q_size = 0, d = -1;
	Long color_nums, c;
	vector<Long> Q, Can_n, X_n, ver_colors;
	Can_n.resize(c_size);
	Q.resize(c_size);
	ver_colors.resize(global_color);
	for (i = 0; i < c_size; ++i) {
		v = Can[i];
		if (degrees[v] > d) {
			d = degrees[v];
			u = v;
		}
	}
	for (i = 0, color_nums = 0; i < c_size; ++i) {
		c = colors[Can[i]];
		if (ver_colors[c]++ == 0)
			++color_nums;
	}
	get_difference(u, Can,c_size, Q, q_size);
	Long cn_size = 0, xn_size = 0;
	X_n.resize(x_size + c_size);
	c_size; i = 0; j = 0 ;
	while (i < q_size)
	{
		if (c_size - i + r_size < Max_Size)
			break;
		if (color_nums + r_size < Max_Size)
			break;
		u = Q[i++];
		while (Can[j] != u && j < c_size)
			++j;
		if (--ver_colors[colors[u]] <= 0)
			--color_nums;
		Can[j] = -1;
		R[r_size] = u;

		get_intersectionP(Can, c_size, u, Can_n, cn_size);
		get_intersectionX(X, x_size, u, X_n, xn_size);

		Search_Clique(Can_n, cn_size, R, r_size + 1, X_n, xn_size);

		X[x_size++] = u;

	}
}

void Decom::get_difference(Long u, vector<Long>& Can, Long c_size, vector<Long>& Q, Long & q_size)
{
	Long us, ut, w, v, i = 0;
	us = vertice[u]; ut = vertice[u + 1];
	q_size = 0;
	while (us < ut && i < c_size)
	{
		w = edges[us];
		v = Can[i];
		if (w < v)
			++us;
		else if (w > v) {
			Q[q_size++] = v;
			++i;
		}
		else {
			++us; ++i;
		}
	}
	while (i < c_size)
		Q[q_size++] = Can[i++];
}

void Decom::get_intersectionX(vector<Long>& X, Long x_size, Long u, vector<Long>& X_n, Long & xn_size)
{
	Long us, ut, i, m, n;
	us = vertice[u]; ut = vertice[u + 1];
	i = 0; xn_size = 0;
	while (us < ut && i < xn_size)
	{
		m = edges[us];
		n = X[i];
		if (m < n)
			++us;
		else if (m > n)
			++i;
		else {
			X_n[xn_size++] = m;
			++us; ++i;
		}
	}
}

void Decom::get_intersectionP(vector<Long>& Can, Long c_size, Long u, vector<Long>& Can_n, Long & cn_size)
{
	Long us, ut, m, n, i;
	us = vertice[u]; ut = vertice[u + 1];
	cn_size = 0; i = 0;
	while (us < ut && i < c_size)
	{
		n = Can[i];
		if (n == -1) {
			++i;
			continue;
		}
		m = edges[us];
		if (m > n)
			++i;
		else if (m < n)
			++us;
		else {
			Can_n[cn_size++] = n;
			++us; ++i;
		}
	}
}

void Decom::print(const char *str, Long alg)
{
	// alg = 1 print maximum core 
	// alg = 2 print maximum truss
	// alg = 3 print top cliques
	Long i, us, ut, v, x, len;
	char Buff[128], p[128], c;
	vector<int> vcore;
	
	FILE *out = NULL;

	len = strlen(str);
	for (i = 0, x = 0; i < len; ++i) {
		c = str[i];
		if (c == '\\' || c == '/')
			x = 0;
		if (c == '.' && x != 0)
			break;
		if ((c <= 'z'&& c >= 'a') || (c <= 'Z'&& c >= 'A')
			|| (c <= '9'&& c >= '0') || c == '-' || c == '_')
			p[x++] = c;
	}
	p[x] = '\0';
	read_graph(str);

	timeval t_start, t_end;
	
	switch (alg)
	{
	case 1:
		gettimeofday(&t_start, NULL);
		core();
		gettimeofday(&t_end, NULL);
		cout << "core times " << (t_end.tv_sec - t_start.tv_sec) * 1000 + (t_end.tv_usec - t_start.tv_usec) / 1000
		<< " ms"<< endl; 

		sprintf(Buff, "%s_core_%ld.txt", p, maxcore);
		cout << Buff << endl;
		out = fopen(Buff, "w");
		if (out == NULL) {
			cout << "maximum core outfile errors" << endl;
			exit(1);
		}
		vcore.resize(V + 1);
		for (i = 0; i <= V; ++i) {
			if (coreness[i] == maxcore)
				vcore[i] = 1;
		}
		for (i = 0; i <= V; ++i) {
			if (vcore[i] == 1) {
				us = vertice[i];
				ut = vertice[i + 1];
				while (us < ut)
				{
					v = edges[us++];
					if (i < v && vcore[v] == 1) {
						fprintf(out, "%ld\t%ld\n", i, v);
					}
				}
			}
		}
		break;
	case 2:
		gettimeofday(&t_start, NULL);
		truss();
		gettimeofday(&t_end, NULL);
		cout << "truss times " << (t_end.tv_sec - t_start.tv_sec) * 1000 + (t_end.tv_usec - t_start.tv_usec) / 1000
		<< " ms"<< endl; 
		sprintf(Buff, "%s_truss_%ld.txt", p, maxtruss);
		cout << Buff << endl;
		out = fopen(Buff, "w");
		if (out == NULL) {
			cout << "maximum truss outfile errors" << endl;
			exit(1);
		}

		for (i = 0; i <= V; ++i) {
			us = vertice[i];
			ut = vertice[i + 1];
			while (us < ut)
			{
				v = edges[us];
				if (i < v && trussness[us] == maxtruss) {
					fprintf(out, "%ld\t%ld\n", i, v);
				}
				us++;
			}
		}
		break;
	case 3:
		if (Top_Size == 0)
			Top_Size = 1;
			
		sprintf(Buff, "%s_Top%ld_clique.txt", p, Top_Size);
		cout << Buff << endl;
		/*out = fopen(Buff, "w");
		if (out == NULL) {
			cout << "Top "<< Top_Size << " clique outfile errors" << endl;
			exit(1);
		}
		*/

		gettimeofday(&t_start, NULL);
		coloring();
		BK();
		gettimeofday(&t_end, NULL);
		cout << "BK times " << (t_end.tv_sec - t_start.tv_sec) * 1000 + (t_end.tv_usec - t_start.tv_usec) / 1000
		<< " ms"<< endl; 
		/*
		for (i = 0; i < Top_Size; ++i) {
			len = Max_C[i].size();
			fprintf(out, "clique %ld size:%ld\n", i + 1, len);
			for (Long j = 0; j < len - 1; ++j) {
				for (Long l = j + 1; l < len; ++l)
					fprintf(out, "%ld\t%ld\n", Max_C[i][j], Max_C[i][l]);
			}
			fprintf(out, "\n");
		}
		*/
		break;
	default:
		cout << "alg != (1, 2, 3)" <<endl;
		break;
	}
}
