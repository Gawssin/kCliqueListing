#include <iostream>
#include <cmath>

#include <time.h>

#include "core_and_truss_decomposition.h"

using namespace std;

int Gcd(int a, int b)
{
	if (b == 0)
		return a;
	return Gcd(b, a % b);
}

int main(int argc, char *argv[])
{
	const char *str = "../datas/dblp.out_180";
	const char *str1 = "../datas/as20000102.txt";
	clock_t  start, end;
	Decom d;

	//parameters:(filePath, alg) 如果是计算clique则可以再加入一个参数topsize
	switch (argc){
	case 4:
		if (atol(argv[2]) != 3){
			cout << "alg " << argv[2] << " no need TopSize" <<endl;
			exit(1);
		}
		d.SetTopSize(atol(argv[3]));
	case 3:
		d.print(argv[1],atol(argv[2]));
		break;
	default:
		cout << "parameter errors" <<endl;
		break;
	}
	return 0;

	//parameters:(filePath, TopSize) 不需要输出版本
	if (argc == 3){
		d.read_graph(argv[1]);
		d.SetTopSize(atol(argv[2]));
	}
	else if (argc == 2)//parameters:(filePath)
		d.read_graph(argv[1]);
	else
		d.read_graph(str1);
	start = clock();
	d.core();
	d.truss();
	d.coloring();
	start = clock();
	d.BK();
	end = clock();
	cout << "time:\t" << end - start <<endl;

	/*int x, y;
	cout << "please in the values" << endl;
	cin >> x >> y;

	int a, b, temp;
	temp = Gcd(x, x % y);
	x /= temp;
	y /= temp;
	a = 0; b = 0;
	temp = y;
	while (temp % 2 == 0){
		++a;
		temp /= 2;
	}
	temp = y;
	while (temp % 5 == 0){
		++b;
		temp /= 5;
	}
	a = a > b ? a : b;
	temp = x;
	while (a-- > 0)
		temp = (temp * 10) % y;
	b = temp;

	a = 10000;
	while (a && temp != 0){
		cout << temp * 10 / y;
		temp = temp * 10 % y;
		--a;
		if (temp == b)
			break;
	}
	cout << endl;
	cout << "len " << 10000 - a << endl;*/

	return 0;
}