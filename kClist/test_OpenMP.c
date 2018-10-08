#include <omp.h>
#include <stdio.h>
int main(int argc, char *argv[])
{
	int nthreads, tid;
	int nprocs;
	char buf[32];
	printf("result = %d\n",-1<-2);
	/*Fork a team of threads*/
#pragma omp parallel private(nthreads,tid)
	{

		/* Obtain and print thread id  */
		tid = omp_get_thread_num();
		printf("Hello, world from OpenMP thread %d\n", tid);   /*Only master thread does this  */
		if (tid == 0)
		{
			nthreads = omp_get_num_threads();
			printf(" Number of threads %d\n", nthreads);
		}
	}return 0;
}
