// ConsoleApplication1.cpp: define el punto de entrada de la aplicación de consola.
//

/***OMP_EX1***/
/*
#include "stdafx.h"
#include <iostream>
#include <cstdio>
#include <omp.h>


using namespace std;

int main() {
#pragma omp parallel
	cout << "Hello world" << endl;
	//return 0;

	std::getchar();
}*/

/***OMP_EX2***/
/*
#include "stdafx.h"
#include <cstdio>
#include <stdio.h>
#include <omp.h>

void main(void)
{
#pragma omp parallel
	{
		int id = omp_get_thread_num();
		int thread_count = omp_get_num_threads();
		printf("Hello from the thread %d of %d \n", id, thread_count);
	}
	std::getchar();
}*/

/***OMP_EX3***/

/*
#include "stdafx.h"
#include <stdio.h>
#include <cstdio>
#include <omp.h>


void Hello(void); // Thread function that prints message

int main(int argc, char* argv[]) {
	// Number of threads
	int thread_count = 10;

	// Get the number of threads from command line
	// int thread_count = strtol(argv[1], NULL, 10);

#pragma omp parallel num_threads(thread_count)
	Hello();

	return 0;
} // Main

void Hello(void) {
	int my_rank = omp_get_thread_num();
	int thread_count = omp_get_num_threads();

	printf("Hello iTA from the thread %d of %d \n", my_rank, thread_count);
	std::getchar();
} // Hello
*/

/***OMP_EX4***/
/*
#include "stdafx.h"
#include <stdio.h>
#include <cstdio>
#include <omp.h>

double a[1000];

void Square(int id, double a[])
{
	int i;
	printf("thread id: %d\n", id);
	for (i = 0; i < 1000; i++)
		a[i] *= a[i];
}

void main(void) {

#pragma omp parallel num_threads(4)
	{
		int id = omp_get_thread_num();
		Square(id, a);
	}
	printf("all done\n");
	std::getchar();
}*/

/***OMP_EX5***/

#include "stdafx.h"
#include <stdio.h>
#include <cstdio>
#include <omp.h>

#define N 10

void main()
{
	double avg = 0.0, a[N];
	int i;
	int j;
	for (j = 0;j < N;j++)
		a[j] = j;

	printf("avg ini = %f \n", avg);

#pragma omp parallel for num_threads(N)
	for (i = 0; i < N; i++) {
		avg += a[i];
		printf("a[%i] = %f \n", i, a[i]);
		printf("parcial = %f \n", avg);
	}


	avg = avg / N;
	printf("Average = %f \n", avg);
	std::getchar();
}


/***OMP_EX6****/
/* Purpose: Estimate definite integral (or area under curve) using the
*          trapezoidal rule.  This version uses a parallel for directive
*/
/*
#include "stdafx.h"
#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

void Usage(char* prog_name);

double f(double x);
double Trap(double a, double b, int n, int thread_count);

int main(int argc, char* argv[]) {
	double  global_result = 0.0;
	double  a, b;
	int     n;
	int     thread_count;
	double aproxArea=0.0;

	if (argc != 2) Usage(argv[0]);
	thread_count = strtol(argv[1], NULL, 10);
	printf("Enter a, b, and n\n");
	scanf_s("%lf %lf %d", &a, &b, &n);

	aproxArea = Trap(a, b, n, thread_count);
	printf("area = %f \n", aproxArea);
	return 0;
	//std::getchar();
}  /* main */

/*
void Usage(char* prog_name) {

	fprintf(stderr, "usage: %s <number of threads>\n", prog_name);
	exit(0);
}


double f(double x) {
	double return_val;

	return_val = x*x;
	return return_val;
}

double Trap(double a, double b, int n, int thread_count) {
	double  h, approx;
	int  i;
	h = (b - a) / n;
	approx = (f(a) + f(b)) / 2.0;
#pragma omp parallel for num_threads(thread_count) reduction(+: approx)
	for (i = 1; i <= n - 1; i++)
		approx += f(a + i*h);
	
	approx = h*approx;

	return approx;
}*/
