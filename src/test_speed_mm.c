#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "functions.h"

#include "mkl.h"

/*
 *	Purpose:
 *		Test runtimes for various matrix sizes used in mm_*.
 *
 *	Comments:
 *		As all of the significantly expensive matrix-matrix/matrix-vector operations will be done using a column-major
 *		orientation, only these tests are performed.
 *		Experiment with column-major vs. row-major layout for the operator matrix and check for any performance
 *		differences. (ToBeDeleted)
 *
 *		Important conclusions from time testing:
 *			ToBeModified.
 *
 *	Notation:
 *
 *  References:
 *
 */

void test_speed_mm_d(void)
{
	printf("test_speed_mm_d:\n\n");

	register unsigned int i, iMax;
	unsigned int m, n, k, NIter;
	double *A, *B, *C;

	clock_t ts, te;
	double tt;

	NIter = 1e4;
	m = 4;
	n = 3;
	k = 2;

	A = malloc(m*k * sizeof *A);
	B = malloc(k*n * sizeof *B);
	C = malloc(m*n * sizeof *C);

	for (i = 0, iMax = m*k; iMax--; i++)
		A[i] = (double)rand()/(double)RAND_MAX;

	for (i = 0, iMax = k*n; iMax--; i++)
		B[i] = (double)rand()/(double)RAND_MAX;

	for (i = 0, iMax = m*n; iMax--; i++)
		C[i] = (double)rand()/(double)RAND_MAX;

	ts = clock();
	for (iMax = NIter; iMax--; ) {
		mm_d(CblasColMajor,CblasTrans,CblasNoTrans,m,n,k,1.0,A,B,C);
	}
	te = clock();
	tt = (te-ts)/(double) CLOCKS_PER_SEC;

	printf("  Time to run %d mm_d calls using m = %d, n = %d, k = %d: %10.6f\n",NIter,m,n,k,tt);

	// Write the mm_d routine using c code only.








}
