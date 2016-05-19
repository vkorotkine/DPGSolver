// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "test.h"
#include "functions.h"
#include "parameters.h"

#include "mkl.h"

/*
 *	Purpose:
 *		Test speed of various implementations (with varying matrix sizes) of mm_CTN_*.
 *
 *	Comments:
 *		As all of the significantly expensive matrix-matrix/matrix-vector operations will be done using a column-major
 *		orientation, only these tests are performed.
 *		Despite being much less reader-friendly, it is certainly advantageous to manually unroll the loops:
 *			Both mv and mm are significantly faster than BLAS calls for "small" A (typical of the code).
 *			Both mv and mm are capable of achieving peak flops (mv BLAS ~50% of peak, mm BLAS ~100% peak)
 *
 *		Important conclusions from time testing:
 *			BLAS_mv breaks even when more than ~120 flops are performed. (ToBeModified)
 *			BLAS_mm breaks even when more than ~330 flops are performed. (ToBeModified)
 *
 *	Notation:
 *
 *	References:
 *
 */

void test_speed_mm_CTN(void)
{
	printf("\ntest_speed_mm_d:\n\n");

	register unsigned int i, iMax;
	unsigned int m, n, k, NIter, flops;
	double *A, *B, *C, *C0;

	clock_t ts, te;
	double tt;

	/*
	 *	mm_CTN_d:
	 *
	 *		Input:
	 *
	 *			A, B, m, n, k, C (only for memory location)
	 */

	// (m)atrix-(v)ector: m = k = 5
	NIter = pow(2,16);
	m = 5;
	n = 1;
	k = m;
	flops = m*n*(2*k-1);

	A  = malloc(m*k * sizeof *A);
	B  = malloc(k*n * sizeof *B);
	C  = malloc(m*n * sizeof *C);
	C0 = malloc(m*n * sizeof *C0);

	for (i = 0, iMax = m*k; iMax--; i++)
		A[i] = (double)rand()/(double)RAND_MAX;

	for (i = 0, iMax = k*n; iMax--; i++)
		B[i] = (double)rand()/(double)RAND_MAX;

	mm_CTN_d(m,n,k,A,B,C0);

	for (i = 0, iMax = m*n; iMax--; i++)
		C[i] = 0.0;
	ts = clock();
	for (iMax = NIter; iMax--; ) {
		mm_CTN_mv1_d(m,n,k,A,B,C);
	}
	te = clock();
	tt = (te-ts)/(double) CLOCKS_PER_SEC;

	printf("\nflops (%d) x NIter (%d) = Total flops (%d)\n",flops,NIter,flops*NIter);
	printf("mm_CTN_mv1_d (m = %d, k = %d): %10.6f\n",m,k,tt);
	if (array_norm_diff_d(m*n,C,C0,"Inf") > EPS*10)
		printf("Incorrect output, error: %e\n",array_norm_diff_d(m*n,C,C0,"Inf"));

	for (i = 0, iMax = m*n; iMax--; i++)
		C[i] = 0.0;
	ts = clock();
	for (iMax = NIter; iMax--; ) {
		mm_CTN_mv2_d(m,n,k,A,B,C);
	}
	te = clock();
	tt = (te-ts)/(double) CLOCKS_PER_SEC;

	printf("mm_CTN_mv2_d               : %10.6f\n",tt);
	if (array_norm_diff_d(m*n,C,C0,"Inf") > EPS*10)
		printf("Incorrect output, error: %e\n",array_norm_diff_d(m*n,C,C0,"Inf"));

	for (i = 0, iMax = m*n; iMax--; i++)
		C[i] = 0.0;
	ts = clock();
	for (iMax = NIter; iMax--; ) {
		mm_CTN_mv3_d(m,n,k,A,B,C);
	}
	te = clock();
	tt = (te-ts)/(double) CLOCKS_PER_SEC;

	printf("mm_CTN_mv3_d               : %10.6f\n",tt);
	if (array_norm_diff_d(m*n,C,C0,"Inf") > EPS*10)
		printf("Incorrect output, error: %e\n",array_norm_diff_d(m*n,C,C0,"Inf"));

	for (i = 0, iMax = m*n; iMax--; i++)
		C[i] = 0.0;
	ts = clock();
	for (iMax = NIter; iMax--; ) {
		mm_CTN_mv4_5d(m,n,k,A,B,C);
	}
	te = clock();
	tt = (te-ts)/(double) CLOCKS_PER_SEC;

	printf("mm_CTN_mv4_5d              : %10.6f\n",tt);
	if (array_norm_diff_d(m*n,C,C0,"Inf") > EPS*10)
		printf("Incorrect output, error: %e\n",array_norm_diff_d(m*n,C,C0,"Inf"));

	for (i = 0, iMax = m*n; iMax--; i++)
		C[i] = 0.0;
	ts = clock();
	for (iMax = NIter; iMax--; ) {
		mm_CTN_mv_unrolled_d(m,n,k,A,B,C);
	}
	te = clock();
	tt = (te-ts)/(double) CLOCKS_PER_SEC;

	printf("mm_CTN_mv_unrolled_d       : %10.6f\n",tt);
	if (array_norm_diff_d(m*n,C,C0,"Inf") > EPS*10)
		printf("Incorrect output, error: %e\n",array_norm_diff_d(m*n,C,C0,"Inf"));

	for (i = 0, iMax = m*n; iMax--; i++)
		C[i] = 0.0;
	ts = clock();
	for (iMax = NIter; iMax--; ) {
		mm_CTN_mv_fully_unrolled_d(m,n,k,A,B,C);
	}
	te = clock();
	tt = (te-ts)/(double) CLOCKS_PER_SEC;

	printf("mm_CTN_mv_fully_unrolled_d : %10.6f\n",tt);
	if (array_norm_diff_d(m*n,C,C0,"Inf") > EPS*10)
		printf("Incorrect output, error: %e\n",array_norm_diff_d(m*n,C,C0,"Inf"));

	for (i = 0, iMax = m*n; iMax--; i++)
		C[i] = 0.0;
	ts = clock();
	for (iMax = NIter; iMax--; ) {
		mm_CTN_fully_unrolled_mv_d(m,n,k,A,B,C);
	}
	te = clock();
	tt = (te-ts)/(double) CLOCKS_PER_SEC;

	printf("mm_CTN_fully_unrolled_mv_d : %10.6f\n",tt);
	if (array_norm_diff_d(m*n,C,C0,"Inf") > EPS*10)
		printf("Incorrect output, error: %e\n",array_norm_diff_d(m*n,C,C0,"Inf"));

	for (i = 0, iMax = m*n; iMax--; i++)
		C[i] = 0.0;
	ts = clock();
	for (iMax = NIter; iMax--; ) {
		mm_CTN_mvBLAS_d(m,n,k,A,B,C);
	}
	te = clock();
	tt = (te-ts)/(double) CLOCKS_PER_SEC;

	printf("mm_CTN_mvBLAS_d            : %10.6f\n",tt);
	if (array_norm_diff_d(m*n,C,C0,"Inf") > EPS*10)
		printf("Incorrect output, error: %e\n",array_norm_diff_d(m*n,C,C0,"Inf"));

	for (i = 0, iMax = m*n; iMax--; i++)
		C[i] = 0.0;
	ts = clock();
	for (iMax = NIter; iMax--; ) {
		mm_CTN_d(m,n,k,A,B,C);
	}
	te = clock();
	tt = (te-ts)/(double) CLOCKS_PER_SEC;

	printf("mm_CTN_d                   : %10.6f\n",tt);
	if (array_norm_diff_d(m*n,C,C0,"Inf") > EPS*10)
		printf("Incorrect output, error: %e\n",array_norm_diff_d(m*n,C,C0,"Inf"));

	free(A);
	free(B);
	free(C);
	free(C0);



	// (m)atrix-(v)ector: m = k = 8
	NIter = pow(2,16);
	m = 8;
	n = 1;
	k = m;
	flops = m*n*(2*k-1);

	A  = malloc(m*k * sizeof *A);
	B  = malloc(k*n * sizeof *B);
	C  = malloc(m*n * sizeof *C);
	C0 = malloc(m*n * sizeof *C0);

	for (i = 0, iMax = m*k; iMax--; i++)
		A[i] = (double)rand()/(double)RAND_MAX;

	for (i = 0, iMax = k*n; iMax--; i++)
		B[i] = (double)rand()/(double)RAND_MAX;

	mm_CTN_mvBLAS_d(m,n,k,A,B,C0);

	for (i = 0, iMax = m*n; iMax--; i++)
		C[i] = 0.0;
	ts = clock();
	for (iMax = NIter; iMax--; ) {
		mm_CTN_mv1_d(m,n,k,A,B,C);
	}
	te = clock();
	tt = (te-ts)/(double) CLOCKS_PER_SEC;

	printf("\nflops (%d) x NIter (%d) = Total flops (%d)\n",flops,NIter,flops*NIter);
	printf("mm_CTN_mv1_d (m = %d, k = %d): %10.6f\n",m,k,tt);
	if (array_norm_diff_d(m*n,C,C0,"Inf") > EPS*10)
		printf("Incorrect output, error: %e\n",array_norm_diff_d(m*n,C,C0,"Inf"));

	for (i = 0, iMax = m*n; iMax--; i++)
		C[i] = 0.0;
	ts = clock();
	for (iMax = NIter; iMax--; ) {
		mm_CTN_mv2_d(m,n,k,A,B,C);
	}
	te = clock();
	tt = (te-ts)/(double) CLOCKS_PER_SEC;

	printf("mm_CTN_mv2_d               : %10.6f\n",tt);
	if (array_norm_diff_d(m*n,C,C0,"Inf") > EPS*10)
		printf("Incorrect output, error: %e\n",array_norm_diff_d(m*n,C,C0,"Inf"));

	for (i = 0, iMax = m*n; iMax--; i++)
		C[i] = 0.0;
	ts = clock();
	for (iMax = NIter; iMax--; ) {
		mm_CTN_mv3_d(m,n,k,A,B,C);
	}
	te = clock();
	tt = (te-ts)/(double) CLOCKS_PER_SEC;

	printf("mm_CTN_mv3_d               : %10.6f\n",tt);
	if (array_norm_diff_d(m*n,C,C0,"Inf") > EPS*10)
		printf("Incorrect output, error: %e\n",array_norm_diff_d(m*n,C,C0,"Inf"));

	for (i = 0, iMax = m*n; iMax--; i++)
		C[i] = 0.0;
	ts = clock();
	for (iMax = NIter; iMax--; ) {
		mm_CTN_mv_unrolled_d(m,n,k,A,B,C);
	}
	te = clock();
	tt = (te-ts)/(double) CLOCKS_PER_SEC;

	printf("mm_CTN_mv_unrolled_d       : %10.6f\n",tt);
	if (array_norm_diff_d(m*n,C,C0,"Inf") > EPS*10)
		printf("Incorrect output, error: %e\n",array_norm_diff_d(m*n,C,C0,"Inf"));

	for (i = 0, iMax = m*n; iMax--; i++)
		C[i] = 0.0;
	ts = clock();
	for (iMax = NIter; iMax--; ) {
		mm_CTN_mv_fully_unrolled_d(m,n,k,A,B,C);
	}
	te = clock();
	tt = (te-ts)/(double) CLOCKS_PER_SEC;

	printf("mm_CTN_mv_fully_unrolled_d : %10.6f\n",tt);
	if (array_norm_diff_d(m*n,C,C0,"Inf") > EPS*10)
		printf("Incorrect output, error: %e\n",array_norm_diff_d(m*n,C,C0,"Inf"));

	for (i = 0, iMax = m*n; iMax--; i++)
		C[i] = 0.0;
	ts = clock();
	for (iMax = NIter; iMax--; ) {
		mm_CTN_fully_unrolled_mv_d(m,n,k,A,B,C);
	}
	te = clock();
	tt = (te-ts)/(double) CLOCKS_PER_SEC;

	printf("mm_CTN_fully_unrolled_mv_d : %10.6f\n",tt);
	if (array_norm_diff_d(m*n,C,C0,"Inf") > EPS*10)
		printf("Incorrect output, error: %e\n",array_norm_diff_d(m*n,C,C0,"Inf"));

	for (i = 0, iMax = m*n; iMax--; i++)
		C[i] = 0.0;
	ts = clock();
	for (iMax = NIter; iMax--; ) {
		mm_CTN_mvBLAS_d(m,n,k,A,B,C);
	}
	te = clock();
	tt = (te-ts)/(double) CLOCKS_PER_SEC;

	printf("mm_CTN_mvBLAS_d            : %10.6f\n",tt);
	if (array_norm_diff_d(m*n,C,C0,"Inf") > EPS*10)
		printf("Incorrect output, error: %e\n",array_norm_diff_d(m*n,C,C0,"Inf"));

	for (i = 0, iMax = m*n; iMax--; i++)
		C[i] = 0.0;
	ts = clock();
	for (iMax = NIter; iMax--; ) {
		mm_CTN_d(m,n,k,A,B,C);
	}
	te = clock();
	tt = (te-ts)/(double) CLOCKS_PER_SEC;

	printf("mm_CTN_d                   : %10.6f\n",tt);
	if (array_norm_diff_d(m*n,C,C0,"Inf") > EPS*10)
		printf("Incorrect output, error: %e\n",array_norm_diff_d(m*n,C,C0,"Inf"));

	free(A);
	free(B);
	free(C);
	free(C0);



	// (m)atrix-(m)atrix: m = k = 8, n = LARGE NUMBER (~Number of elements)
	NIter = pow(2,10);
	m = 8;
	n = 100;
	k = m;
	flops = m*n*(2*k-1);

	A  = malloc(m*k * sizeof *A);
	B  = malloc(k*n * sizeof *B);
	C  = malloc(m*n * sizeof *C);
	C0 = malloc(m*n * sizeof *C0);

	for (i = 0, iMax = m*k; iMax--; i++)
		A[i] = (double)rand()/(double)RAND_MAX;

	for (i = 0, iMax = k*n; iMax--; i++)
		B[i] = (double)rand()/(double)RAND_MAX;

	mm_CTN_d(m,n,k,A,B,C0);

	for (i = 0, iMax = m*n; iMax--; i++)
		C[i] = 0.0;
	ts = clock();
	for (iMax = NIter; iMax--; ) {
		mm_CTN_mm1_d(m,n,k,A,B,C);
	}
	te = clock();
	tt = (te-ts)/(double) CLOCKS_PER_SEC;

	printf("\nflops (%d) x NIter (%d) = Total flops (%d)\n",flops,NIter,flops*NIter);
	printf("mm_CTN_mm1_d (m = %d, k = %d, n = %d): %10.6f\n",m,k,n,tt);
	if (array_norm_diff_d(m*n,C,C0,"Inf") > EPS*10)
		printf("Incorrect output, error: %e\n",array_norm_diff_d(m*n,C,C0,"Inf"));

	for (i = 0, iMax = m*n; iMax--; i++)
		C[i] = 0.0;
	ts = clock();
	for (iMax = NIter; iMax--; ) {
		mm_CTN_mm2_d(m,n,k,A,B,C);
	}
	te = clock();
	tt = (te-ts)/(double) CLOCKS_PER_SEC;

	printf("mm_CTN_mm2_d                        : %10.6f\n",tt);
	if (array_norm_diff_d(m*n,C,C0,"Inf") > EPS*10)
		printf("Incorrect output, error: %e\n",array_norm_diff_d(m*n,C,C0,"Inf"));

	for (i = 0, iMax = m*n; iMax--; i++)
		C[i] = 0.0;
	ts = clock();
	for (iMax = NIter; iMax--; ) {
		mm_CTN_mm3_d(m,n,k,A,B,C);
	}
	te = clock();
	tt = (te-ts)/(double) CLOCKS_PER_SEC;

	printf("mm_CTN_mm3_d                        : %10.6f\n",tt);
	if (array_norm_diff_d(m*n,C,C0,"Inf") > EPS*10)
		printf("Incorrect output, error: %e\n",array_norm_diff_d(m*n,C,C0,"Inf"));

	for (i = 0, iMax = m*n; iMax--; i++)
		C[i] = 0.0;
	ts = clock();
	for (iMax = NIter; iMax--; ) {
		mm_CTN_fully_unrolled_mv_d(m,n,k,A,B,C);
	}
	te = clock();
	tt = (te-ts)/(double) CLOCKS_PER_SEC;

	printf("mm_CTN_fully_unrolled_mv_d          : %10.6f\n",tt);
	if (array_norm_diff_d(m*n,C,C0,"Inf") > EPS*10)
		printf("Incorrect output, error: %e\n",array_norm_diff_d(m*n,C,C0,"Inf"));


	for (i = 0, iMax = m*n; iMax--; i++)
		C[i] = 0.0;
	ts = clock();
	for (iMax = NIter; iMax--; ) {
		mm_CTN_mmBLAS_d(m,n,k,A,B,C);
	}
	te = clock();
	tt = (te-ts)/(double) CLOCKS_PER_SEC;

	printf("mm_CTN_mmBLAS_d                     : %10.6f\n",tt);
	if (array_norm_diff_d(m*n,C,C0,"Inf") > EPS*10)
		printf("Incorrect output, error: %e\n",array_norm_diff_d(m*n,C,C0,"Inf"));

	for (i = 0, iMax = m*n; iMax--; i++)
		C[i] = 0.0;
	ts = clock();
	for (iMax = NIter; iMax--; ) {
		mm_CTN_d(m,n,k,A,B,C);
	}
	te = clock();
	tt = (te-ts)/(double) CLOCKS_PER_SEC;

	printf("mm_CTN_d                            : %10.6f\n",tt);

	free(A);
	free(B);
	free(C);
	free(C0);

}
