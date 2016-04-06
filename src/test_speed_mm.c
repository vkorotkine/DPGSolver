#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "functions.h"

#include "mkl.h"

/*
 *	Purpose:
 *		Test runtimes for various matrix sizes used in mm_*.
 *
 *	Comments:
 *
 *		To run this function, modify the arguments of mm_CTN_d and uncomment calls below. This is necessary as the
 *		'useBLAS' parameter was subsequently included as part of the function.
 *
 *		As all of the significantly expensive matrix-matrix/matrix-vector operations will be done using a column-major
 *		orientation, only these tests are performed.
 *		Experiment with column-major vs. row-major layout for the operator matrix and check for any performance
 *		differences. (ToBeDeleted)
 *
 *		Important conclusions from time testing:
 *			No significant difference between using A*B (A in column-major layout) and A'*B (A in row-major layout)
 *			ToBeModified.
 *			BLAS_mv breaks even when more than ~120 flops are performed.
 *			BLAS_mm breaks even when more than ~330 flops are performed.
 *
 *	Notation:
 *
 *  References:
 *
 */

void test_speed_mm_d(void)
{
	printf("test_speed_mm_d:\n\n");

	register unsigned int i, iMax, unroll;
	unsigned int m, n, k, NIter, useBLAS;
	double *A, *B, *C;

	MKL_INT m_MKL, n_MKL, k_MKL, ldA, ldB, ldC;

	clock_t ts, te;
	double tt;

	NIter = pow(2,18);
	//NIter = pow(2,0);
	m = 6;
	//n = 2;
	n = m-1;
	k = m;
	//k = 4;
	unroll = 4;

	A = malloc(m*k * sizeof *A);
	B = malloc(k*n * sizeof *B);
	C = malloc(m*n * sizeof *C);

	for (i = 0, iMax = m*k; iMax--; i++)
		A[i] = (double)rand()/(double)RAND_MAX;

	for (i = 0, iMax = k*n; iMax--; i++)
		B[i] = (double)rand()/(double)RAND_MAX;

//	for (i = 0, iMax = m*n; iMax--; i++)
//		C[i] = (double)rand()/(double)RAND_MAX;

	// Note: From intel's site the first call causes some setup and is slower

	// C = A'*B (A: row-major, B: col-major)
	mm_d(CblasColMajor,CblasTrans,CblasNoTrans,m,n,k,1.0,A,B,C);
	ts = clock();
	for (iMax = NIter; iMax--; ) {
		mm_d(CblasColMajor,CblasTrans,CblasNoTrans,m,n,k,1.0,A,B,C);
	}
	te = clock();
	tt = (te-ts)/(double) CLOCKS_PER_SEC;

	printf("  Time to run %d mm_d     calls of ColMajor A'*B using m = %d, n = %d, k = %d: %10.6f\n",NIter,m,n,k,tt);

	// C = A*B (A: col-major, B: col-major)
/*	 No significant difference with C = A'*B (A: row-major, B: col-major)

	mm_d(CblasColMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,A,B,C);
	ts = clock();
	for (iMax = NIter; iMax--; ) {
		mm_d(CblasColMajor,CblasNoTrans,CblasNoTrans,m,n,k,1.0,A,B,C);
	}
	te = clock();
	tt = (te-ts)/(double) CLOCKS_PER_SEC;

	printf("  Time to run %d mm_d calls of ColMajor A*B using m = %d, n = %d, k = %d:  %10.6f\n",NIter,m,n,k,tt);
*/

	// C = A'*B (A: col-major, B: col-major (*gemm))
/*	Potentially, very small difference between this and calling mm_d (resulting in 2 function calls per iteration) ~5%.
	m_MKL = m, n_MKL = n, k_MKL = k;
	ldA = k, ldB = k, ldC = m;

	cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,m_MKL,n_MKL,k_MKL,1.0,A,ldA,B,ldB,0.0,C,ldC);
	ts = clock();
	for (iMax = NIter; iMax--; ) {
		cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,m_MKL,n_MKL,k_MKL,1.0,A,ldA,B,ldB,0.0,C,ldC);
	}
	te = clock();
	tt = (te-ts)/(double) CLOCKS_PER_SEC;

	printf("  Time to run %d dgemm calls of ColMajor A'*B using m = %d, n = %d, k = %d: %10.6f\n",NIter,m,n,k,tt);
*/

	// C = A'*B (A: row-major, B: col-major, loop unrolling)
/*	No significant difference when unrolling is added.
	mm_d(CblasColMajor,CblasTrans,CblasNoTrans,m,n,k,1.0,A,B,C);
	ts = clock();
	switch (unroll) {
	case 1:
		for (iMax = NIter; iMax--; ) {
			mm_d(CblasColMajor,CblasTrans,CblasNoTrans,m,n,k,1.0,A,B,C);
		}
		break;
	case 2:
		for (iMax = NIter/2; iMax--; ) {
			mm_d(CblasColMajor,CblasTrans,CblasNoTrans,m,n,k,1.0,A,B,C);
			mm_d(CblasColMajor,CblasTrans,CblasNoTrans,m,n,k,1.0,A,B,C);
		}
		break;
	case 4:
		for (iMax = NIter/4; iMax--; ) {
			mm_d(CblasColMajor,CblasTrans,CblasNoTrans,m,n,k,1.0,A,B,C);
			mm_d(CblasColMajor,CblasTrans,CblasNoTrans,m,n,k,1.0,A,B,C);
			mm_d(CblasColMajor,CblasTrans,CblasNoTrans,m,n,k,1.0,A,B,C);
			mm_d(CblasColMajor,CblasTrans,CblasNoTrans,m,n,k,1.0,A,B,C);
		}
		break;
	default:
		printf("Error: Unsupported value of 'unroll' entered.\n"), exit(1);
		break;
	}
	if (NIter % unroll != 0)
		printf("Implement something to handle this.\n");

	te = clock();
	tt = (te-ts)/(double) CLOCKS_PER_SEC;

	printf("  Time to run %d mm_d  calls of ColMajor A'*B using m = %d, n = %d, k = %d: %10.6f (unroll = %d)\n",
	       NIter,m,n,k,tt,unroll);
*/

	// C = A'*B (A: row-major, B: col-major)
	useBLAS = 1;
//	mm_CTN_d(m,n,k,A,B,C,useBLAS);
	ts = clock();
	for (iMax = NIter; iMax--; ) {
//		mm_CTN_d(m,n,k,A,B,C,useBLAS);
	}
	te = clock();
	tt = (te-ts)/(double) CLOCKS_PER_SEC;

	printf("  Time to run %d mm_CTN_d calls of ColMajor A'*B using m = %d, n = %d, k = %d: %10.6f (useBLAS = %d)\n",
	       NIter,m,n,k,tt,useBLAS);

//array_print_d(m,n,C,'C');


	// C = A'*B (A: row-major, B: col-major (*gemv))
	if (n == 1) {
		// Pretending that A is in column-major storage and will be transposed => reverse dimensions
		// NOTE: DIFFERENT CONVENTION FROM cblas_*gemm WHERE DIMENSIONS ARE OF op(A); here dimensions are of A.
		MKL_INT m_MKL = k, n_MKL = m, incx = 1, incy = 1;

		cblas_dgemv(CblasColMajor,CblasTrans,m_MKL,n_MKL,1.0,A,m_MKL,B,incx,0.0,C,incy);
		ts = clock();
		for (iMax = NIter; iMax--; ) {
			cblas_dgemv(CblasColMajor,CblasTrans,m_MKL,n_MKL,1.0,A,m_MKL,B,incx,0.0,C,incy);
		}
		te = clock();
		tt = (te-ts)/(double) CLOCKS_PER_SEC;

		printf("  Time to run %d dgemv    calls of ColMajor A'*B using m = %d, n = %d, k = %d: %10.6f\n",NIter,m,n,k,tt);
	}

	// C = A'*B (A: row-major, B: col-major), No BLAS
	useBLAS = 0;
//	mm_CTN_d(m,n,k,A,B,C,useBLAS);
	ts = clock();
	for (iMax = NIter; iMax--; ) {
//		mm_CTN_d(m,n,k,A,B,C,useBLAS);
	}
	te = clock();
	tt = (te-ts)/(double) CLOCKS_PER_SEC;

	printf("  Time to run %d mm_CTN_d calls of ColMajor A'*B using m = %d, n = %d, k = %d: %10.6f (useBLAS = %d)\n",
	       NIter,m,n,k,tt,useBLAS);

//array_print_d(m,n,C,'C');


	// Try with some loop unrolling for the function calls. If this is effective, change the mm_d routine to put the
	// loop inside of the function and unroll using binary progression switch statements (Perhaps don't actually do this
	// until the code is profiled, but write it as a potential optimization strategy).

	// Write the mm_d routine using c code only.

	printf("\nFlops: %d \n",m*n*(2*k-1));

	free(A);
	free(B);
	free(C);
}
