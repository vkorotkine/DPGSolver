// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>

#include "mkl.h"

/*
 *	Purpose:
 *		Provide matrix functions (ToBeModified):
 *			double *identity_d(const unsigned int N);
 *			double *inverse_d(int N, int NRHS, double A, double b);
 *			double mm_Alloc_d(const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const int m, const int n,
 *			                  const int k, const double alpha, const double *A, const double *B)
 *			void mm_d(const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const int m, const int n,
 *			          const int k, const double alpha, const double *A, const double *B)
 *
 *	Comments:
 *		The functions using LAPACKE do not use const arguments for consistency with the LAPACKE function definitions.
 *
 *	Notation:
 *
 *	References:
 *
 */

double *diag_d(const double *x, const unsigned int N)
{
	/*
	 *	Comments:
	 *		The returned array requires an external free.
	 */

	unsigned int i, j;
	double *X;

	X = malloc(N*N * sizeof *X); // keep (requires external free)
	for (i = 0; i < N; i++) {
	for (j = 0; j < N; j++) {
		if (i == j) X[i*N+j] = x[i];
		else        X[i*N+j] = 0.0;
	}}
	return X;
}

double *identity_d(const unsigned int N)
{
	/*
	 *	Comments:
	 *		The returned array requires an external free.
	 */

	unsigned int i, j;
	double *I;

	I = malloc(N*N * sizeof *I); // keep (requires external free)
	for (i = 0; i < N; i++) {
	for (j = 0; j < N; j++) {
		if (i == j) I[i*N+j] = 1.0;
		else        I[i*N+j] = 0.0;
	}}
	return I;
}

double *inverse_d(const unsigned int N, const unsigned int NRHS, const double *AIn, const double *b)
{
	/*
	 *	Comments:
	 *		The LAPACKE_dgesv modifies A while solving for x (which overwrites b's memory). Thus, new memory is
	 *		allocated such that both AIn and b remain unchanged after solving in this routine.
	 *
	 *		The returned array requires an external free.
	 */

	unsigned int i, iMax;
	double *A, *x;

	lapack_int N_LA, NRHS_LA, ipiv[N], info;

	N_LA    = (lapack_int) N;
	NRHS_LA = (lapack_int) NRHS;

	A = malloc(N*N    * sizeof *A); // free
	x = malloc(N*NRHS * sizeof *x); // keep (requires external free)
	for (i = 0, iMax = N*NRHS; i < iMax; i++) {
		A[i] = AIn[i];
		x[i] = b[i];
	}

	info = LAPACKE_dgesv(LAPACK_ROW_MAJOR,N_LA,NRHS_LA,A,N_LA,ipiv,x,NRHS_LA);
	if (info > 0) {
		printf("The diagonal element of the triangular factor of A,\n");
		printf("U(%i,%i) is zero, so that A is singular;\n", info, info);
		printf("the solution could not be computed.\n");
		exit(1);
	}

	free(A);
	return x;
}

double *mm_Alloc_d(const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb,
                   const int m, const int n, const int k, const double alpha, const double *A, const double *B)
{
	/*
	 *	Purpose:
	 *		Returns: C = alpha*op(A)*op(B)
	 *
	 *	Comments:
	 *		The returned array requires an external free.
	 *
	 *	Notation:
	 *		m : Number of rows of matrix op(A)
	 *		n : Number of columns of matrix op(B)
	 *		k : Number of columns of matrix op(A)
	 *
	 *		op() : Transpose, Conjugate Transpose
	 *
	 *	Options: (ToBeModified - update as usage is made)
	 *		transa/transb: CblasNoTrans, CblasTrans, CblasConjTrans
	 */

	double *C;
	MKL_INT m_MKL, n_MKL, k_MKL, ldA, ldB, ldC;

	m_MKL = (MKL_INT) m;
	n_MKL = (MKL_INT) n;
	k_MKL = (MKL_INT) k;

	C = malloc(m*n * sizeof *C); // keep (requires external free)
	if (layout == CblasColMajor) {
		if (transa == CblasNoTrans) ldA = m_MKL;
		else                        ldA = k_MKL;

		if (transb == CblasNoTrans) ldB = k_MKL;
		else                        ldB = n_MKL;

		ldC = m_MKL;
		cblas_dgemm(CblasColMajor,transa,transb,m_MKL,n_MKL,k_MKL,alpha,A,ldA,B,ldB,0.0,C,ldC);
	} else if (layout == CblasRowMajor) {
		if (transa == CblasNoTrans) ldA = k_MKL;
		else                        ldA = m_MKL;

		if (transb == CblasNoTrans) ldB = n_MKL;
		else                        ldB = k_MKL;

		ldC = n_MKL;
		cblas_dgemm(CblasRowMajor,transa,transb,m_MKL,n_MKL,k_MKL,alpha,A,ldA,B,ldB,0.0,C,ldC);
	} else {
		printf("Error: Invalid layout in mm_*.\n"), exit(1);
	}

	return C;
}

void mm_d(const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const int m,
          const int n, const int k, const double alpha, const double *A, const double *B, double *C)
{
	/*
	 *	Purpose:
	 *		Returns: C = alpha*op(A)*op(B) with memory already allocated in calling function.
	 *
	 *	Comments:
	 *		The 'C' array pointer is not declared 'const' as this prefix is discarded by the cblas_dgemm call.
	 *
	 *	Notation:
	 *		m : Number of rows of matrix op(A)
	 *		n : Number of columns of matrix op(B)
	 *		k : Number of columns of matrix op(A)
	 *
	 *		op() : Transpose, Conjugate Transpose
	 *
	 *	Options: (ToBeModified - update as usage is made)
	 *		transa/transb: CblasNoTrans, CblasTrans, CblasConjTrans
	 */

	MKL_INT m_MKL, n_MKL, k_MKL, ldA, ldB, ldC;

	m_MKL = (MKL_INT) m;
	n_MKL = (MKL_INT) n;
	k_MKL = (MKL_INT) k;

	if (layout == CblasColMajor) {
		if (transa == CblasNoTrans) ldA = m_MKL;
		else                        ldA = k_MKL;

		if (transb == CblasNoTrans) ldB = k_MKL;
		else                        ldB = n_MKL;

		ldC = m_MKL;
		cblas_dgemm(CblasColMajor,transa,transb,m_MKL,n_MKL,k_MKL,alpha,A,ldA,B,ldB,0.0,C,ldC);
	} else if (layout == CblasRowMajor) {
		if (transa == CblasNoTrans) ldA = k_MKL;
		else                        ldA = m_MKL;

		if (transb == CblasNoTrans) ldB = n_MKL;
		else                        ldB = k_MKL;

		ldC = n_MKL;
		cblas_dgemm(CblasRowMajor,transa,transb,m_MKL,n_MKL,k_MKL,alpha,A,ldA,B,ldB,0.0,C,ldC);
	} else {
		printf("Error: Invalid layout in mm_*.\n"), exit(1);
	}
}

void mm_CTN_d(const int m, const int n, const int k, double *A, double *B, double *C)
{
	/*
	 *	Purpose:
	 *		Compute matrix-matrix and matrix-vector products:
	 *			0) return C = A'*B
	 *			1) Row-major storage of A
	 *			2) Column-major storage of B and C
	 *
	 *	Comments:
	 *		The 'CTN' function name refers to:
	 *			C : Column-major output
	 *			T : Transpose of A
	 *			N : No transpose of B
	 *		This implementation/memory ordering is optimal to reduce memory stride when multiplying small A matrices
	 *		(relative to B) by B as in vectorized implementations of the code. (ToBeModified)
	 *
	 *		The function allows for the usage of either custom c-code or BLAS calls depending on the value of 'useBLAS'.
	 *		The 'A', 'B' and 'C' array pointers are not declared as 'const' as the qualifiers are discarded during the
	 *		initialization/assignement of pointer variables in the custom implementation.
	 *
	 *		Note that the cblas calls to gemv and gemm have slightly different conventions for specifying the dimensions
	 *		of A:
	 *			gemm: m_MKL = m, k_MKL = k provide the dimensions of op(A) = A' (m rows, k columns).
	 *			gemv: m_MKL = k, k_MKL = m provide the dimensions of A (k rows, m columns)
	 *		In the custom implementations of mv and mm below, the pointer version was very slightly faster than the
	 *		array indexed version (similarly to array_swap).
	 *
	 *		If the breakeven values are changed, ensure that test_imp_matrix_mm is updated. (ToBeModified)
	 *
	 *		Conclusions from test_speed_mm:
	 *			Testing done using -O3.
	 *			Without any manual loop unrolling (potentially to be done after profiling; ToBeModified):
	 *				BLAS_mv breaks even when more than ~120 flops are performed.
	 *				BLAS_mm breaks even when more than ~330 flops are performed.
	 *			With matrix-vector product unrolling:
	 *				The custom implementation is generally equivalent to or faster than BLAS, even with a huge number of
	 *				switch statements. INVESTIGATE FURTHER WHEN PROFILING THE CODE (ToBeDeleted)
	 *
	 *	Notation:
	 *		m : Number of rows of A'
	 *		n : Number of columns of B
	 *		k : Number of columns of A'
	 *
	 *	References:
	 *
	 */

	int useBLAS;

	// Note: Flops = m*n*(2*k-1) ~= 2*m*n*k
	switch(n) {
	case 1:
		if (m*k > 60) // 60 = (breakeven flops)/2
			useBLAS = 1;
		else
			useBLAS = 0;

		break;
	default:
		if (m*n*k > 165) // 165 = (breakeven flops)/2
			useBLAS = 1;
		else
			useBLAS = 0;

		break;
	}


	switch (useBLAS) {
	case 0:
		switch (n) {
			case 1: { // matrix-vector
				register unsigned int mMax, kMax;
				register double *pA = A, *pB = B, *pC = C;

				// First row of A
				*pC = (*pA)*(*pB);
				for (kMax = k-1; kMax--; ) { // loop over columns of A/rows of B/C
					pA++;
					pB++;

					*pC += (*pA)*(*pB);
				}

				// Remaining rows of A
				for (mMax = m-1; mMax--; ) { // loop over rows of A
					pA++;
					pC++;
					pB = B;

					*pC = (*pA)*(*pB);
					for (kMax = k-1; kMax--; ) { // loop over columns of A/rows of B/C
						pA++;
						pB++;

						*pC += (*pA)*(*pB);
					}
				}

/* Very slightly slower than pointer version above (1-2%)

				register unsigned int i, j, iMax, jMax, IndA;
				for (i = 0, iMax = m; i < iMax; i++) {
					C[i] = 0.0;
					IndA = i*k;
					for (j = 0, jMax = k; j < jMax; j++) {
						C[i] += A[IndA+j]*B[j];
					}
				}
*/
				break;
			}
			default: { // matrix-matrix
				register unsigned int mMax, nMax, kMax;
				register double *pA = A, *pB = B, *pC = C;

				// First column of B/C
				// First row of A
				*pC = (*pA)*(*pB);
				for (kMax = k-1; kMax--; ) { // loop over columns of A/rows of B/C
					pA++;
					pB++;

					*pC += (*pA)*(*pB);
				}

				// Remaining rows of A
				for (mMax = m-1; mMax--; ) { // loop over rows of A
					pA++;
					pC++;
					pB = B;

					*pC = (*pA)*(*pB);
					for (kMax = k-1; kMax--; ) { // loop over columns of A/rows of B/C
						pA++;
						pB++;

						*pC += (*pA)*(*pB);
					}
				}

				// Remaining columns of B/C
				for (nMax = n-1; nMax--; ) { // loop over columns of B/C
					pA = A;
					pB = B+(n-nMax-1)*k;
					pC++;
					*pC = (*pA)*(*pB);

					for (kMax = k-1; kMax--; ) { // loop over columns of A/rows of B/C
						pA++;
						pB++;

						*pC += (*pA)*(*pB);
					}

					// Remaining rows of A
					for (mMax = m-1; mMax--; ) { // loop over rows of A
						pA++;
						pC++;
						pB = B+(n-nMax-1)*k;

						*pC = (*pA)*(*pB);
						for (kMax = k-1; kMax--; ) { // loop over columns of A/rows of B/C
							pA++;
							pB++;

							*pC += (*pA)*(*pB);
						}
					}
				}

				break;
			}
		}
		break;
	case 1:
		switch (n) {
			case 1: { // matrix-vector
				MKL_INT m_MKL   = k,
				        k_MKL   = m,
				        inc_MKL = 1;

				cblas_dgemv(CblasColMajor,CblasTrans,m_MKL,k_MKL,1.0,A,m_MKL,B,inc_MKL,0.0,C,inc_MKL);
				break;
			}
			default: { // matrix-matrix
				MKL_INT m_MKL = m,
				        n_MKL = n,
				        k_MKL = k;

				cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,m_MKL,n_MKL,k_MKL,1.0,A,k_MKL,B,k_MKL,0.0,C,m_MKL);
				break;
			}
		}
		break;
	default:
		printf("Error: Unsupported value of 'useBLAS' passed to mm_d_CTN.\n"), exit(1);
		break;
	}
}
