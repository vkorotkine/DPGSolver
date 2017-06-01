// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "matrix_functions.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>

#include "mkl.h"

#include "Parameters.h"
#include "Macros.h"
#include "S_OpCSR.h"
#include "matrix_structs.h"

#include "array_print.h"

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
	unsigned int i, j;
	double *A;

	A = malloc(N*N * sizeof *A); // keep (requires external free)
	for (i = 0; i < N; i++) {
	for (j = 0; j < N; j++) {
		if (i == j) A[i*N+j] = 1.0;
		else        A[i*N+j] = 0.0;
	}}
	return A;
}

double *inverse_d(const unsigned int N, const unsigned int NRHS, const double *AIn, const double *b)
{
	/*
	 *	Comments:
	 *		The LAPACKE_dgesv modifies A while solving for x (which overwrites b's memory). Thus, new memory is
	 *		allocated such that both AIn and b remain unchanged after solving in this routine.
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
		printf("A diagonal element of the triangular factor of A is zero, so that A is singular.\n");
		EXIT_MSG;
	}

	free(A);
//	mkl_free_buffers(); // Potentially only works with 2017 version of IntelMKL (ToBeModified)
	return x;
}

void mm_diag_d(const unsigned int NRows, const unsigned int NCols, double const *const a, double const *const A,
               double *const Output, const double alpha, const double beta, const char side, const char layout)
{
	unsigned int i, j, iMax;

	if (beta != 1.0) {
		if (beta == 0.0) {
			for (i = 0, iMax = NRows*NCols; i < iMax; i++)
				Output[i] = 0.0;
		} else {
			for (i = 0, iMax = NRows*NCols; i < iMax; i++)
				Output[i] = beta*Output[i];
		}
	}

	if (layout == 'R') {
		if (side == 'L') {
			for (i = 0; i < NRows; i++) {
			for (j = 0; j < NCols; j++) {
				Output[i*NCols+j] += alpha*a[i]*A[i*NCols+j];
			}}
		} else if (side == 'R') {
			for (i = 0; i < NRows; i++) {
			for (j = 0; j < NCols; j++) {
				Output[i*NCols+j] += alpha*a[j]*A[i*NCols+j];
			}}
		} else {
			printf("Error: Unsupported.\n"), EXIT_MSG;
		}
	} else if (layout == 'C') {
		if (side == 'L') {
			for (j = 0; j < NCols; j++) {
			for (i = 0; i < NRows; i++) {
				Output[i+j*NRows] += alpha*a[i]*A[i+j*NRows];
			}}
		} else if (side == 'R') {
			for (j = 0; j < NCols; j++) {
			for (i = 0; i < NRows; i++) {
				Output[i+j*NRows] += alpha*a[j]*A[i+j*NRows];
			}}
		} else {
			printf("Error: Unsupported.\n"), EXIT_MSG;
		}
	}
}

double *mm_Alloc_d(const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb,
                   const int m, const int n, const int k, const double alpha, const double *A, const double *B)
{
	/*
	 *	Purpose:
	 *		Returns: C = alpha*op(A)*op(B)
	 *
	 *	Comments:
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

	if (m == 0 || n == 0 || k == 0)
		printf("Error: %d %d %d\n",m,n,k), EXIT_MSG;

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
		printf("Error: Invalid layout in mm_*.\n"), EXIT_MSG;
	}

	return C;
}

void mm_d(const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const int m,
          const int n, const int k, const double alpha, const double beta, const double *A, const double *B, double *C)
{
	/*
	 *	Purpose:
	 *		Returns: C = alpha*op(A)*op(B)+beta*C with memory already allocated in calling function.
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

	if (m == 0 || n == 0 || k == 0)
		printf("Error: %d %d %d\n",m,n,k), EXIT_MSG;

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
		cblas_dgemm(CblasColMajor,transa,transb,m_MKL,n_MKL,k_MKL,alpha,A,ldA,B,ldB,beta,C,ldC);
	} else if (layout == CblasRowMajor) {
		if (transa == CblasNoTrans) ldA = k_MKL;
		else                        ldA = m_MKL;

		if (transb == CblasNoTrans) ldB = n_MKL;
		else                        ldB = k_MKL;

		ldC = n_MKL;
		cblas_dgemm(CblasRowMajor,transa,transb,m_MKL,n_MKL,k_MKL,alpha,A,ldA,B,ldB,beta,C,ldC);
	} else {
		printf("Error: Invalid layout in mm_*.\n"), EXIT_MSG;
	}
}

void mm_dcc(const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const int m,
            const int n, const int k, const double alpha, const double beta, const double *A, const void *B, void *C)
{
	/*
	 *	Purpose:
	 *		Returns: C = alpha*op(A)*op(B)+beta*C with memory already allocated in calling function.
	 *
	 *	Comments:
	 *		The 'C' array pointer is not declared 'const' as this prefix is discarded by the cblas_zgemm call.
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

	if (m == 0 || n == 0 || k == 0)
		printf("Error: %d %d %d\n",m,n,k), EXIT_MSG;

	unsigned int  i, iMax;
	MKL_INT       m_MKL, n_MKL, k_MKL, ldA, ldB, ldC;
	MKL_Complex16 alpha_c, beta_c, *A_c, *C_c;
	const MKL_Complex16 *B_c;

	m_MKL = (MKL_INT) m;
	n_MKL = (MKL_INT) n;
	k_MKL = (MKL_INT) k;

	alpha_c.real = alpha;
	alpha_c.imag = 0.0;
	beta_c.real  = beta;
	beta_c.imag  = 0.0;

	A_c = malloc(m*k * sizeof *A_c); // free
	for (i = 0, iMax = m*k; i < iMax; i++) {
		A_c[i].real = A[i];
		A_c[i].imag = 0.0;
	}

	B_c = B;
	C_c = C;

	if (layout == CblasColMajor) {
		if (transa == CblasNoTrans) ldA = m_MKL;
		else                        ldA = k_MKL;

		if (transb == CblasNoTrans) ldB = k_MKL;
		else                        ldB = n_MKL;

		ldC = m_MKL;
		cblas_zgemm(CblasColMajor,transa,transb,m_MKL,n_MKL,k_MKL,&alpha_c,A_c,ldA,B_c,ldB,&beta_c,C_c,ldC);
	} else if (layout == CblasRowMajor) {
		if (transa == CblasNoTrans) ldA = k_MKL;
		else                        ldA = m_MKL;

		if (transb == CblasNoTrans) ldB = n_MKL;
		else                        ldB = k_MKL;

		ldC = n_MKL;
		cblas_zgemm(CblasRowMajor,transa,transb,m_MKL,n_MKL,k_MKL,&alpha_c,A_c,ldA,B_c,ldB,&beta_c,C_c,ldC);
	} else {
		printf("Error: Invalid layout.\n"), EXIT_MSG;
	}

	free(A_c);
}

void mm_CTN_d(const int m, const int n, const int k, const double *A, const double *B, double *C)
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
	 *			With matrix-vector product unrolling:
	 *				The custom implementation is generally equivalent to or faster than BLAS, even with a huge number of
	 *				switch statements only on OSX. On Ubuntu, BLAS breaks-even with the custom implementation for
	 *				A_{3x3} and becomes significantly faster for larger A. In both cases, only between 30-40% of peak
	 *				flops are achieved. (Investigate further, ToBeModified)
	 *
	 *	Notation:
	 *		m : Number of rows of A'
	 *		n : Number of columns of B
	 *		k : Number of columns of A'
	 *
	 *	References:
	 *
	 */

	if (m == 0 || n == 0 || k == 0)
		printf("Error: %d %d %d\n",m,n,k), EXIT_MSG;

	unsigned int useBLAS;

	if (m > 8 || k > 8)
		useBLAS = 1;
	else
		useBLAS = 0;

	switch (useBLAS) {
	case 0:
		switch(m) {
		case 1: {
			switch(k) {
			case 1: {
				register const double *a0  = A   ,
				                *b0  = B   ;

				*C     = (*a0 )*(*b0 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k;

					*(++C) = (*a0 )*(*b0 );
				}
				break; // m1k1
			} case 2: {
				register const double *a0  = A   , *a1  = A+1 ,
				                *b0  = B   , *b1  = B+1 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 );
				}
				break; // m1k2
			} case 3: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 ,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 );
				}
				break; // m1k3
			} case 4: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 ,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 );
				}
				break; // m1k4
			} case 5: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 ,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k, b4  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 );
				}
				break; // m1k5
			} case 6: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 , *a5  = A+5 ,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 , *b5  = B+5 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k, b4  += k, b5  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 );
				}
				break; // m1k6
			} case 7: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 , *a5  = A+5 , *a6  = A+6 ,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 , *b5  = B+5 , *b6  = B+6 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 ) + (*a6 )*(*b6 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k, b4  += k, b5  += k, b6  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 ) + (*a6 )*(*b6 );
				}
				break; // m1k7
			} case 8: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 , *a5  = A+5 , *a6  = A+6 , *a7  = A+7 ,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 , *b5  = B+5 , *b6  = B+6 , *b7  = B+7 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 ) + (*a6 )*(*b6 ) + (*a7 )*(*b7 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k, b4  += k, b5  += k, b6  += k, b7  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 ) + (*a6 )*(*b6 ) + (*a7 )*(*b7 );
				}
				break; // m1k8
			} default: {
				printf("Error: Unsupported m = %d, k = %d, combination in mm_CTN_d (useBlas = 0).\n",m,k), EXIT_MSG;
				break; // m1
			}}
			break;
		} case 2: {
			switch(k) {
			case 1: {
				register const double *a0  = A   ,
				                *a1  = A+1 ,
				                *b0  = B   ;

				*C     = (*a0 )*(*b0 );
				*(++C) = (*a1 )*(*b0 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k;

					*(++C) = (*a0 )*(*b0 );
					*(++C) = (*a1 )*(*b0 );
				}
				break; // m2k1
			} case 2: {
				register const double *a0  = A   , *a1  = A+1 ,
				                *a2  = A+2 , *a3  = A+3 ,
				                *b0  = B   , *b1  = B+1 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 );
				*(++C) = (*a2 )*(*b0 ) + (*a3 )*(*b1 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 );
					*(++C) = (*a2 )*(*b0 ) + (*a3 )*(*b1 );
				}
				break; // m2k2
			} case 3: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 ,
				                *a3  = A+3 , *a4  = A+4 , *a5  = A+5 ,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 );
				*(++C) = (*a3 )*(*b0 ) + (*a4 )*(*b1 ) + (*a5 )*(*b2 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 );
					*(++C) = (*a3 )*(*b0 ) + (*a4 )*(*b1 ) + (*a5 )*(*b2 );
				}
				break; // m2k3
			} case 4: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 ,
				                *a4  = A+4 , *a5  = A+5 , *a6  = A+6 , *a7  = A+7 ,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 );
				*(++C) = (*a4 )*(*b0 ) + (*a5 )*(*b1 ) + (*a6 )*(*b2 ) + (*a7 )*(*b3 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 );
					*(++C) = (*a4 )*(*b0 ) + (*a5 )*(*b1 ) + (*a6 )*(*b2 ) + (*a7 )*(*b3 );
				}
				break; // m2k4
			} case 5: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 ,
				                *a5  = A+5 , *a6  = A+6 , *a7  = A+7 , *a8  = A+8 , *a9  = A+9 ,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 );
				*(++C) = (*a5 )*(*b0 ) + (*a6 )*(*b1 ) + (*a7 )*(*b2 ) + (*a8 )*(*b3 ) + (*a9 )*(*b4 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k, b4  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 );
					*(++C) = (*a5 )*(*b0 ) + (*a6 )*(*b1 ) + (*a7 )*(*b2 ) + (*a8 )*(*b3 ) + (*a9 )*(*b4 );
				}
				break; // m2k5
			} case 6: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 , *a5  = A+5 ,
				                *a6  = A+6 , *a7  = A+7 , *a8  = A+8 , *a9  = A+9 , *a10 = A+10, *a11 = A+11,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 , *b5  = B+5 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 );
				*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 ) + (*a8 )*(*b2 ) + (*a9 )*(*b3 ) + (*a10)*(*b4 ) + (*a11)*(*b5 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k, b4  += k, b5  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 );
					*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 ) + (*a8 )*(*b2 ) + (*a9 )*(*b3 ) + (*a10)*(*b4 ) + (*a11)*(*b5 );
				}
				break; // m2k6
			} case 7: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 , *a5  = A+5 , *a6  = A+6 ,
				                *a7  = A+7 , *a8  = A+8 , *a9  = A+9 , *a10 = A+10, *a11 = A+11, *a12 = A+12, *a13 = A+13,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 , *b5  = B+5 , *b6  = B+6 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 ) + (*a6 )*(*b6 );
				*(++C) = (*a7 )*(*b0 ) + (*a8 )*(*b1 ) + (*a9 )*(*b2 ) + (*a10)*(*b3 ) + (*a11)*(*b4 ) + (*a12)*(*b5 ) + (*a13)*(*b6 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k, b4  += k, b5  += k, b6  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 ) + (*a6 )*(*b6 );
					*(++C) = (*a7 )*(*b0 ) + (*a8 )*(*b1 ) + (*a9 )*(*b2 ) + (*a10)*(*b3 ) + (*a11)*(*b4 ) + (*a12)*(*b5 ) + (*a13)*(*b6 );
				}
				break; // m2k7
			} case 8: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 , *a5  = A+5 , *a6  = A+6 , *a7  = A+7 ,
				                *a8  = A+8 , *a9  = A+9 , *a10 = A+10, *a11 = A+11, *a12 = A+12, *a13 = A+13, *a14 = A+14, *a15 = A+15,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 , *b5  = B+5 , *b6  = B+6 , *b7  = B+7 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 ) + (*a6 )*(*b6 ) + (*a7 )*(*b7 );
				*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 ) + (*a10)*(*b2 ) + (*a11)*(*b3 ) + (*a12)*(*b4 ) + (*a13)*(*b5 ) + (*a14)*(*b6 ) + (*a15)*(*b7 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k, b4  += k, b5  += k, b6  += k, b7  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 ) + (*a6 )*(*b6 ) + (*a7 )*(*b7 );
					*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 ) + (*a10)*(*b2 ) + (*a11)*(*b3 ) + (*a12)*(*b4 ) + (*a13)*(*b5 ) + (*a14)*(*b6 ) + (*a15)*(*b7 );
				}
				break; // m2k8
			} default: {
				printf("Error: Unsupported m = %d, k = %d, combination in mm_CTN_d (useBlas = 0).\n",m,k), EXIT_MSG;
				break; // m2
			}}
			break;
		} case 3: {
			switch(k) {
			case 1: {
				register const double *a0  = A   ,
				                *a1  = A+1 ,
				                *a2  = A+2 ,
				                *b0  = B   ;

				*C     = (*a0 )*(*b0 );
				*(++C) = (*a1 )*(*b0 );
				*(++C) = (*a2 )*(*b0 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k;

					*(++C) = (*a0 )*(*b0 );
					*(++C) = (*a1 )*(*b0 );
					*(++C) = (*a2 )*(*b0 );
				}
				break; // m3k1
			} case 2: {
				register const double *a0  = A   , *a1  = A+1 ,
				                *a2  = A+2 , *a3  = A+3 ,
				                *a4  = A+4 , *a5  = A+5 ,
				                *b0  = B   , *b1  = B+1 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 );
				*(++C) = (*a2 )*(*b0 ) + (*a3 )*(*b1 );
				*(++C) = (*a4 )*(*b0 ) + (*a5 )*(*b1 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 );
					*(++C) = (*a2 )*(*b0 ) + (*a3 )*(*b1 );
					*(++C) = (*a4 )*(*b0 ) + (*a5 )*(*b1 );
				}
				break; // m3k2
			} case 3: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 ,
				                *a3  = A+3 , *a4  = A+4 , *a5  = A+5 ,
				                *a6  = A+6 , *a7  = A+7 , *a8  = A+8 ,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 );
				*(++C) = (*a3 )*(*b0 ) + (*a4 )*(*b1 ) + (*a5 )*(*b2 );
				*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 ) + (*a8 )*(*b2 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 );
					*(++C) = (*a3 )*(*b0 ) + (*a4 )*(*b1 ) + (*a5 )*(*b2 );
					*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 ) + (*a8 )*(*b2 );
				}
				break; // m3k3
			} case 4: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 ,
				                *a4  = A+4 , *a5  = A+5 , *a6  = A+6 , *a7  = A+7 ,
				                *a8  = A+8 , *a9  = A+9 , *a10 = A+10, *a11 = A+11,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 );
				*(++C) = (*a4 )*(*b0 ) + (*a5 )*(*b1 ) + (*a6 )*(*b2 ) + (*a7 )*(*b3 );
				*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 ) + (*a10)*(*b2 ) + (*a11)*(*b3 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 );
					*(++C) = (*a4 )*(*b0 ) + (*a5 )*(*b1 ) + (*a6 )*(*b2 ) + (*a7 )*(*b3 );
					*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 ) + (*a10)*(*b2 ) + (*a11)*(*b3 );
				}
				break; // m3k4
			} case 5: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 ,
				                *a5  = A+5 , *a6  = A+6 , *a7  = A+7 , *a8  = A+8 , *a9  = A+9 ,
				                *a10 = A+10, *a11 = A+11, *a12 = A+12, *a13 = A+13, *a14 = A+14,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 );
				*(++C) = (*a5 )*(*b0 ) + (*a6 )*(*b1 ) + (*a7 )*(*b2 ) + (*a8 )*(*b3 ) + (*a9 )*(*b4 );
				*(++C) = (*a10)*(*b0 ) + (*a11)*(*b1 ) + (*a12)*(*b2 ) + (*a13)*(*b3 ) + (*a14)*(*b4 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k, b4  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 );
					*(++C) = (*a5 )*(*b0 ) + (*a6 )*(*b1 ) + (*a7 )*(*b2 ) + (*a8 )*(*b3 ) + (*a9 )*(*b4 );
					*(++C) = (*a10)*(*b0 ) + (*a11)*(*b1 ) + (*a12)*(*b2 ) + (*a13)*(*b3 ) + (*a14)*(*b4 );
				}
				break; // m3k5
			} case 6: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 , *a5  = A+5 ,
				                *a6  = A+6 , *a7  = A+7 , *a8  = A+8 , *a9  = A+9 , *a10 = A+10, *a11 = A+11,
				                *a12 = A+12, *a13 = A+13, *a14 = A+14, *a15 = A+15, *a16 = A+16, *a17 = A+17,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 , *b5  = B+5 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 );
				*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 ) + (*a8 )*(*b2 ) + (*a9 )*(*b3 ) + (*a10)*(*b4 ) + (*a11)*(*b5 );
				*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 ) + (*a14)*(*b2 ) + (*a15)*(*b3 ) + (*a16)*(*b4 ) + (*a17)*(*b5 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k, b4  += k, b5  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 );
					*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 ) + (*a8 )*(*b2 ) + (*a9 )*(*b3 ) + (*a10)*(*b4 ) + (*a11)*(*b5 );
					*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 ) + (*a14)*(*b2 ) + (*a15)*(*b3 ) + (*a16)*(*b4 ) + (*a17)*(*b5 );
				}
				break; // m3k6
			} case 7: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 , *a5  = A+5 , *a6  = A+6 ,
				                *a7  = A+7 , *a8  = A+8 , *a9  = A+9 , *a10 = A+10, *a11 = A+11, *a12 = A+12, *a13 = A+13,
				                *a14 = A+14, *a15 = A+15, *a16 = A+16, *a17 = A+17, *a18 = A+18, *a19 = A+19, *a20 = A+20,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 , *b5  = B+5 , *b6  = B+6 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 ) + (*a6 )*(*b6 );
				*(++C) = (*a7 )*(*b0 ) + (*a8 )*(*b1 ) + (*a9 )*(*b2 ) + (*a10)*(*b3 ) + (*a11)*(*b4 ) + (*a12)*(*b5 ) + (*a13)*(*b6 );
				*(++C) = (*a14)*(*b0 ) + (*a15)*(*b1 ) + (*a16)*(*b2 ) + (*a17)*(*b3 ) + (*a18)*(*b4 ) + (*a19)*(*b5 ) + (*a20)*(*b6 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k, b4  += k, b5  += k, b6  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 ) + (*a6 )*(*b6 );
					*(++C) = (*a7 )*(*b0 ) + (*a8 )*(*b1 ) + (*a9 )*(*b2 ) + (*a10)*(*b3 ) + (*a11)*(*b4 ) + (*a12)*(*b5 ) + (*a13)*(*b6 );
					*(++C) = (*a14)*(*b0 ) + (*a15)*(*b1 ) + (*a16)*(*b2 ) + (*a17)*(*b3 ) + (*a18)*(*b4 ) + (*a19)*(*b5 ) + (*a20)*(*b6 );
				}
				break; // m3k7
			} case 8: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 , *a5  = A+5 , *a6  = A+6 , *a7  = A+7 ,
				                *a8  = A+8 , *a9  = A+9 , *a10 = A+10, *a11 = A+11, *a12 = A+12, *a13 = A+13, *a14 = A+14, *a15 = A+15,
				                *a16 = A+16, *a17 = A+17, *a18 = A+18, *a19 = A+19, *a20 = A+20, *a21 = A+21, *a22 = A+22, *a23 = A+23,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 , *b5  = B+5 , *b6  = B+6 , *b7  = B+7 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 ) + (*a6 )*(*b6 ) + (*a7 )*(*b7 );
				*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 ) + (*a10)*(*b2 ) + (*a11)*(*b3 ) + (*a12)*(*b4 ) + (*a13)*(*b5 ) + (*a14)*(*b6 ) + (*a15)*(*b7 );
				*(++C) = (*a16)*(*b0 ) + (*a17)*(*b1 ) + (*a18)*(*b2 ) + (*a19)*(*b3 ) + (*a20)*(*b4 ) + (*a21)*(*b5 ) + (*a22)*(*b6 ) + (*a23)*(*b7 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k, b4  += k, b5  += k, b6  += k, b7  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 ) + (*a6 )*(*b6 ) + (*a7 )*(*b7 );
					*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 ) + (*a10)*(*b2 ) + (*a11)*(*b3 ) + (*a12)*(*b4 ) + (*a13)*(*b5 ) + (*a14)*(*b6 ) + (*a15)*(*b7 );
					*(++C) = (*a16)*(*b0 ) + (*a17)*(*b1 ) + (*a18)*(*b2 ) + (*a19)*(*b3 ) + (*a20)*(*b4 ) + (*a21)*(*b5 ) + (*a22)*(*b6 ) + (*a23)*(*b7 );
				}
				break; // m3k8
			} default: {
				printf("Error: Unsupported m = %d, k = %d, combination in mm_CTN_d (useBlas = 0).\n",m,k), EXIT_MSG;
				break; // m3
			}}
			break;
		} case 4: {
			switch(k) {
			case 1: {
				register const double *a0  = A   ,
				                *a1  = A+1 ,
				                *a2  = A+2 ,
				                *a3  = A+3 ,
				                *b0  = B   ;

				*C     = (*a0 )*(*b0 );
				*(++C) = (*a1 )*(*b0 );
				*(++C) = (*a2 )*(*b0 );
				*(++C) = (*a3 )*(*b0 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k;

					*(++C) = (*a0 )*(*b0 );
					*(++C) = (*a1 )*(*b0 );
					*(++C) = (*a2 )*(*b0 );
					*(++C) = (*a3 )*(*b0 );
				}
				break; // m4k1
			} case 2: {
				register const double *a0  = A   , *a1  = A+1 ,
				                *a2  = A+2 , *a3  = A+3 ,
				                *a4  = A+4 , *a5  = A+5 ,
				                *a6  = A+6 , *a7  = A+7 ,
				                *b0  = B   , *b1  = B+1 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 );
				*(++C) = (*a2 )*(*b0 ) + (*a3 )*(*b1 );
				*(++C) = (*a4 )*(*b0 ) + (*a5 )*(*b1 );
				*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 );
					*(++C) = (*a2 )*(*b0 ) + (*a3 )*(*b1 );
					*(++C) = (*a4 )*(*b0 ) + (*a5 )*(*b1 );
					*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 );
				}
				break; // m4k2
			} case 3: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 ,
				                *a3  = A+3 , *a4  = A+4 , *a5  = A+5 ,
				                *a6  = A+6 , *a7  = A+7 , *a8  = A+8 ,
				                *a9  = A+9 , *a10 = A+10, *a11 = A+11,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 );
				*(++C) = (*a3 )*(*b0 ) + (*a4 )*(*b1 ) + (*a5 )*(*b2 );
				*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 ) + (*a8 )*(*b2 );
				*(++C) = (*a9 )*(*b0 ) + (*a10)*(*b1 ) + (*a11)*(*b2 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 );
					*(++C) = (*a3 )*(*b0 ) + (*a4 )*(*b1 ) + (*a5 )*(*b2 );
					*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 ) + (*a8 )*(*b2 );
					*(++C) = (*a9 )*(*b0 ) + (*a10)*(*b1 ) + (*a11)*(*b2 );
				}
				break; // m4k3
			} case 4: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 ,
				                *a4  = A+4 , *a5  = A+5 , *a6  = A+6 , *a7  = A+7 ,
				                *a8  = A+8 , *a9  = A+9 , *a10 = A+10, *a11 = A+11,
				                *a12 = A+12, *a13 = A+13, *a14 = A+14, *a15 = A+15,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 );
				*(++C) = (*a4 )*(*b0 ) + (*a5 )*(*b1 ) + (*a6 )*(*b2 ) + (*a7 )*(*b3 );
				*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 ) + (*a10)*(*b2 ) + (*a11)*(*b3 );
				*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 ) + (*a14)*(*b2 ) + (*a15)*(*b3 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 );
					*(++C) = (*a4 )*(*b0 ) + (*a5 )*(*b1 ) + (*a6 )*(*b2 ) + (*a7 )*(*b3 );
					*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 ) + (*a10)*(*b2 ) + (*a11)*(*b3 );
					*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 ) + (*a14)*(*b2 ) + (*a15)*(*b3 );
				}
				break; // m4k4
			} case 5: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 ,
				                *a5  = A+5 , *a6  = A+6 , *a7  = A+7 , *a8  = A+8 , *a9  = A+9 ,
				                *a10 = A+10, *a11 = A+11, *a12 = A+12, *a13 = A+13, *a14 = A+14,
				                *a15 = A+15, *a16 = A+16, *a17 = A+17, *a18 = A+18, *a19 = A+19,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 );
				*(++C) = (*a5 )*(*b0 ) + (*a6 )*(*b1 ) + (*a7 )*(*b2 ) + (*a8 )*(*b3 ) + (*a9 )*(*b4 );
				*(++C) = (*a10)*(*b0 ) + (*a11)*(*b1 ) + (*a12)*(*b2 ) + (*a13)*(*b3 ) + (*a14)*(*b4 );
				*(++C) = (*a15)*(*b0 ) + (*a16)*(*b1 ) + (*a17)*(*b2 ) + (*a18)*(*b3 ) + (*a19)*(*b4 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k, b4  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 );
					*(++C) = (*a5 )*(*b0 ) + (*a6 )*(*b1 ) + (*a7 )*(*b2 ) + (*a8 )*(*b3 ) + (*a9 )*(*b4 );
					*(++C) = (*a10)*(*b0 ) + (*a11)*(*b1 ) + (*a12)*(*b2 ) + (*a13)*(*b3 ) + (*a14)*(*b4 );
					*(++C) = (*a15)*(*b0 ) + (*a16)*(*b1 ) + (*a17)*(*b2 ) + (*a18)*(*b3 ) + (*a19)*(*b4 );
				}
				break; // m4k5
			} case 6: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 , *a5  = A+5 ,
				                *a6  = A+6 , *a7  = A+7 , *a8  = A+8 , *a9  = A+9 , *a10 = A+10, *a11 = A+11,
				                *a12 = A+12, *a13 = A+13, *a14 = A+14, *a15 = A+15, *a16 = A+16, *a17 = A+17,
				                *a18 = A+18, *a19 = A+19, *a20 = A+20, *a21 = A+21, *a22 = A+22, *a23 = A+23,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 , *b5  = B+5 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 );
				*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 ) + (*a8 )*(*b2 ) + (*a9 )*(*b3 ) + (*a10)*(*b4 ) + (*a11)*(*b5 );
				*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 ) + (*a14)*(*b2 ) + (*a15)*(*b3 ) + (*a16)*(*b4 ) + (*a17)*(*b5 );
				*(++C) = (*a18)*(*b0 ) + (*a19)*(*b1 ) + (*a20)*(*b2 ) + (*a21)*(*b3 ) + (*a22)*(*b4 ) + (*a23)*(*b5 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k, b4  += k, b5  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 );
					*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 ) + (*a8 )*(*b2 ) + (*a9 )*(*b3 ) + (*a10)*(*b4 ) + (*a11)*(*b5 );
					*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 ) + (*a14)*(*b2 ) + (*a15)*(*b3 ) + (*a16)*(*b4 ) + (*a17)*(*b5 );
					*(++C) = (*a18)*(*b0 ) + (*a19)*(*b1 ) + (*a20)*(*b2 ) + (*a21)*(*b3 ) + (*a22)*(*b4 ) + (*a23)*(*b5 );
				}
				break; // m4k6
			} case 7: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 , *a5  = A+5 , *a6  = A+6 ,
				                *a7  = A+7 , *a8  = A+8 , *a9  = A+9 , *a10 = A+10, *a11 = A+11, *a12 = A+12, *a13 = A+13,
				                *a14 = A+14, *a15 = A+15, *a16 = A+16, *a17 = A+17, *a18 = A+18, *a19 = A+19, *a20 = A+20,
				                *a21 = A+21, *a22 = A+22, *a23 = A+23, *a24 = A+24, *a25 = A+25, *a26 = A+26, *a27 = A+27,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 , *b5  = B+5 , *b6  = B+6 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 ) + (*a6 )*(*b6 );
				*(++C) = (*a7 )*(*b0 ) + (*a8 )*(*b1 ) + (*a9 )*(*b2 ) + (*a10)*(*b3 ) + (*a11)*(*b4 ) + (*a12)*(*b5 ) + (*a13)*(*b6 );
				*(++C) = (*a14)*(*b0 ) + (*a15)*(*b1 ) + (*a16)*(*b2 ) + (*a17)*(*b3 ) + (*a18)*(*b4 ) + (*a19)*(*b5 ) + (*a20)*(*b6 );
				*(++C) = (*a21)*(*b0 ) + (*a22)*(*b1 ) + (*a23)*(*b2 ) + (*a24)*(*b3 ) + (*a25)*(*b4 ) + (*a26)*(*b5 ) + (*a27)*(*b6 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k, b4  += k, b5  += k, b6  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 ) + (*a6 )*(*b6 );
					*(++C) = (*a7 )*(*b0 ) + (*a8 )*(*b1 ) + (*a9 )*(*b2 ) + (*a10)*(*b3 ) + (*a11)*(*b4 ) + (*a12)*(*b5 ) + (*a13)*(*b6 );
					*(++C) = (*a14)*(*b0 ) + (*a15)*(*b1 ) + (*a16)*(*b2 ) + (*a17)*(*b3 ) + (*a18)*(*b4 ) + (*a19)*(*b5 ) + (*a20)*(*b6 );
					*(++C) = (*a21)*(*b0 ) + (*a22)*(*b1 ) + (*a23)*(*b2 ) + (*a24)*(*b3 ) + (*a25)*(*b4 ) + (*a26)*(*b5 ) + (*a27)*(*b6 );
				}
				break; // m4k7
			} case 8: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 , *a5  = A+5 , *a6  = A+6 , *a7  = A+7 ,
				                *a8  = A+8 , *a9  = A+9 , *a10 = A+10, *a11 = A+11, *a12 = A+12, *a13 = A+13, *a14 = A+14, *a15 = A+15,
				                *a16 = A+16, *a17 = A+17, *a18 = A+18, *a19 = A+19, *a20 = A+20, *a21 = A+21, *a22 = A+22, *a23 = A+23,
				                *a24 = A+24, *a25 = A+25, *a26 = A+26, *a27 = A+27, *a28 = A+28, *a29 = A+29, *a30 = A+30, *a31 = A+31,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 , *b5  = B+5 , *b6  = B+6 , *b7  = B+7 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 ) + (*a6 )*(*b6 ) + (*a7 )*(*b7 );
				*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 ) + (*a10)*(*b2 ) + (*a11)*(*b3 ) + (*a12)*(*b4 ) + (*a13)*(*b5 ) + (*a14)*(*b6 ) + (*a15)*(*b7 );
				*(++C) = (*a16)*(*b0 ) + (*a17)*(*b1 ) + (*a18)*(*b2 ) + (*a19)*(*b3 ) + (*a20)*(*b4 ) + (*a21)*(*b5 ) + (*a22)*(*b6 ) + (*a23)*(*b7 );
				*(++C) = (*a24)*(*b0 ) + (*a25)*(*b1 ) + (*a26)*(*b2 ) + (*a27)*(*b3 ) + (*a28)*(*b4 ) + (*a29)*(*b5 ) + (*a30)*(*b6 ) + (*a31)*(*b7 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k, b4  += k, b5  += k, b6  += k, b7  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 ) + (*a6 )*(*b6 ) + (*a7 )*(*b7 );
					*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 ) + (*a10)*(*b2 ) + (*a11)*(*b3 ) + (*a12)*(*b4 ) + (*a13)*(*b5 ) + (*a14)*(*b6 ) + (*a15)*(*b7 );
					*(++C) = (*a16)*(*b0 ) + (*a17)*(*b1 ) + (*a18)*(*b2 ) + (*a19)*(*b3 ) + (*a20)*(*b4 ) + (*a21)*(*b5 ) + (*a22)*(*b6 ) + (*a23)*(*b7 );
					*(++C) = (*a24)*(*b0 ) + (*a25)*(*b1 ) + (*a26)*(*b2 ) + (*a27)*(*b3 ) + (*a28)*(*b4 ) + (*a29)*(*b5 ) + (*a30)*(*b6 ) + (*a31)*(*b7 );
				}
				break; // m4k8
			} default: {
				printf("Error: Unsupported m = %d, k = %d, combination in mm_CTN_d (useBlas = 0).\n",m,k), EXIT_MSG;
				break; // m4
			}}
			break;
		} case 5: {
			switch(k) {
			case 1: {
				register const double *a0  = A   ,
				                *a1  = A+1 ,
				                *a2  = A+2 ,
				                *a3  = A+3 ,
				                *a4  = A+4 ,
				                *b0  = B   ;

				*C     = (*a0 )*(*b0 );
				*(++C) = (*a1 )*(*b0 );
				*(++C) = (*a2 )*(*b0 );
				*(++C) = (*a3 )*(*b0 );
				*(++C) = (*a4 )*(*b0 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k;

					*(++C) = (*a0 )*(*b0 );
					*(++C) = (*a1 )*(*b0 );
					*(++C) = (*a2 )*(*b0 );
					*(++C) = (*a3 )*(*b0 );
					*(++C) = (*a4 )*(*b0 );
				}
				break; // m5k1
			} case 2: {
				register const double *a0  = A   , *a1  = A+1 ,
				                *a2  = A+2 , *a3  = A+3 ,
				                *a4  = A+4 , *a5  = A+5 ,
				                *a6  = A+6 , *a7  = A+7 ,
				                *a8  = A+8 , *a9  = A+9 ,
				                *b0  = B   , *b1  = B+1 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 );
				*(++C) = (*a2 )*(*b0 ) + (*a3 )*(*b1 );
				*(++C) = (*a4 )*(*b0 ) + (*a5 )*(*b1 );
				*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 );
				*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 );
					*(++C) = (*a2 )*(*b0 ) + (*a3 )*(*b1 );
					*(++C) = (*a4 )*(*b0 ) + (*a5 )*(*b1 );
					*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 );
					*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 );
				}
				break; // m5k2
			} case 3: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 ,
				                *a3  = A+3 , *a4  = A+4 , *a5  = A+5 ,
				                *a6  = A+6 , *a7  = A+7 , *a8  = A+8 ,
				                *a9  = A+9 , *a10 = A+10, *a11 = A+11,
				                *a12 = A+12, *a13 = A+13, *a14 = A+14,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 );
				*(++C) = (*a3 )*(*b0 ) + (*a4 )*(*b1 ) + (*a5 )*(*b2 );
				*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 ) + (*a8 )*(*b2 );
				*(++C) = (*a9 )*(*b0 ) + (*a10)*(*b1 ) + (*a11)*(*b2 );
				*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 ) + (*a14)*(*b2 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 );
					*(++C) = (*a3 )*(*b0 ) + (*a4 )*(*b1 ) + (*a5 )*(*b2 );
					*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 ) + (*a8 )*(*b2 );
					*(++C) = (*a9 )*(*b0 ) + (*a10)*(*b1 ) + (*a11)*(*b2 );
					*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 ) + (*a14)*(*b2 );
				}
				break; // m5k3
			} case 4: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 ,
				                *a4  = A+4 , *a5  = A+5 , *a6  = A+6 , *a7  = A+7 ,
				                *a8  = A+8 , *a9  = A+9 , *a10 = A+10, *a11 = A+11,
				                *a12 = A+12, *a13 = A+13, *a14 = A+14, *a15 = A+15,
				                *a16 = A+16, *a17 = A+17, *a18 = A+18, *a19 = A+19,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 );
				*(++C) = (*a4 )*(*b0 ) + (*a5 )*(*b1 ) + (*a6 )*(*b2 ) + (*a7 )*(*b3 );
				*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 ) + (*a10)*(*b2 ) + (*a11)*(*b3 );
				*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 ) + (*a14)*(*b2 ) + (*a15)*(*b3 );
				*(++C) = (*a16)*(*b0 ) + (*a17)*(*b1 ) + (*a18)*(*b2 ) + (*a19)*(*b3 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 );
					*(++C) = (*a4 )*(*b0 ) + (*a5 )*(*b1 ) + (*a6 )*(*b2 ) + (*a7 )*(*b3 );
					*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 ) + (*a10)*(*b2 ) + (*a11)*(*b3 );
					*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 ) + (*a14)*(*b2 ) + (*a15)*(*b3 );
					*(++C) = (*a16)*(*b0 ) + (*a17)*(*b1 ) + (*a18)*(*b2 ) + (*a19)*(*b3 );
				}
				break; // m5k4
			} case 5: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 ,
				                *a5  = A+5 , *a6  = A+6 , *a7  = A+7 , *a8  = A+8 , *a9  = A+9 ,
				                *a10 = A+10, *a11 = A+11, *a12 = A+12, *a13 = A+13, *a14 = A+14,
				                *a15 = A+15, *a16 = A+16, *a17 = A+17, *a18 = A+18, *a19 = A+19,
				                *a20 = A+20, *a21 = A+21, *a22 = A+22, *a23 = A+23, *a24 = A+24,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 );
				*(++C) = (*a5 )*(*b0 ) + (*a6 )*(*b1 ) + (*a7 )*(*b2 ) + (*a8 )*(*b3 ) + (*a9 )*(*b4 );
				*(++C) = (*a10)*(*b0 ) + (*a11)*(*b1 ) + (*a12)*(*b2 ) + (*a13)*(*b3 ) + (*a14)*(*b4 );
				*(++C) = (*a15)*(*b0 ) + (*a16)*(*b1 ) + (*a17)*(*b2 ) + (*a18)*(*b3 ) + (*a19)*(*b4 );
				*(++C) = (*a20)*(*b0 ) + (*a21)*(*b1 ) + (*a22)*(*b2 ) + (*a23)*(*b3 ) + (*a24)*(*b4 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k, b4  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 );
					*(++C) = (*a5 )*(*b0 ) + (*a6 )*(*b1 ) + (*a7 )*(*b2 ) + (*a8 )*(*b3 ) + (*a9 )*(*b4 );
					*(++C) = (*a10)*(*b0 ) + (*a11)*(*b1 ) + (*a12)*(*b2 ) + (*a13)*(*b3 ) + (*a14)*(*b4 );
					*(++C) = (*a15)*(*b0 ) + (*a16)*(*b1 ) + (*a17)*(*b2 ) + (*a18)*(*b3 ) + (*a19)*(*b4 );
					*(++C) = (*a20)*(*b0 ) + (*a21)*(*b1 ) + (*a22)*(*b2 ) + (*a23)*(*b3 ) + (*a24)*(*b4 );
				}
				break; // m5k5
			} case 6: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 , *a5  = A+5 ,
				                *a6  = A+6 , *a7  = A+7 , *a8  = A+8 , *a9  = A+9 , *a10 = A+10, *a11 = A+11,
				                *a12 = A+12, *a13 = A+13, *a14 = A+14, *a15 = A+15, *a16 = A+16, *a17 = A+17,
				                *a18 = A+18, *a19 = A+19, *a20 = A+20, *a21 = A+21, *a22 = A+22, *a23 = A+23,
				                *a24 = A+24, *a25 = A+25, *a26 = A+26, *a27 = A+27, *a28 = A+28, *a29 = A+29,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 , *b5  = B+5 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 );
				*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 ) + (*a8 )*(*b2 ) + (*a9 )*(*b3 ) + (*a10)*(*b4 ) + (*a11)*(*b5 );
				*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 ) + (*a14)*(*b2 ) + (*a15)*(*b3 ) + (*a16)*(*b4 ) + (*a17)*(*b5 );
				*(++C) = (*a18)*(*b0 ) + (*a19)*(*b1 ) + (*a20)*(*b2 ) + (*a21)*(*b3 ) + (*a22)*(*b4 ) + (*a23)*(*b5 );
				*(++C) = (*a24)*(*b0 ) + (*a25)*(*b1 ) + (*a26)*(*b2 ) + (*a27)*(*b3 ) + (*a28)*(*b4 ) + (*a29)*(*b5 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k, b4  += k, b5  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 );
					*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 ) + (*a8 )*(*b2 ) + (*a9 )*(*b3 ) + (*a10)*(*b4 ) + (*a11)*(*b5 );
					*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 ) + (*a14)*(*b2 ) + (*a15)*(*b3 ) + (*a16)*(*b4 ) + (*a17)*(*b5 );
					*(++C) = (*a18)*(*b0 ) + (*a19)*(*b1 ) + (*a20)*(*b2 ) + (*a21)*(*b3 ) + (*a22)*(*b4 ) + (*a23)*(*b5 );
					*(++C) = (*a24)*(*b0 ) + (*a25)*(*b1 ) + (*a26)*(*b2 ) + (*a27)*(*b3 ) + (*a28)*(*b4 ) + (*a29)*(*b5 );
				}
				break; // m5k6
			} case 7: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 , *a5  = A+5 , *a6  = A+6 ,
				                *a7  = A+7 , *a8  = A+8 , *a9  = A+9 , *a10 = A+10, *a11 = A+11, *a12 = A+12, *a13 = A+13,
				                *a14 = A+14, *a15 = A+15, *a16 = A+16, *a17 = A+17, *a18 = A+18, *a19 = A+19, *a20 = A+20,
				                *a21 = A+21, *a22 = A+22, *a23 = A+23, *a24 = A+24, *a25 = A+25, *a26 = A+26, *a27 = A+27,
				                *a28 = A+28, *a29 = A+29, *a30 = A+30, *a31 = A+31, *a32 = A+32, *a33 = A+33, *a34 = A+34,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 , *b5  = B+5 , *b6  = B+6 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 ) + (*a6 )*(*b6 );
				*(++C) = (*a7 )*(*b0 ) + (*a8 )*(*b1 ) + (*a9 )*(*b2 ) + (*a10)*(*b3 ) + (*a11)*(*b4 ) + (*a12)*(*b5 ) + (*a13)*(*b6 );
				*(++C) = (*a14)*(*b0 ) + (*a15)*(*b1 ) + (*a16)*(*b2 ) + (*a17)*(*b3 ) + (*a18)*(*b4 ) + (*a19)*(*b5 ) + (*a20)*(*b6 );
				*(++C) = (*a21)*(*b0 ) + (*a22)*(*b1 ) + (*a23)*(*b2 ) + (*a24)*(*b3 ) + (*a25)*(*b4 ) + (*a26)*(*b5 ) + (*a27)*(*b6 );
				*(++C) = (*a28)*(*b0 ) + (*a29)*(*b1 ) + (*a30)*(*b2 ) + (*a31)*(*b3 ) + (*a32)*(*b4 ) + (*a33)*(*b5 ) + (*a34)*(*b6 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k, b4  += k, b5  += k, b6  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 ) + (*a6 )*(*b6 );
					*(++C) = (*a7 )*(*b0 ) + (*a8 )*(*b1 ) + (*a9 )*(*b2 ) + (*a10)*(*b3 ) + (*a11)*(*b4 ) + (*a12)*(*b5 ) + (*a13)*(*b6 );
					*(++C) = (*a14)*(*b0 ) + (*a15)*(*b1 ) + (*a16)*(*b2 ) + (*a17)*(*b3 ) + (*a18)*(*b4 ) + (*a19)*(*b5 ) + (*a20)*(*b6 );
					*(++C) = (*a21)*(*b0 ) + (*a22)*(*b1 ) + (*a23)*(*b2 ) + (*a24)*(*b3 ) + (*a25)*(*b4 ) + (*a26)*(*b5 ) + (*a27)*(*b6 );
					*(++C) = (*a28)*(*b0 ) + (*a29)*(*b1 ) + (*a30)*(*b2 ) + (*a31)*(*b3 ) + (*a32)*(*b4 ) + (*a33)*(*b5 ) + (*a34)*(*b6 );
				}
				break; // m5k7
			} case 8: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 , *a5  = A+5 , *a6  = A+6 , *a7  = A+7 ,
				                *a8  = A+8 , *a9  = A+9 , *a10 = A+10, *a11 = A+11, *a12 = A+12, *a13 = A+13, *a14 = A+14, *a15 = A+15,
				                *a16 = A+16, *a17 = A+17, *a18 = A+18, *a19 = A+19, *a20 = A+20, *a21 = A+21, *a22 = A+22, *a23 = A+23,
				                *a24 = A+24, *a25 = A+25, *a26 = A+26, *a27 = A+27, *a28 = A+28, *a29 = A+29, *a30 = A+30, *a31 = A+31,
				                *a32 = A+32, *a33 = A+33, *a34 = A+34, *a35 = A+35, *a36 = A+36, *a37 = A+37, *a38 = A+38, *a39 = A+39,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 , *b5  = B+5 , *b6  = B+6 , *b7  = B+7 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 ) + (*a6 )*(*b6 ) + (*a7 )*(*b7 );
				*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 ) + (*a10)*(*b2 ) + (*a11)*(*b3 ) + (*a12)*(*b4 ) + (*a13)*(*b5 ) + (*a14)*(*b6 ) + (*a15)*(*b7 );
				*(++C) = (*a16)*(*b0 ) + (*a17)*(*b1 ) + (*a18)*(*b2 ) + (*a19)*(*b3 ) + (*a20)*(*b4 ) + (*a21)*(*b5 ) + (*a22)*(*b6 ) + (*a23)*(*b7 );
				*(++C) = (*a24)*(*b0 ) + (*a25)*(*b1 ) + (*a26)*(*b2 ) + (*a27)*(*b3 ) + (*a28)*(*b4 ) + (*a29)*(*b5 ) + (*a30)*(*b6 ) + (*a31)*(*b7 );
				*(++C) = (*a32)*(*b0 ) + (*a33)*(*b1 ) + (*a34)*(*b2 ) + (*a35)*(*b3 ) + (*a36)*(*b4 ) + (*a37)*(*b5 ) + (*a38)*(*b6 ) + (*a39)*(*b7 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k, b4  += k, b5  += k, b6  += k, b7  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 ) + (*a6 )*(*b6 ) + (*a7 )*(*b7 );
					*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 ) + (*a10)*(*b2 ) + (*a11)*(*b3 ) + (*a12)*(*b4 ) + (*a13)*(*b5 ) + (*a14)*(*b6 ) + (*a15)*(*b7 );
					*(++C) = (*a16)*(*b0 ) + (*a17)*(*b1 ) + (*a18)*(*b2 ) + (*a19)*(*b3 ) + (*a20)*(*b4 ) + (*a21)*(*b5 ) + (*a22)*(*b6 ) + (*a23)*(*b7 );
					*(++C) = (*a24)*(*b0 ) + (*a25)*(*b1 ) + (*a26)*(*b2 ) + (*a27)*(*b3 ) + (*a28)*(*b4 ) + (*a29)*(*b5 ) + (*a30)*(*b6 ) + (*a31)*(*b7 );
					*(++C) = (*a32)*(*b0 ) + (*a33)*(*b1 ) + (*a34)*(*b2 ) + (*a35)*(*b3 ) + (*a36)*(*b4 ) + (*a37)*(*b5 ) + (*a38)*(*b6 ) + (*a39)*(*b7 );
				}
				break; // m5k8
			} default: {
				printf("Error: Unsupported m = %d, k = %d, combination in mm_CTN_d (useBlas = 0).\n",m,k), EXIT_MSG;
				break; // m5
			}}
			break;
		} case 6: {
			switch(k) {
			case 1: {
				register const double *a0  = A   ,
				                *a1  = A+1 ,
				                *a2  = A+2 ,
				                *a3  = A+3 ,
				                *a4  = A+4 ,
				                *a5  = A+5 ,
				                *b0  = B   ;

				*C     = (*a0 )*(*b0 );
				*(++C) = (*a1 )*(*b0 );
				*(++C) = (*a2 )*(*b0 );
				*(++C) = (*a3 )*(*b0 );
				*(++C) = (*a4 )*(*b0 );
				*(++C) = (*a5 )*(*b0 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k;

					*(++C) = (*a0 )*(*b0 );
					*(++C) = (*a1 )*(*b0 );
					*(++C) = (*a2 )*(*b0 );
					*(++C) = (*a3 )*(*b0 );
					*(++C) = (*a4 )*(*b0 );
					*(++C) = (*a5 )*(*b0 );
				}
				break; // m6k1
			} case 2: {
				register const double *a0  = A   , *a1  = A+1 ,
				                *a2  = A+2 , *a3  = A+3 ,
				                *a4  = A+4 , *a5  = A+5 ,
				                *a6  = A+6 , *a7  = A+7 ,
				                *a8  = A+8 , *a9  = A+9 ,
				                *a10 = A+10, *a11 = A+11,
				                *b0  = B   , *b1  = B+1 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 );
				*(++C) = (*a2 )*(*b0 ) + (*a3 )*(*b1 );
				*(++C) = (*a4 )*(*b0 ) + (*a5 )*(*b1 );
				*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 );
				*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 );
				*(++C) = (*a10)*(*b0 ) + (*a11)*(*b1 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 );
					*(++C) = (*a2 )*(*b0 ) + (*a3 )*(*b1 );
					*(++C) = (*a4 )*(*b0 ) + (*a5 )*(*b1 );
					*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 );
					*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 );
					*(++C) = (*a10)*(*b0 ) + (*a11)*(*b1 );
				}
				break; // m6k2
			} case 3: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 ,
				                *a3  = A+3 , *a4  = A+4 , *a5  = A+5 ,
				                *a6  = A+6 , *a7  = A+7 , *a8  = A+8 ,
				                *a9  = A+9 , *a10 = A+10, *a11 = A+11,
				                *a12 = A+12, *a13 = A+13, *a14 = A+14,
				                *a15 = A+15, *a16 = A+16, *a17 = A+17,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 );
				*(++C) = (*a3 )*(*b0 ) + (*a4 )*(*b1 ) + (*a5 )*(*b2 );
				*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 ) + (*a8 )*(*b2 );
				*(++C) = (*a9 )*(*b0 ) + (*a10)*(*b1 ) + (*a11)*(*b2 );
				*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 ) + (*a14)*(*b2 );
				*(++C) = (*a15)*(*b0 ) + (*a16)*(*b1 ) + (*a17)*(*b2 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 );
					*(++C) = (*a3 )*(*b0 ) + (*a4 )*(*b1 ) + (*a5 )*(*b2 );
					*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 ) + (*a8 )*(*b2 );
					*(++C) = (*a9 )*(*b0 ) + (*a10)*(*b1 ) + (*a11)*(*b2 );
					*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 ) + (*a14)*(*b2 );
					*(++C) = (*a15)*(*b0 ) + (*a16)*(*b1 ) + (*a17)*(*b2 );
				}
				break; // m6k3
			} case 4: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 ,
				                *a4  = A+4 , *a5  = A+5 , *a6  = A+6 , *a7  = A+7 ,
				                *a8  = A+8 , *a9  = A+9 , *a10 = A+10, *a11 = A+11,
				                *a12 = A+12, *a13 = A+13, *a14 = A+14, *a15 = A+15,
				                *a16 = A+16, *a17 = A+17, *a18 = A+18, *a19 = A+19,
				                *a20 = A+20, *a21 = A+21, *a22 = A+22, *a23 = A+23,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 );
				*(++C) = (*a4 )*(*b0 ) + (*a5 )*(*b1 ) + (*a6 )*(*b2 ) + (*a7 )*(*b3 );
				*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 ) + (*a10)*(*b2 ) + (*a11)*(*b3 );
				*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 ) + (*a14)*(*b2 ) + (*a15)*(*b3 );
				*(++C) = (*a16)*(*b0 ) + (*a17)*(*b1 ) + (*a18)*(*b2 ) + (*a19)*(*b3 );
				*(++C) = (*a20)*(*b0 ) + (*a21)*(*b1 ) + (*a22)*(*b2 ) + (*a23)*(*b3 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 );
					*(++C) = (*a4 )*(*b0 ) + (*a5 )*(*b1 ) + (*a6 )*(*b2 ) + (*a7 )*(*b3 );
					*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 ) + (*a10)*(*b2 ) + (*a11)*(*b3 );
					*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 ) + (*a14)*(*b2 ) + (*a15)*(*b3 );
					*(++C) = (*a16)*(*b0 ) + (*a17)*(*b1 ) + (*a18)*(*b2 ) + (*a19)*(*b3 );
					*(++C) = (*a20)*(*b0 ) + (*a21)*(*b1 ) + (*a22)*(*b2 ) + (*a23)*(*b3 );
				}
				break; // m6k4
			} case 5: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 ,
				                *a5  = A+5 , *a6  = A+6 , *a7  = A+7 , *a8  = A+8 , *a9  = A+9 ,
				                *a10 = A+10, *a11 = A+11, *a12 = A+12, *a13 = A+13, *a14 = A+14,
				                *a15 = A+15, *a16 = A+16, *a17 = A+17, *a18 = A+18, *a19 = A+19,
				                *a20 = A+20, *a21 = A+21, *a22 = A+22, *a23 = A+23, *a24 = A+24,
				                *a25 = A+25, *a26 = A+26, *a27 = A+27, *a28 = A+28, *a29 = A+29,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 );
				*(++C) = (*a5 )*(*b0 ) + (*a6 )*(*b1 ) + (*a7 )*(*b2 ) + (*a8 )*(*b3 ) + (*a9 )*(*b4 );
				*(++C) = (*a10)*(*b0 ) + (*a11)*(*b1 ) + (*a12)*(*b2 ) + (*a13)*(*b3 ) + (*a14)*(*b4 );
				*(++C) = (*a15)*(*b0 ) + (*a16)*(*b1 ) + (*a17)*(*b2 ) + (*a18)*(*b3 ) + (*a19)*(*b4 );
				*(++C) = (*a20)*(*b0 ) + (*a21)*(*b1 ) + (*a22)*(*b2 ) + (*a23)*(*b3 ) + (*a24)*(*b4 );
				*(++C) = (*a25)*(*b0 ) + (*a26)*(*b1 ) + (*a27)*(*b2 ) + (*a28)*(*b3 ) + (*a29)*(*b4 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k, b4  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 );
					*(++C) = (*a5 )*(*b0 ) + (*a6 )*(*b1 ) + (*a7 )*(*b2 ) + (*a8 )*(*b3 ) + (*a9 )*(*b4 );
					*(++C) = (*a10)*(*b0 ) + (*a11)*(*b1 ) + (*a12)*(*b2 ) + (*a13)*(*b3 ) + (*a14)*(*b4 );
					*(++C) = (*a15)*(*b0 ) + (*a16)*(*b1 ) + (*a17)*(*b2 ) + (*a18)*(*b3 ) + (*a19)*(*b4 );
					*(++C) = (*a20)*(*b0 ) + (*a21)*(*b1 ) + (*a22)*(*b2 ) + (*a23)*(*b3 ) + (*a24)*(*b4 );
					*(++C) = (*a25)*(*b0 ) + (*a26)*(*b1 ) + (*a27)*(*b2 ) + (*a28)*(*b3 ) + (*a29)*(*b4 );
				}
				break; // m6k5
			} case 6: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 , *a5  = A+5 ,
				                *a6  = A+6 , *a7  = A+7 , *a8  = A+8 , *a9  = A+9 , *a10 = A+10, *a11 = A+11,
				                *a12 = A+12, *a13 = A+13, *a14 = A+14, *a15 = A+15, *a16 = A+16, *a17 = A+17,
				                *a18 = A+18, *a19 = A+19, *a20 = A+20, *a21 = A+21, *a22 = A+22, *a23 = A+23,
				                *a24 = A+24, *a25 = A+25, *a26 = A+26, *a27 = A+27, *a28 = A+28, *a29 = A+29,
				                *a30 = A+30, *a31 = A+31, *a32 = A+32, *a33 = A+33, *a34 = A+34, *a35 = A+35,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 , *b5  = B+5 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 );
				*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 ) + (*a8 )*(*b2 ) + (*a9 )*(*b3 ) + (*a10)*(*b4 ) + (*a11)*(*b5 );
				*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 ) + (*a14)*(*b2 ) + (*a15)*(*b3 ) + (*a16)*(*b4 ) + (*a17)*(*b5 );
				*(++C) = (*a18)*(*b0 ) + (*a19)*(*b1 ) + (*a20)*(*b2 ) + (*a21)*(*b3 ) + (*a22)*(*b4 ) + (*a23)*(*b5 );
				*(++C) = (*a24)*(*b0 ) + (*a25)*(*b1 ) + (*a26)*(*b2 ) + (*a27)*(*b3 ) + (*a28)*(*b4 ) + (*a29)*(*b5 );
				*(++C) = (*a30)*(*b0 ) + (*a31)*(*b1 ) + (*a32)*(*b2 ) + (*a33)*(*b3 ) + (*a34)*(*b4 ) + (*a35)*(*b5 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k, b4  += k, b5  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 );
					*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 ) + (*a8 )*(*b2 ) + (*a9 )*(*b3 ) + (*a10)*(*b4 ) + (*a11)*(*b5 );
					*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 ) + (*a14)*(*b2 ) + (*a15)*(*b3 ) + (*a16)*(*b4 ) + (*a17)*(*b5 );
					*(++C) = (*a18)*(*b0 ) + (*a19)*(*b1 ) + (*a20)*(*b2 ) + (*a21)*(*b3 ) + (*a22)*(*b4 ) + (*a23)*(*b5 );
					*(++C) = (*a24)*(*b0 ) + (*a25)*(*b1 ) + (*a26)*(*b2 ) + (*a27)*(*b3 ) + (*a28)*(*b4 ) + (*a29)*(*b5 );
					*(++C) = (*a30)*(*b0 ) + (*a31)*(*b1 ) + (*a32)*(*b2 ) + (*a33)*(*b3 ) + (*a34)*(*b4 ) + (*a35)*(*b5 );
				}
				break; // m6k6
			} case 7: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 , *a5  = A+5 , *a6  = A+6 ,
				                *a7  = A+7 , *a8  = A+8 , *a9  = A+9 , *a10 = A+10, *a11 = A+11, *a12 = A+12, *a13 = A+13,
				                *a14 = A+14, *a15 = A+15, *a16 = A+16, *a17 = A+17, *a18 = A+18, *a19 = A+19, *a20 = A+20,
				                *a21 = A+21, *a22 = A+22, *a23 = A+23, *a24 = A+24, *a25 = A+25, *a26 = A+26, *a27 = A+27,
				                *a28 = A+28, *a29 = A+29, *a30 = A+30, *a31 = A+31, *a32 = A+32, *a33 = A+33, *a34 = A+34,
				                *a35 = A+35, *a36 = A+36, *a37 = A+37, *a38 = A+38, *a39 = A+39, *a40 = A+40, *a41 = A+41,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 , *b5  = B+5 , *b6  = B+6 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 ) + (*a6 )*(*b6 );
				*(++C) = (*a7 )*(*b0 ) + (*a8 )*(*b1 ) + (*a9 )*(*b2 ) + (*a10)*(*b3 ) + (*a11)*(*b4 ) + (*a12)*(*b5 ) + (*a13)*(*b6 );
				*(++C) = (*a14)*(*b0 ) + (*a15)*(*b1 ) + (*a16)*(*b2 ) + (*a17)*(*b3 ) + (*a18)*(*b4 ) + (*a19)*(*b5 ) + (*a20)*(*b6 );
				*(++C) = (*a21)*(*b0 ) + (*a22)*(*b1 ) + (*a23)*(*b2 ) + (*a24)*(*b3 ) + (*a25)*(*b4 ) + (*a26)*(*b5 ) + (*a27)*(*b6 );
				*(++C) = (*a28)*(*b0 ) + (*a29)*(*b1 ) + (*a30)*(*b2 ) + (*a31)*(*b3 ) + (*a32)*(*b4 ) + (*a33)*(*b5 ) + (*a34)*(*b6 );
				*(++C) = (*a35)*(*b0 ) + (*a36)*(*b1 ) + (*a37)*(*b2 ) + (*a38)*(*b3 ) + (*a39)*(*b4 ) + (*a40)*(*b5 ) + (*a41)*(*b6 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k, b4  += k, b5  += k, b6  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 ) + (*a6 )*(*b6 );
					*(++C) = (*a7 )*(*b0 ) + (*a8 )*(*b1 ) + (*a9 )*(*b2 ) + (*a10)*(*b3 ) + (*a11)*(*b4 ) + (*a12)*(*b5 ) + (*a13)*(*b6 );
					*(++C) = (*a14)*(*b0 ) + (*a15)*(*b1 ) + (*a16)*(*b2 ) + (*a17)*(*b3 ) + (*a18)*(*b4 ) + (*a19)*(*b5 ) + (*a20)*(*b6 );
					*(++C) = (*a21)*(*b0 ) + (*a22)*(*b1 ) + (*a23)*(*b2 ) + (*a24)*(*b3 ) + (*a25)*(*b4 ) + (*a26)*(*b5 ) + (*a27)*(*b6 );
					*(++C) = (*a28)*(*b0 ) + (*a29)*(*b1 ) + (*a30)*(*b2 ) + (*a31)*(*b3 ) + (*a32)*(*b4 ) + (*a33)*(*b5 ) + (*a34)*(*b6 );
					*(++C) = (*a35)*(*b0 ) + (*a36)*(*b1 ) + (*a37)*(*b2 ) + (*a38)*(*b3 ) + (*a39)*(*b4 ) + (*a40)*(*b5 ) + (*a41)*(*b6 );
				}
				break; // m6k7
			} case 8: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 , *a5  = A+5 , *a6  = A+6 , *a7  = A+7 ,
				                *a8  = A+8 , *a9  = A+9 , *a10 = A+10, *a11 = A+11, *a12 = A+12, *a13 = A+13, *a14 = A+14, *a15 = A+15,
				                *a16 = A+16, *a17 = A+17, *a18 = A+18, *a19 = A+19, *a20 = A+20, *a21 = A+21, *a22 = A+22, *a23 = A+23,
				                *a24 = A+24, *a25 = A+25, *a26 = A+26, *a27 = A+27, *a28 = A+28, *a29 = A+29, *a30 = A+30, *a31 = A+31,
				                *a32 = A+32, *a33 = A+33, *a34 = A+34, *a35 = A+35, *a36 = A+36, *a37 = A+37, *a38 = A+38, *a39 = A+39,
				                *a40 = A+40, *a41 = A+41, *a42 = A+42, *a43 = A+43, *a44 = A+44, *a45 = A+45, *a46 = A+46, *a47 = A+47,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 , *b5  = B+5 , *b6  = B+6 , *b7  = B+7 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 ) + (*a6 )*(*b6 ) + (*a7 )*(*b7 );
				*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 ) + (*a10)*(*b2 ) + (*a11)*(*b3 ) + (*a12)*(*b4 ) + (*a13)*(*b5 ) + (*a14)*(*b6 ) + (*a15)*(*b7 );
				*(++C) = (*a16)*(*b0 ) + (*a17)*(*b1 ) + (*a18)*(*b2 ) + (*a19)*(*b3 ) + (*a20)*(*b4 ) + (*a21)*(*b5 ) + (*a22)*(*b6 ) + (*a23)*(*b7 );
				*(++C) = (*a24)*(*b0 ) + (*a25)*(*b1 ) + (*a26)*(*b2 ) + (*a27)*(*b3 ) + (*a28)*(*b4 ) + (*a29)*(*b5 ) + (*a30)*(*b6 ) + (*a31)*(*b7 );
				*(++C) = (*a32)*(*b0 ) + (*a33)*(*b1 ) + (*a34)*(*b2 ) + (*a35)*(*b3 ) + (*a36)*(*b4 ) + (*a37)*(*b5 ) + (*a38)*(*b6 ) + (*a39)*(*b7 );
				*(++C) = (*a40)*(*b0 ) + (*a41)*(*b1 ) + (*a42)*(*b2 ) + (*a43)*(*b3 ) + (*a44)*(*b4 ) + (*a45)*(*b5 ) + (*a46)*(*b6 ) + (*a47)*(*b7 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k, b4  += k, b5  += k, b6  += k, b7  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 ) + (*a6 )*(*b6 ) + (*a7 )*(*b7 );
					*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 ) + (*a10)*(*b2 ) + (*a11)*(*b3 ) + (*a12)*(*b4 ) + (*a13)*(*b5 ) + (*a14)*(*b6 ) + (*a15)*(*b7 );
					*(++C) = (*a16)*(*b0 ) + (*a17)*(*b1 ) + (*a18)*(*b2 ) + (*a19)*(*b3 ) + (*a20)*(*b4 ) + (*a21)*(*b5 ) + (*a22)*(*b6 ) + (*a23)*(*b7 );
					*(++C) = (*a24)*(*b0 ) + (*a25)*(*b1 ) + (*a26)*(*b2 ) + (*a27)*(*b3 ) + (*a28)*(*b4 ) + (*a29)*(*b5 ) + (*a30)*(*b6 ) + (*a31)*(*b7 );
					*(++C) = (*a32)*(*b0 ) + (*a33)*(*b1 ) + (*a34)*(*b2 ) + (*a35)*(*b3 ) + (*a36)*(*b4 ) + (*a37)*(*b5 ) + (*a38)*(*b6 ) + (*a39)*(*b7 );
					*(++C) = (*a40)*(*b0 ) + (*a41)*(*b1 ) + (*a42)*(*b2 ) + (*a43)*(*b3 ) + (*a44)*(*b4 ) + (*a45)*(*b5 ) + (*a46)*(*b6 ) + (*a47)*(*b7 );
				}
				break; // m6k8
			} default: {
				printf("Error: Unsupported m = %d, k = %d, combination in mm_CTN_d (useBlas = 0).\n",m,k), EXIT_MSG;
				break; // m6
			}}
			break;
		} case 7: {
			switch(k) {
			case 1: {
				register const double *a0  = A   ,
				                *a1  = A+1 ,
				                *a2  = A+2 ,
				                *a3  = A+3 ,
				                *a4  = A+4 ,
				                *a5  = A+5 ,
				                *a6  = A+6 ,
				                *b0  = B   ;

				*C     = (*a0 )*(*b0 );
				*(++C) = (*a1 )*(*b0 );
				*(++C) = (*a2 )*(*b0 );
				*(++C) = (*a3 )*(*b0 );
				*(++C) = (*a4 )*(*b0 );
				*(++C) = (*a5 )*(*b0 );
				*(++C) = (*a6 )*(*b0 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k;

					*(++C) = (*a0 )*(*b0 );
					*(++C) = (*a1 )*(*b0 );
					*(++C) = (*a2 )*(*b0 );
					*(++C) = (*a3 )*(*b0 );
					*(++C) = (*a4 )*(*b0 );
					*(++C) = (*a5 )*(*b0 );
					*(++C) = (*a6 )*(*b0 );
				}
				break; // m7k1
			} case 2: {
				register const double *a0  = A   , *a1  = A+1 ,
				                *a2  = A+2 , *a3  = A+3 ,
				                *a4  = A+4 , *a5  = A+5 ,
				                *a6  = A+6 , *a7  = A+7 ,
				                *a8  = A+8 , *a9  = A+9 ,
				                *a10 = A+10, *a11 = A+11,
				                *a12 = A+12, *a13 = A+13,
				                *b0  = B   , *b1  = B+1 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 );
				*(++C) = (*a2 )*(*b0 ) + (*a3 )*(*b1 );
				*(++C) = (*a4 )*(*b0 ) + (*a5 )*(*b1 );
				*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 );
				*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 );
				*(++C) = (*a10)*(*b0 ) + (*a11)*(*b1 );
				*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 );
					*(++C) = (*a2 )*(*b0 ) + (*a3 )*(*b1 );
					*(++C) = (*a4 )*(*b0 ) + (*a5 )*(*b1 );
					*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 );
					*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 );
					*(++C) = (*a10)*(*b0 ) + (*a11)*(*b1 );
					*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 );
				}
				break; // m7k2
			} case 3: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 ,
				                *a3  = A+3 , *a4  = A+4 , *a5  = A+5 ,
				                *a6  = A+6 , *a7  = A+7 , *a8  = A+8 ,
				                *a9  = A+9 , *a10 = A+10, *a11 = A+11,
				                *a12 = A+12, *a13 = A+13, *a14 = A+14,
				                *a15 = A+15, *a16 = A+16, *a17 = A+17,
				                *a18 = A+18, *a19 = A+19, *a20 = A+20,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 );
				*(++C) = (*a3 )*(*b0 ) + (*a4 )*(*b1 ) + (*a5 )*(*b2 );
				*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 ) + (*a8 )*(*b2 );
				*(++C) = (*a9 )*(*b0 ) + (*a10)*(*b1 ) + (*a11)*(*b2 );
				*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 ) + (*a14)*(*b2 );
				*(++C) = (*a15)*(*b0 ) + (*a16)*(*b1 ) + (*a17)*(*b2 );
				*(++C) = (*a18)*(*b0 ) + (*a19)*(*b1 ) + (*a20)*(*b2 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 );
					*(++C) = (*a3 )*(*b0 ) + (*a4 )*(*b1 ) + (*a5 )*(*b2 );
					*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 ) + (*a8 )*(*b2 );
					*(++C) = (*a9 )*(*b0 ) + (*a10)*(*b1 ) + (*a11)*(*b2 );
					*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 ) + (*a14)*(*b2 );
					*(++C) = (*a15)*(*b0 ) + (*a16)*(*b1 ) + (*a17)*(*b2 );
					*(++C) = (*a18)*(*b0 ) + (*a19)*(*b1 ) + (*a20)*(*b2 );
				}
				break; // m7k3
			} case 4: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 ,
				                *a4  = A+4 , *a5  = A+5 , *a6  = A+6 , *a7  = A+7 ,
				                *a8  = A+8 , *a9  = A+9 , *a10 = A+10, *a11 = A+11,
				                *a12 = A+12, *a13 = A+13, *a14 = A+14, *a15 = A+15,
				                *a16 = A+16, *a17 = A+17, *a18 = A+18, *a19 = A+19,
				                *a20 = A+20, *a21 = A+21, *a22 = A+22, *a23 = A+23,
				                *a24 = A+24, *a25 = A+25, *a26 = A+26, *a27 = A+27,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 );
				*(++C) = (*a4 )*(*b0 ) + (*a5 )*(*b1 ) + (*a6 )*(*b2 ) + (*a7 )*(*b3 );
				*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 ) + (*a10)*(*b2 ) + (*a11)*(*b3 );
				*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 ) + (*a14)*(*b2 ) + (*a15)*(*b3 );
				*(++C) = (*a16)*(*b0 ) + (*a17)*(*b1 ) + (*a18)*(*b2 ) + (*a19)*(*b3 );
				*(++C) = (*a20)*(*b0 ) + (*a21)*(*b1 ) + (*a22)*(*b2 ) + (*a23)*(*b3 );
				*(++C) = (*a24)*(*b0 ) + (*a25)*(*b1 ) + (*a26)*(*b2 ) + (*a27)*(*b3 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 );
					*(++C) = (*a4 )*(*b0 ) + (*a5 )*(*b1 ) + (*a6 )*(*b2 ) + (*a7 )*(*b3 );
					*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 ) + (*a10)*(*b2 ) + (*a11)*(*b3 );
					*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 ) + (*a14)*(*b2 ) + (*a15)*(*b3 );
					*(++C) = (*a16)*(*b0 ) + (*a17)*(*b1 ) + (*a18)*(*b2 ) + (*a19)*(*b3 );
					*(++C) = (*a20)*(*b0 ) + (*a21)*(*b1 ) + (*a22)*(*b2 ) + (*a23)*(*b3 );
					*(++C) = (*a24)*(*b0 ) + (*a25)*(*b1 ) + (*a26)*(*b2 ) + (*a27)*(*b3 );
				}
				break; // m7k4
			} case 5: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 ,
				                *a5  = A+5 , *a6  = A+6 , *a7  = A+7 , *a8  = A+8 , *a9  = A+9 ,
				                *a10 = A+10, *a11 = A+11, *a12 = A+12, *a13 = A+13, *a14 = A+14,
				                *a15 = A+15, *a16 = A+16, *a17 = A+17, *a18 = A+18, *a19 = A+19,
				                *a20 = A+20, *a21 = A+21, *a22 = A+22, *a23 = A+23, *a24 = A+24,
				                *a25 = A+25, *a26 = A+26, *a27 = A+27, *a28 = A+28, *a29 = A+29,
				                *a30 = A+30, *a31 = A+31, *a32 = A+32, *a33 = A+33, *a34 = A+34,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 );
				*(++C) = (*a5 )*(*b0 ) + (*a6 )*(*b1 ) + (*a7 )*(*b2 ) + (*a8 )*(*b3 ) + (*a9 )*(*b4 );
				*(++C) = (*a10)*(*b0 ) + (*a11)*(*b1 ) + (*a12)*(*b2 ) + (*a13)*(*b3 ) + (*a14)*(*b4 );
				*(++C) = (*a15)*(*b0 ) + (*a16)*(*b1 ) + (*a17)*(*b2 ) + (*a18)*(*b3 ) + (*a19)*(*b4 );
				*(++C) = (*a20)*(*b0 ) + (*a21)*(*b1 ) + (*a22)*(*b2 ) + (*a23)*(*b3 ) + (*a24)*(*b4 );
				*(++C) = (*a25)*(*b0 ) + (*a26)*(*b1 ) + (*a27)*(*b2 ) + (*a28)*(*b3 ) + (*a29)*(*b4 );
				*(++C) = (*a30)*(*b0 ) + (*a31)*(*b1 ) + (*a32)*(*b2 ) + (*a33)*(*b3 ) + (*a34)*(*b4 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k, b4  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 );
					*(++C) = (*a5 )*(*b0 ) + (*a6 )*(*b1 ) + (*a7 )*(*b2 ) + (*a8 )*(*b3 ) + (*a9 )*(*b4 );
					*(++C) = (*a10)*(*b0 ) + (*a11)*(*b1 ) + (*a12)*(*b2 ) + (*a13)*(*b3 ) + (*a14)*(*b4 );
					*(++C) = (*a15)*(*b0 ) + (*a16)*(*b1 ) + (*a17)*(*b2 ) + (*a18)*(*b3 ) + (*a19)*(*b4 );
					*(++C) = (*a20)*(*b0 ) + (*a21)*(*b1 ) + (*a22)*(*b2 ) + (*a23)*(*b3 ) + (*a24)*(*b4 );
					*(++C) = (*a25)*(*b0 ) + (*a26)*(*b1 ) + (*a27)*(*b2 ) + (*a28)*(*b3 ) + (*a29)*(*b4 );
					*(++C) = (*a30)*(*b0 ) + (*a31)*(*b1 ) + (*a32)*(*b2 ) + (*a33)*(*b3 ) + (*a34)*(*b4 );
				}
				break; // m7k5
			} case 6: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 , *a5  = A+5 ,
				                *a6  = A+6 , *a7  = A+7 , *a8  = A+8 , *a9  = A+9 , *a10 = A+10, *a11 = A+11,
				                *a12 = A+12, *a13 = A+13, *a14 = A+14, *a15 = A+15, *a16 = A+16, *a17 = A+17,
				                *a18 = A+18, *a19 = A+19, *a20 = A+20, *a21 = A+21, *a22 = A+22, *a23 = A+23,
				                *a24 = A+24, *a25 = A+25, *a26 = A+26, *a27 = A+27, *a28 = A+28, *a29 = A+29,
				                *a30 = A+30, *a31 = A+31, *a32 = A+32, *a33 = A+33, *a34 = A+34, *a35 = A+35,
				                *a36 = A+36, *a37 = A+37, *a38 = A+38, *a39 = A+39, *a40 = A+40, *a41 = A+41,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 , *b5  = B+5 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 );
				*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 ) + (*a8 )*(*b2 ) + (*a9 )*(*b3 ) + (*a10)*(*b4 ) + (*a11)*(*b5 );
				*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 ) + (*a14)*(*b2 ) + (*a15)*(*b3 ) + (*a16)*(*b4 ) + (*a17)*(*b5 );
				*(++C) = (*a18)*(*b0 ) + (*a19)*(*b1 ) + (*a20)*(*b2 ) + (*a21)*(*b3 ) + (*a22)*(*b4 ) + (*a23)*(*b5 );
				*(++C) = (*a24)*(*b0 ) + (*a25)*(*b1 ) + (*a26)*(*b2 ) + (*a27)*(*b3 ) + (*a28)*(*b4 ) + (*a29)*(*b5 );
				*(++C) = (*a30)*(*b0 ) + (*a31)*(*b1 ) + (*a32)*(*b2 ) + (*a33)*(*b3 ) + (*a34)*(*b4 ) + (*a35)*(*b5 );
				*(++C) = (*a36)*(*b0 ) + (*a37)*(*b1 ) + (*a38)*(*b2 ) + (*a39)*(*b3 ) + (*a40)*(*b4 ) + (*a41)*(*b5 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k, b4  += k, b5  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 );
					*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 ) + (*a8 )*(*b2 ) + (*a9 )*(*b3 ) + (*a10)*(*b4 ) + (*a11)*(*b5 );
					*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 ) + (*a14)*(*b2 ) + (*a15)*(*b3 ) + (*a16)*(*b4 ) + (*a17)*(*b5 );
					*(++C) = (*a18)*(*b0 ) + (*a19)*(*b1 ) + (*a20)*(*b2 ) + (*a21)*(*b3 ) + (*a22)*(*b4 ) + (*a23)*(*b5 );
					*(++C) = (*a24)*(*b0 ) + (*a25)*(*b1 ) + (*a26)*(*b2 ) + (*a27)*(*b3 ) + (*a28)*(*b4 ) + (*a29)*(*b5 );
					*(++C) = (*a30)*(*b0 ) + (*a31)*(*b1 ) + (*a32)*(*b2 ) + (*a33)*(*b3 ) + (*a34)*(*b4 ) + (*a35)*(*b5 );
					*(++C) = (*a36)*(*b0 ) + (*a37)*(*b1 ) + (*a38)*(*b2 ) + (*a39)*(*b3 ) + (*a40)*(*b4 ) + (*a41)*(*b5 );
				}
				break; // m7k6
			} case 7: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 , *a5  = A+5 , *a6  = A+6 ,
				                *a7  = A+7 , *a8  = A+8 , *a9  = A+9 , *a10 = A+10, *a11 = A+11, *a12 = A+12, *a13 = A+13,
				                *a14 = A+14, *a15 = A+15, *a16 = A+16, *a17 = A+17, *a18 = A+18, *a19 = A+19, *a20 = A+20,
				                *a21 = A+21, *a22 = A+22, *a23 = A+23, *a24 = A+24, *a25 = A+25, *a26 = A+26, *a27 = A+27,
				                *a28 = A+28, *a29 = A+29, *a30 = A+30, *a31 = A+31, *a32 = A+32, *a33 = A+33, *a34 = A+34,
				                *a35 = A+35, *a36 = A+36, *a37 = A+37, *a38 = A+38, *a39 = A+39, *a40 = A+40, *a41 = A+41,
				                *a42 = A+42, *a43 = A+43, *a44 = A+44, *a45 = A+45, *a46 = A+46, *a47 = A+47, *a48 = A+48,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 , *b5  = B+5 , *b6  = B+6 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 ) + (*a6 )*(*b6 );
				*(++C) = (*a7 )*(*b0 ) + (*a8 )*(*b1 ) + (*a9 )*(*b2 ) + (*a10)*(*b3 ) + (*a11)*(*b4 ) + (*a12)*(*b5 ) + (*a13)*(*b6 );
				*(++C) = (*a14)*(*b0 ) + (*a15)*(*b1 ) + (*a16)*(*b2 ) + (*a17)*(*b3 ) + (*a18)*(*b4 ) + (*a19)*(*b5 ) + (*a20)*(*b6 );
				*(++C) = (*a21)*(*b0 ) + (*a22)*(*b1 ) + (*a23)*(*b2 ) + (*a24)*(*b3 ) + (*a25)*(*b4 ) + (*a26)*(*b5 ) + (*a27)*(*b6 );
				*(++C) = (*a28)*(*b0 ) + (*a29)*(*b1 ) + (*a30)*(*b2 ) + (*a31)*(*b3 ) + (*a32)*(*b4 ) + (*a33)*(*b5 ) + (*a34)*(*b6 );
				*(++C) = (*a35)*(*b0 ) + (*a36)*(*b1 ) + (*a37)*(*b2 ) + (*a38)*(*b3 ) + (*a39)*(*b4 ) + (*a40)*(*b5 ) + (*a41)*(*b6 );
				*(++C) = (*a42)*(*b0 ) + (*a43)*(*b1 ) + (*a44)*(*b2 ) + (*a45)*(*b3 ) + (*a46)*(*b4 ) + (*a47)*(*b5 ) + (*a48)*(*b6 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k, b4  += k, b5  += k, b6  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 ) + (*a6 )*(*b6 );
					*(++C) = (*a7 )*(*b0 ) + (*a8 )*(*b1 ) + (*a9 )*(*b2 ) + (*a10)*(*b3 ) + (*a11)*(*b4 ) + (*a12)*(*b5 ) + (*a13)*(*b6 );
					*(++C) = (*a14)*(*b0 ) + (*a15)*(*b1 ) + (*a16)*(*b2 ) + (*a17)*(*b3 ) + (*a18)*(*b4 ) + (*a19)*(*b5 ) + (*a20)*(*b6 );
					*(++C) = (*a21)*(*b0 ) + (*a22)*(*b1 ) + (*a23)*(*b2 ) + (*a24)*(*b3 ) + (*a25)*(*b4 ) + (*a26)*(*b5 ) + (*a27)*(*b6 );
					*(++C) = (*a28)*(*b0 ) + (*a29)*(*b1 ) + (*a30)*(*b2 ) + (*a31)*(*b3 ) + (*a32)*(*b4 ) + (*a33)*(*b5 ) + (*a34)*(*b6 );
					*(++C) = (*a35)*(*b0 ) + (*a36)*(*b1 ) + (*a37)*(*b2 ) + (*a38)*(*b3 ) + (*a39)*(*b4 ) + (*a40)*(*b5 ) + (*a41)*(*b6 );
					*(++C) = (*a42)*(*b0 ) + (*a43)*(*b1 ) + (*a44)*(*b2 ) + (*a45)*(*b3 ) + (*a46)*(*b4 ) + (*a47)*(*b5 ) + (*a48)*(*b6 );
				}
				break; // m7k7
			} case 8: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 , *a5  = A+5 , *a6  = A+6 , *a7  = A+7 ,
				                *a8  = A+8 , *a9  = A+9 , *a10 = A+10, *a11 = A+11, *a12 = A+12, *a13 = A+13, *a14 = A+14, *a15 = A+15,
				                *a16 = A+16, *a17 = A+17, *a18 = A+18, *a19 = A+19, *a20 = A+20, *a21 = A+21, *a22 = A+22, *a23 = A+23,
				                *a24 = A+24, *a25 = A+25, *a26 = A+26, *a27 = A+27, *a28 = A+28, *a29 = A+29, *a30 = A+30, *a31 = A+31,
				                *a32 = A+32, *a33 = A+33, *a34 = A+34, *a35 = A+35, *a36 = A+36, *a37 = A+37, *a38 = A+38, *a39 = A+39,
				                *a40 = A+40, *a41 = A+41, *a42 = A+42, *a43 = A+43, *a44 = A+44, *a45 = A+45, *a46 = A+46, *a47 = A+47,
				                *a48 = A+48, *a49 = A+49, *a50 = A+50, *a51 = A+51, *a52 = A+52, *a53 = A+53, *a54 = A+54, *a55 = A+55,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 , *b5  = B+5 , *b6  = B+6 , *b7  = B+7 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 ) + (*a6 )*(*b6 ) + (*a7 )*(*b7 );
				*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 ) + (*a10)*(*b2 ) + (*a11)*(*b3 ) + (*a12)*(*b4 ) + (*a13)*(*b5 ) + (*a14)*(*b6 ) + (*a15)*(*b7 );
				*(++C) = (*a16)*(*b0 ) + (*a17)*(*b1 ) + (*a18)*(*b2 ) + (*a19)*(*b3 ) + (*a20)*(*b4 ) + (*a21)*(*b5 ) + (*a22)*(*b6 ) + (*a23)*(*b7 );
				*(++C) = (*a24)*(*b0 ) + (*a25)*(*b1 ) + (*a26)*(*b2 ) + (*a27)*(*b3 ) + (*a28)*(*b4 ) + (*a29)*(*b5 ) + (*a30)*(*b6 ) + (*a31)*(*b7 );
				*(++C) = (*a32)*(*b0 ) + (*a33)*(*b1 ) + (*a34)*(*b2 ) + (*a35)*(*b3 ) + (*a36)*(*b4 ) + (*a37)*(*b5 ) + (*a38)*(*b6 ) + (*a39)*(*b7 );
				*(++C) = (*a40)*(*b0 ) + (*a41)*(*b1 ) + (*a42)*(*b2 ) + (*a43)*(*b3 ) + (*a44)*(*b4 ) + (*a45)*(*b5 ) + (*a46)*(*b6 ) + (*a47)*(*b7 );
				*(++C) = (*a48)*(*b0 ) + (*a49)*(*b1 ) + (*a50)*(*b2 ) + (*a51)*(*b3 ) + (*a52)*(*b4 ) + (*a53)*(*b5 ) + (*a54)*(*b6 ) + (*a55)*(*b7 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k, b4  += k, b5  += k, b6  += k, b7  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 ) + (*a6 )*(*b6 ) + (*a7 )*(*b7 );
					*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 ) + (*a10)*(*b2 ) + (*a11)*(*b3 ) + (*a12)*(*b4 ) + (*a13)*(*b5 ) + (*a14)*(*b6 ) + (*a15)*(*b7 );
					*(++C) = (*a16)*(*b0 ) + (*a17)*(*b1 ) + (*a18)*(*b2 ) + (*a19)*(*b3 ) + (*a20)*(*b4 ) + (*a21)*(*b5 ) + (*a22)*(*b6 ) + (*a23)*(*b7 );
					*(++C) = (*a24)*(*b0 ) + (*a25)*(*b1 ) + (*a26)*(*b2 ) + (*a27)*(*b3 ) + (*a28)*(*b4 ) + (*a29)*(*b5 ) + (*a30)*(*b6 ) + (*a31)*(*b7 );
					*(++C) = (*a32)*(*b0 ) + (*a33)*(*b1 ) + (*a34)*(*b2 ) + (*a35)*(*b3 ) + (*a36)*(*b4 ) + (*a37)*(*b5 ) + (*a38)*(*b6 ) + (*a39)*(*b7 );
					*(++C) = (*a40)*(*b0 ) + (*a41)*(*b1 ) + (*a42)*(*b2 ) + (*a43)*(*b3 ) + (*a44)*(*b4 ) + (*a45)*(*b5 ) + (*a46)*(*b6 ) + (*a47)*(*b7 );
					*(++C) = (*a48)*(*b0 ) + (*a49)*(*b1 ) + (*a50)*(*b2 ) + (*a51)*(*b3 ) + (*a52)*(*b4 ) + (*a53)*(*b5 ) + (*a54)*(*b6 ) + (*a55)*(*b7 );
				}
				break; // m7k8
			} default: {
				printf("Error: Unsupported m = %d, k = %d, combination in mm_CTN_d (useBlas = 0).\n",m,k), EXIT_MSG;
				break; // m7
			}}
			break;
		} case 8: {
			switch(k) {
			case 1: {
				register const double *a0  = A   ,
				                *a1  = A+1 ,
				                *a2  = A+2 ,
				                *a3  = A+3 ,
				                *a4  = A+4 ,
				                *a5  = A+5 ,
				                *a6  = A+6 ,
				                *a7  = A+7 ,
				                *b0  = B   ;

				*C     = (*a0 )*(*b0 );
				*(++C) = (*a1 )*(*b0 );
				*(++C) = (*a2 )*(*b0 );
				*(++C) = (*a3 )*(*b0 );
				*(++C) = (*a4 )*(*b0 );
				*(++C) = (*a5 )*(*b0 );
				*(++C) = (*a6 )*(*b0 );
				*(++C) = (*a7 )*(*b0 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k;

					*(++C) = (*a0 )*(*b0 );
					*(++C) = (*a1 )*(*b0 );
					*(++C) = (*a2 )*(*b0 );
					*(++C) = (*a3 )*(*b0 );
					*(++C) = (*a4 )*(*b0 );
					*(++C) = (*a5 )*(*b0 );
					*(++C) = (*a6 )*(*b0 );
					*(++C) = (*a7 )*(*b0 );
				}
				break; // m8k1
			} case 2: {
				register const double *a0  = A   , *a1  = A+1 ,
				                *a2  = A+2 , *a3  = A+3 ,
				                *a4  = A+4 , *a5  = A+5 ,
				                *a6  = A+6 , *a7  = A+7 ,
				                *a8  = A+8 , *a9  = A+9 ,
				                *a10 = A+10, *a11 = A+11,
				                *a12 = A+12, *a13 = A+13,
				                *a14 = A+14, *a15 = A+15,
				                *b0  = B   , *b1  = B+1 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 );
				*(++C) = (*a2 )*(*b0 ) + (*a3 )*(*b1 );
				*(++C) = (*a4 )*(*b0 ) + (*a5 )*(*b1 );
				*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 );
				*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 );
				*(++C) = (*a10)*(*b0 ) + (*a11)*(*b1 );
				*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 );
				*(++C) = (*a14)*(*b0 ) + (*a15)*(*b1 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 );
					*(++C) = (*a2 )*(*b0 ) + (*a3 )*(*b1 );
					*(++C) = (*a4 )*(*b0 ) + (*a5 )*(*b1 );
					*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 );
					*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 );
					*(++C) = (*a10)*(*b0 ) + (*a11)*(*b1 );
					*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 );
					*(++C) = (*a14)*(*b0 ) + (*a15)*(*b1 );
				}
				break; // m8k2
			} case 3: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 ,
				                *a3  = A+3 , *a4  = A+4 , *a5  = A+5 ,
				                *a6  = A+6 , *a7  = A+7 , *a8  = A+8 ,
				                *a9  = A+9 , *a10 = A+10, *a11 = A+11,
				                *a12 = A+12, *a13 = A+13, *a14 = A+14,
				                *a15 = A+15, *a16 = A+16, *a17 = A+17,
				                *a18 = A+18, *a19 = A+19, *a20 = A+20,
				                *a21 = A+21, *a22 = A+22, *a23 = A+23,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 );
				*(++C) = (*a3 )*(*b0 ) + (*a4 )*(*b1 ) + (*a5 )*(*b2 );
				*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 ) + (*a8 )*(*b2 );
				*(++C) = (*a9 )*(*b0 ) + (*a10)*(*b1 ) + (*a11)*(*b2 );
				*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 ) + (*a14)*(*b2 );
				*(++C) = (*a15)*(*b0 ) + (*a16)*(*b1 ) + (*a17)*(*b2 );
				*(++C) = (*a18)*(*b0 ) + (*a19)*(*b1 ) + (*a20)*(*b2 );
				*(++C) = (*a21)*(*b0 ) + (*a22)*(*b1 ) + (*a23)*(*b2 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 );
					*(++C) = (*a3 )*(*b0 ) + (*a4 )*(*b1 ) + (*a5 )*(*b2 );
					*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 ) + (*a8 )*(*b2 );
					*(++C) = (*a9 )*(*b0 ) + (*a10)*(*b1 ) + (*a11)*(*b2 );
					*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 ) + (*a14)*(*b2 );
					*(++C) = (*a15)*(*b0 ) + (*a16)*(*b1 ) + (*a17)*(*b2 );
					*(++C) = (*a18)*(*b0 ) + (*a19)*(*b1 ) + (*a20)*(*b2 );
					*(++C) = (*a21)*(*b0 ) + (*a22)*(*b1 ) + (*a23)*(*b2 );
				}
				break; // m8k3
			} case 4: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 ,
				                *a4  = A+4 , *a5  = A+5 , *a6  = A+6 , *a7  = A+7 ,
				                *a8  = A+8 , *a9  = A+9 , *a10 = A+10, *a11 = A+11,
				                *a12 = A+12, *a13 = A+13, *a14 = A+14, *a15 = A+15,
				                *a16 = A+16, *a17 = A+17, *a18 = A+18, *a19 = A+19,
				                *a20 = A+20, *a21 = A+21, *a22 = A+22, *a23 = A+23,
				                *a24 = A+24, *a25 = A+25, *a26 = A+26, *a27 = A+27,
				                *a28 = A+28, *a29 = A+29, *a30 = A+30, *a31 = A+31,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 );
				*(++C) = (*a4 )*(*b0 ) + (*a5 )*(*b1 ) + (*a6 )*(*b2 ) + (*a7 )*(*b3 );
				*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 ) + (*a10)*(*b2 ) + (*a11)*(*b3 );
				*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 ) + (*a14)*(*b2 ) + (*a15)*(*b3 );
				*(++C) = (*a16)*(*b0 ) + (*a17)*(*b1 ) + (*a18)*(*b2 ) + (*a19)*(*b3 );
				*(++C) = (*a20)*(*b0 ) + (*a21)*(*b1 ) + (*a22)*(*b2 ) + (*a23)*(*b3 );
				*(++C) = (*a24)*(*b0 ) + (*a25)*(*b1 ) + (*a26)*(*b2 ) + (*a27)*(*b3 );
				*(++C) = (*a28)*(*b0 ) + (*a29)*(*b1 ) + (*a30)*(*b2 ) + (*a31)*(*b3 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 );
					*(++C) = (*a4 )*(*b0 ) + (*a5 )*(*b1 ) + (*a6 )*(*b2 ) + (*a7 )*(*b3 );
					*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 ) + (*a10)*(*b2 ) + (*a11)*(*b3 );
					*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 ) + (*a14)*(*b2 ) + (*a15)*(*b3 );
					*(++C) = (*a16)*(*b0 ) + (*a17)*(*b1 ) + (*a18)*(*b2 ) + (*a19)*(*b3 );
					*(++C) = (*a20)*(*b0 ) + (*a21)*(*b1 ) + (*a22)*(*b2 ) + (*a23)*(*b3 );
					*(++C) = (*a24)*(*b0 ) + (*a25)*(*b1 ) + (*a26)*(*b2 ) + (*a27)*(*b3 );
					*(++C) = (*a28)*(*b0 ) + (*a29)*(*b1 ) + (*a30)*(*b2 ) + (*a31)*(*b3 );
				}
				break; // m8k4
			} case 5: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 ,
				                *a5  = A+5 , *a6  = A+6 , *a7  = A+7 , *a8  = A+8 , *a9  = A+9 ,
				                *a10 = A+10, *a11 = A+11, *a12 = A+12, *a13 = A+13, *a14 = A+14,
				                *a15 = A+15, *a16 = A+16, *a17 = A+17, *a18 = A+18, *a19 = A+19,
				                *a20 = A+20, *a21 = A+21, *a22 = A+22, *a23 = A+23, *a24 = A+24,
				                *a25 = A+25, *a26 = A+26, *a27 = A+27, *a28 = A+28, *a29 = A+29,
				                *a30 = A+30, *a31 = A+31, *a32 = A+32, *a33 = A+33, *a34 = A+34,
				                *a35 = A+35, *a36 = A+36, *a37 = A+37, *a38 = A+38, *a39 = A+39,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 );
				*(++C) = (*a5 )*(*b0 ) + (*a6 )*(*b1 ) + (*a7 )*(*b2 ) + (*a8 )*(*b3 ) + (*a9 )*(*b4 );
				*(++C) = (*a10)*(*b0 ) + (*a11)*(*b1 ) + (*a12)*(*b2 ) + (*a13)*(*b3 ) + (*a14)*(*b4 );
				*(++C) = (*a15)*(*b0 ) + (*a16)*(*b1 ) + (*a17)*(*b2 ) + (*a18)*(*b3 ) + (*a19)*(*b4 );
				*(++C) = (*a20)*(*b0 ) + (*a21)*(*b1 ) + (*a22)*(*b2 ) + (*a23)*(*b3 ) + (*a24)*(*b4 );
				*(++C) = (*a25)*(*b0 ) + (*a26)*(*b1 ) + (*a27)*(*b2 ) + (*a28)*(*b3 ) + (*a29)*(*b4 );
				*(++C) = (*a30)*(*b0 ) + (*a31)*(*b1 ) + (*a32)*(*b2 ) + (*a33)*(*b3 ) + (*a34)*(*b4 );
				*(++C) = (*a35)*(*b0 ) + (*a36)*(*b1 ) + (*a37)*(*b2 ) + (*a38)*(*b3 ) + (*a39)*(*b4 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k, b4  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 );
					*(++C) = (*a5 )*(*b0 ) + (*a6 )*(*b1 ) + (*a7 )*(*b2 ) + (*a8 )*(*b3 ) + (*a9 )*(*b4 );
					*(++C) = (*a10)*(*b0 ) + (*a11)*(*b1 ) + (*a12)*(*b2 ) + (*a13)*(*b3 ) + (*a14)*(*b4 );
					*(++C) = (*a15)*(*b0 ) + (*a16)*(*b1 ) + (*a17)*(*b2 ) + (*a18)*(*b3 ) + (*a19)*(*b4 );
					*(++C) = (*a20)*(*b0 ) + (*a21)*(*b1 ) + (*a22)*(*b2 ) + (*a23)*(*b3 ) + (*a24)*(*b4 );
					*(++C) = (*a25)*(*b0 ) + (*a26)*(*b1 ) + (*a27)*(*b2 ) + (*a28)*(*b3 ) + (*a29)*(*b4 );
					*(++C) = (*a30)*(*b0 ) + (*a31)*(*b1 ) + (*a32)*(*b2 ) + (*a33)*(*b3 ) + (*a34)*(*b4 );
					*(++C) = (*a35)*(*b0 ) + (*a36)*(*b1 ) + (*a37)*(*b2 ) + (*a38)*(*b3 ) + (*a39)*(*b4 );
				}
				break; // m8k5
			} case 6: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 , *a5  = A+5 ,
				                *a6  = A+6 , *a7  = A+7 , *a8  = A+8 , *a9  = A+9 , *a10 = A+10, *a11 = A+11,
				                *a12 = A+12, *a13 = A+13, *a14 = A+14, *a15 = A+15, *a16 = A+16, *a17 = A+17,
				                *a18 = A+18, *a19 = A+19, *a20 = A+20, *a21 = A+21, *a22 = A+22, *a23 = A+23,
				                *a24 = A+24, *a25 = A+25, *a26 = A+26, *a27 = A+27, *a28 = A+28, *a29 = A+29,
				                *a30 = A+30, *a31 = A+31, *a32 = A+32, *a33 = A+33, *a34 = A+34, *a35 = A+35,
				                *a36 = A+36, *a37 = A+37, *a38 = A+38, *a39 = A+39, *a40 = A+40, *a41 = A+41,
				                *a42 = A+42, *a43 = A+43, *a44 = A+44, *a45 = A+45, *a46 = A+46, *a47 = A+47,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 , *b5  = B+5 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 );
				*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 ) + (*a8 )*(*b2 ) + (*a9 )*(*b3 ) + (*a10)*(*b4 ) + (*a11)*(*b5 );
				*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 ) + (*a14)*(*b2 ) + (*a15)*(*b3 ) + (*a16)*(*b4 ) + (*a17)*(*b5 );
				*(++C) = (*a18)*(*b0 ) + (*a19)*(*b1 ) + (*a20)*(*b2 ) + (*a21)*(*b3 ) + (*a22)*(*b4 ) + (*a23)*(*b5 );
				*(++C) = (*a24)*(*b0 ) + (*a25)*(*b1 ) + (*a26)*(*b2 ) + (*a27)*(*b3 ) + (*a28)*(*b4 ) + (*a29)*(*b5 );
				*(++C) = (*a30)*(*b0 ) + (*a31)*(*b1 ) + (*a32)*(*b2 ) + (*a33)*(*b3 ) + (*a34)*(*b4 ) + (*a35)*(*b5 );
				*(++C) = (*a36)*(*b0 ) + (*a37)*(*b1 ) + (*a38)*(*b2 ) + (*a39)*(*b3 ) + (*a40)*(*b4 ) + (*a41)*(*b5 );
				*(++C) = (*a42)*(*b0 ) + (*a43)*(*b1 ) + (*a44)*(*b2 ) + (*a45)*(*b3 ) + (*a46)*(*b4 ) + (*a47)*(*b5 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k, b4  += k, b5  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 );
					*(++C) = (*a6 )*(*b0 ) + (*a7 )*(*b1 ) + (*a8 )*(*b2 ) + (*a9 )*(*b3 ) + (*a10)*(*b4 ) + (*a11)*(*b5 );
					*(++C) = (*a12)*(*b0 ) + (*a13)*(*b1 ) + (*a14)*(*b2 ) + (*a15)*(*b3 ) + (*a16)*(*b4 ) + (*a17)*(*b5 );
					*(++C) = (*a18)*(*b0 ) + (*a19)*(*b1 ) + (*a20)*(*b2 ) + (*a21)*(*b3 ) + (*a22)*(*b4 ) + (*a23)*(*b5 );
					*(++C) = (*a24)*(*b0 ) + (*a25)*(*b1 ) + (*a26)*(*b2 ) + (*a27)*(*b3 ) + (*a28)*(*b4 ) + (*a29)*(*b5 );
					*(++C) = (*a30)*(*b0 ) + (*a31)*(*b1 ) + (*a32)*(*b2 ) + (*a33)*(*b3 ) + (*a34)*(*b4 ) + (*a35)*(*b5 );
					*(++C) = (*a36)*(*b0 ) + (*a37)*(*b1 ) + (*a38)*(*b2 ) + (*a39)*(*b3 ) + (*a40)*(*b4 ) + (*a41)*(*b5 );
					*(++C) = (*a42)*(*b0 ) + (*a43)*(*b1 ) + (*a44)*(*b2 ) + (*a45)*(*b3 ) + (*a46)*(*b4 ) + (*a47)*(*b5 );
				}
				break; // m8k6
			} case 7: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 , *a5  = A+5 , *a6  = A+6 ,
				                *a7  = A+7 , *a8  = A+8 , *a9  = A+9 , *a10 = A+10, *a11 = A+11, *a12 = A+12, *a13 = A+13,
				                *a14 = A+14, *a15 = A+15, *a16 = A+16, *a17 = A+17, *a18 = A+18, *a19 = A+19, *a20 = A+20,
				                *a21 = A+21, *a22 = A+22, *a23 = A+23, *a24 = A+24, *a25 = A+25, *a26 = A+26, *a27 = A+27,
				                *a28 = A+28, *a29 = A+29, *a30 = A+30, *a31 = A+31, *a32 = A+32, *a33 = A+33, *a34 = A+34,
				                *a35 = A+35, *a36 = A+36, *a37 = A+37, *a38 = A+38, *a39 = A+39, *a40 = A+40, *a41 = A+41,
				                *a42 = A+42, *a43 = A+43, *a44 = A+44, *a45 = A+45, *a46 = A+46, *a47 = A+47, *a48 = A+48,
				                *a49 = A+49, *a50 = A+50, *a51 = A+51, *a52 = A+52, *a53 = A+53, *a54 = A+54, *a55 = A+55,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 , *b5  = B+5 , *b6  = B+6 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 ) + (*a6 )*(*b6 );
				*(++C) = (*a7 )*(*b0 ) + (*a8 )*(*b1 ) + (*a9 )*(*b2 ) + (*a10)*(*b3 ) + (*a11)*(*b4 ) + (*a12)*(*b5 ) + (*a13)*(*b6 );
				*(++C) = (*a14)*(*b0 ) + (*a15)*(*b1 ) + (*a16)*(*b2 ) + (*a17)*(*b3 ) + (*a18)*(*b4 ) + (*a19)*(*b5 ) + (*a20)*(*b6 );
				*(++C) = (*a21)*(*b0 ) + (*a22)*(*b1 ) + (*a23)*(*b2 ) + (*a24)*(*b3 ) + (*a25)*(*b4 ) + (*a26)*(*b5 ) + (*a27)*(*b6 );
				*(++C) = (*a28)*(*b0 ) + (*a29)*(*b1 ) + (*a30)*(*b2 ) + (*a31)*(*b3 ) + (*a32)*(*b4 ) + (*a33)*(*b5 ) + (*a34)*(*b6 );
				*(++C) = (*a35)*(*b0 ) + (*a36)*(*b1 ) + (*a37)*(*b2 ) + (*a38)*(*b3 ) + (*a39)*(*b4 ) + (*a40)*(*b5 ) + (*a41)*(*b6 );
				*(++C) = (*a42)*(*b0 ) + (*a43)*(*b1 ) + (*a44)*(*b2 ) + (*a45)*(*b3 ) + (*a46)*(*b4 ) + (*a47)*(*b5 ) + (*a48)*(*b6 );
				*(++C) = (*a49)*(*b0 ) + (*a50)*(*b1 ) + (*a51)*(*b2 ) + (*a52)*(*b3 ) + (*a53)*(*b4 ) + (*a54)*(*b5 ) + (*a55)*(*b6 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k, b4  += k, b5  += k, b6  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 ) + (*a6 )*(*b6 );
					*(++C) = (*a7 )*(*b0 ) + (*a8 )*(*b1 ) + (*a9 )*(*b2 ) + (*a10)*(*b3 ) + (*a11)*(*b4 ) + (*a12)*(*b5 ) + (*a13)*(*b6 );
					*(++C) = (*a14)*(*b0 ) + (*a15)*(*b1 ) + (*a16)*(*b2 ) + (*a17)*(*b3 ) + (*a18)*(*b4 ) + (*a19)*(*b5 ) + (*a20)*(*b6 );
					*(++C) = (*a21)*(*b0 ) + (*a22)*(*b1 ) + (*a23)*(*b2 ) + (*a24)*(*b3 ) + (*a25)*(*b4 ) + (*a26)*(*b5 ) + (*a27)*(*b6 );
					*(++C) = (*a28)*(*b0 ) + (*a29)*(*b1 ) + (*a30)*(*b2 ) + (*a31)*(*b3 ) + (*a32)*(*b4 ) + (*a33)*(*b5 ) + (*a34)*(*b6 );
					*(++C) = (*a35)*(*b0 ) + (*a36)*(*b1 ) + (*a37)*(*b2 ) + (*a38)*(*b3 ) + (*a39)*(*b4 ) + (*a40)*(*b5 ) + (*a41)*(*b6 );
					*(++C) = (*a42)*(*b0 ) + (*a43)*(*b1 ) + (*a44)*(*b2 ) + (*a45)*(*b3 ) + (*a46)*(*b4 ) + (*a47)*(*b5 ) + (*a48)*(*b6 );
					*(++C) = (*a49)*(*b0 ) + (*a50)*(*b1 ) + (*a51)*(*b2 ) + (*a52)*(*b3 ) + (*a53)*(*b4 ) + (*a54)*(*b5 ) + (*a55)*(*b6 );
				}
				break; // m8k7
			} case 8: {
				register const double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 , *a5  = A+5 , *a6  = A+6 , *a7  = A+7 ,
				                *a8  = A+8 , *a9  = A+9 , *a10 = A+10, *a11 = A+11, *a12 = A+12, *a13 = A+13, *a14 = A+14, *a15 = A+15,
				                *a16 = A+16, *a17 = A+17, *a18 = A+18, *a19 = A+19, *a20 = A+20, *a21 = A+21, *a22 = A+22, *a23 = A+23,
				                *a24 = A+24, *a25 = A+25, *a26 = A+26, *a27 = A+27, *a28 = A+28, *a29 = A+29, *a30 = A+30, *a31 = A+31,
				                *a32 = A+32, *a33 = A+33, *a34 = A+34, *a35 = A+35, *a36 = A+36, *a37 = A+37, *a38 = A+38, *a39 = A+39,
				                *a40 = A+40, *a41 = A+41, *a42 = A+42, *a43 = A+43, *a44 = A+44, *a45 = A+45, *a46 = A+46, *a47 = A+47,
				                *a48 = A+48, *a49 = A+49, *a50 = A+50, *a51 = A+51, *a52 = A+52, *a53 = A+53, *a54 = A+54, *a55 = A+55,
				                *a56 = A+56, *a57 = A+57, *a58 = A+58, *a59 = A+59, *a60 = A+60, *a61 = A+61, *a62 = A+62, *a63 = A+63,
				                *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 , *b5  = B+5 , *b6  = B+6 , *b7  = B+7 ;

				*C     = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 ) + (*a6 )*(*b6 ) + (*a7 )*(*b7 );
				*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 ) + (*a10)*(*b2 ) + (*a11)*(*b3 ) + (*a12)*(*b4 ) + (*a13)*(*b5 ) + (*a14)*(*b6 ) + (*a15)*(*b7 );
				*(++C) = (*a16)*(*b0 ) + (*a17)*(*b1 ) + (*a18)*(*b2 ) + (*a19)*(*b3 ) + (*a20)*(*b4 ) + (*a21)*(*b5 ) + (*a22)*(*b6 ) + (*a23)*(*b7 );
				*(++C) = (*a24)*(*b0 ) + (*a25)*(*b1 ) + (*a26)*(*b2 ) + (*a27)*(*b3 ) + (*a28)*(*b4 ) + (*a29)*(*b5 ) + (*a30)*(*b6 ) + (*a31)*(*b7 );
				*(++C) = (*a32)*(*b0 ) + (*a33)*(*b1 ) + (*a34)*(*b2 ) + (*a35)*(*b3 ) + (*a36)*(*b4 ) + (*a37)*(*b5 ) + (*a38)*(*b6 ) + (*a39)*(*b7 );
				*(++C) = (*a40)*(*b0 ) + (*a41)*(*b1 ) + (*a42)*(*b2 ) + (*a43)*(*b3 ) + (*a44)*(*b4 ) + (*a45)*(*b5 ) + (*a46)*(*b6 ) + (*a47)*(*b7 );
				*(++C) = (*a48)*(*b0 ) + (*a49)*(*b1 ) + (*a50)*(*b2 ) + (*a51)*(*b3 ) + (*a52)*(*b4 ) + (*a53)*(*b5 ) + (*a54)*(*b6 ) + (*a55)*(*b7 );
				*(++C) = (*a56)*(*b0 ) + (*a57)*(*b1 ) + (*a58)*(*b2 ) + (*a59)*(*b3 ) + (*a60)*(*b4 ) + (*a61)*(*b5 ) + (*a62)*(*b6 ) + (*a63)*(*b7 );

				for (register unsigned int nRem = n-1; nRem--; ) {
					b0  += k, b1  += k, b2  += k, b3  += k, b4  += k, b5  += k, b6  += k, b7  += k;

					*(++C) = (*a0 )*(*b0 ) + (*a1 )*(*b1 ) + (*a2 )*(*b2 ) + (*a3 )*(*b3 ) + (*a4 )*(*b4 ) + (*a5 )*(*b5 ) + (*a6 )*(*b6 ) + (*a7 )*(*b7 );
					*(++C) = (*a8 )*(*b0 ) + (*a9 )*(*b1 ) + (*a10)*(*b2 ) + (*a11)*(*b3 ) + (*a12)*(*b4 ) + (*a13)*(*b5 ) + (*a14)*(*b6 ) + (*a15)*(*b7 );
					*(++C) = (*a16)*(*b0 ) + (*a17)*(*b1 ) + (*a18)*(*b2 ) + (*a19)*(*b3 ) + (*a20)*(*b4 ) + (*a21)*(*b5 ) + (*a22)*(*b6 ) + (*a23)*(*b7 );
					*(++C) = (*a24)*(*b0 ) + (*a25)*(*b1 ) + (*a26)*(*b2 ) + (*a27)*(*b3 ) + (*a28)*(*b4 ) + (*a29)*(*b5 ) + (*a30)*(*b6 ) + (*a31)*(*b7 );
					*(++C) = (*a32)*(*b0 ) + (*a33)*(*b1 ) + (*a34)*(*b2 ) + (*a35)*(*b3 ) + (*a36)*(*b4 ) + (*a37)*(*b5 ) + (*a38)*(*b6 ) + (*a39)*(*b7 );
					*(++C) = (*a40)*(*b0 ) + (*a41)*(*b1 ) + (*a42)*(*b2 ) + (*a43)*(*b3 ) + (*a44)*(*b4 ) + (*a45)*(*b5 ) + (*a46)*(*b6 ) + (*a47)*(*b7 );
					*(++C) = (*a48)*(*b0 ) + (*a49)*(*b1 ) + (*a50)*(*b2 ) + (*a51)*(*b3 ) + (*a52)*(*b4 ) + (*a53)*(*b5 ) + (*a54)*(*b6 ) + (*a55)*(*b7 );
					*(++C) = (*a56)*(*b0 ) + (*a57)*(*b1 ) + (*a58)*(*b2 ) + (*a59)*(*b3 ) + (*a60)*(*b4 ) + (*a61)*(*b5 ) + (*a62)*(*b6 ) + (*a63)*(*b7 );
				}
				break; // m8k8
			} default: {
				printf("Error: Unsupported m = %d, k = %d, combination in mm_CTN_d (useBlas = 0).\n",m,k), EXIT_MSG;
				break; // m8
			}}
			break;
		} default: {
			printf("Error: Unsupported m = %d, in mm_CTN_d (useBlas = 0).\n",m), EXIT_MSG;
			break;
		}}
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
		printf("Error: Unsupported value of useBLAS (%d) in mm_CTN_d.\n",useBLAS), EXIT_MSG;
		break;
	}
}

void mm_CTN_CSR_d(const int m, const int n, const int k, const struct S_OpCSR *A, const double *B, double *C)
{
	/*
	 *	Purpose:
	 *		Compute matrix-matrix and matrix-vector products:
	 *			0) return C = A'*B
	 *			1) Row-major storage of A
	 *			2) Column-major storage of B and C
	 *
	 *	Comments:
	 *		Comparison with mkl_dcsrmm showed that the c code was faster for all cases considered, becoming increasingly
	 *		better as the size of A was increased. This was surprising and will be further investigated in the future.
	 *		As a result of the tests, only the c code is used here. (ToBeModified)
	 *
	 *		Note that entries 2 and 3 of matdescra are not used in this function. Further, entries 5 and 6 were reserved
	 *		for future use by the intel MKL developers at the time that this function was written.
	 *		Transposing of input B and then outputs B and C is required when using mkl_dcsrmm which seems to expect
	 *		row-major ordering of all entries.
	 */

	// Standard datatypes
	register unsigned int *rowIndex, *columns;
	register const double       *values;

	register double *C_ptr = C;
	for (register unsigned int l = 0, lMax = n; l < lMax; l++) {
		register const double *B_ptr = &B[l*k];
		rowIndex = A->rowIndex;
		columns  = A->columns;
		values   = A->values;

		for (register unsigned int i = 0, iMax = m; i < iMax; i++) {
			*C_ptr = 0.0;
			for (register unsigned int jInd = *rowIndex, jIndMax = *(++rowIndex); jInd < jIndMax; jInd++)
				*C_ptr += (*values++)*B_ptr[*columns++];
			C_ptr++;
		}
	}

/*
	// Code for mkl_dcsrmm function call
	char transa = 'N', matdescra[6] = {'G','L','N','C'};
	double alpha = 1.0, beta = 0.0;

	rowIndex = A->rowIndex;
	columns  = A->columns;
	values   = A->values;

	mkl_dimatcopy('C','T',k,n,1.0,B,k,n);
	mkl_dcsrmm(&transa,&m,&n,&k,&alpha,matdescra,values,(int*)columns,(int*)rowIndex,(int*)&(rowIndex[1]),B,&n,&beta,C,&n);
	mkl_dimatcopy('C','T',n,k,1.0,B,n,k);
	mkl_dimatcopy('C','T',n,m,1.0,C,n,m);
*/
}

void convert_to_CSR_d(const unsigned int NRows, const unsigned int NCols, const double *Input, struct S_OpCSR **Output)
{
	/*
	 *	Purpose:
	 *		Return rowIndex, columns, and values of sparse operator in compressed sparse row (CSR) format.
	 *
	 *	Comments:
	 *
	 *	References:
	 *		Intel MKL Sparse BLAS CSR Matrix Storage Format: https://software.intel.com/en-us/node/599835
	 */

	unsigned int i, iInd, j, Nval,
	             *rowIndex, *columns, *columnsOver;
	double       *values, *valuesOver, tmp_d;

	*Output = malloc(sizeof **Output); // keep

	Nval = 0;
	rowIndex    = malloc((NRows+1)   * sizeof *rowIndex);    // keep
	columnsOver = malloc(NRows*NCols * sizeof *columnsOver); // free
	valuesOver  = malloc(NRows*NCols * sizeof *valuesOver);  // free

	rowIndex[0] = 0;
	for (i = 0; i < NRows; i++) {
		iInd = i*NCols;
		for (j = 0; j < NCols; j++) {
			tmp_d = Input[iInd+j];
			if (fabs(tmp_d) > EPS) {
				columnsOver[Nval] = j;
				valuesOver[Nval] = tmp_d;
				Nval++;
			}
		}
		rowIndex[i+1] = Nval;
	}

	columns = malloc(Nval * sizeof *columns); // keep
	values  = malloc(Nval * sizeof *values);  // keep

	for (i = 0; i < Nval; i++) {
		columns[i] = columnsOver[i];
		values[i]  = valuesOver[i];
	}

	free(columnsOver);
	free(valuesOver);

	(*Output)->NRows = NRows;
	(*Output)->NVals = Nval;
	(*Output)->rowIndex = rowIndex;
	(*Output)->columns = columns;
	(*Output)->values = values;
}

struct S_MATRIX *mm_alloc_mat (char const layout, char const opA, char const opB, struct S_MATRIX const *const A,
                               struct S_MATRIX const *const B)
{
	/*
	 *	Purpose:
	 *		Compute C = opA(A)*opB(B) in the appropriate layout using the matrix struct.
	 *
	 *	Comments:
	 *		Note that a matrix in the alternate layout is simply interpreted as being transposed in the current layout.
	 *
	 *	Notation:
	 *		Op(A/B) : 'N'(o transpose), 'T'(ranspose)
	 */

	if (A->format != 'D' || B->format != 'D')
		EXIT_UNSUPPORTED;

	CBLAS_LAYOUT    const CBlayout = ( layout == 'R' ? CBRM : CBCM );
	CBLAS_TRANSPOSE const transa   = ( (layout == A->layout) == (opA == 'N') ? CBNT : CBT ),
	                      transb   = ( (layout == B->layout) == (opB == 'N') ? CBNT : CBT );
	size_t          const m        = ( transa == CBNT ? A->NRows : A->NCols ),
	                      n        = ( transb == CBNT ? B->NCols : B->NRows ),
	                      k        = ( transa == CBNT ? A->NCols : A->NRows );
	double          const alpha    = 1.0;

	// Check matching internal length
	if (k != ( transb == CBNT ? B->NRows : B->NCols ))
		EXIT_UNSUPPORTED;

	struct S_MATRIX *C = malloc(sizeof *C); // keep

	C->layout = layout;
	C->format = 'D';
	C->NRows  = m;
	C->NCols  = n;
	C->values = mm_Alloc_d(CBlayout,transa,transb,m,n,k,alpha,A->values,B->values); // keep

	return C;
}

static struct S_MATRIX *mat_constructor_dense (const char layout, const size_t NRows, const size_t NCols)
{
	struct S_MATRIX *A = malloc(sizeof *A); // keep

	A->layout = layout;
	A->format = 'D';
	A->NRows  = NRows;
	A->NCols  = NCols;
	A->values = calloc(NRows*NCols , sizeof *(A->values)); // keep

	return A;
}

struct S_MATRIX *mat_copy (struct S_MATRIX const *const A)
{
	struct S_MATRIX *B = malloc(sizeof *B); // keep

	B->layout = A->layout;
	B->format = A->format;
	B->NRows  = A->NRows;
	B->NCols  = A->NCols;

	size_t const Nvals = (B->NRows)*(B->NCols);

	B->values = malloc(Nvals * sizeof *(B->values)); // keep
	for (size_t i = 0; i < Nvals; i++)
		B->values[i] = A->values[i];

	return B;
}

void mat_transpose (struct S_MATRIX *const A)
{
	if (A->format != 'D')
		EXIT_UNSUPPORTED;

	mkl_dimatcopy(A->layout,'T',A->NRows,A->NCols,1.0,A->values,A->NCols,A->NRows);

	size_t const swap = A->NRows;
	A->NRows = A->NCols;
	A->NCols = swap;

	A->layout = (A->layout == 'R' ? 'C' : 'R');
}

struct S_MATRIX *mm_diag_mat_alloc (char const layout, char const opA, char const side, struct S_MATRIX const *const A,
                                    struct S_VECTOR const *const a)
{
	/*
	 *	Comments:
	 *		B is computed with the layout of A and transposed if necessary.
	 */

	if (A->format != 'D')
		EXIT_UNSUPPORTED;

	struct S_MATRIX *B = mat_constructor_dense(A->layout,A->NRows,A->NCols); // keep

	char const side_to_use = ( opA == 'N' ? side : ( side == 'L' ? 'R' : 'L' ) );

	if (a->NRows != ( side_to_use == 'L' ? A->NRows : A->NCols))
		EXIT_UNSUPPORTED;

	mm_diag_d(B->NRows,B->NCols,a->values,A->values,B->values,1.0,0.0,side_to_use,A->layout);

	// Transpose B if necessary
	if ((layout == A->layout) != (opA == 'N'))
		mat_transpose(B);

	return B;
}

struct S_MATRIX *inverse_mat (struct S_MATRIX const *const A)
{
	if (A->NRows != A->NCols)
		EXIT_UNSUPPORTED;

	struct S_MATRIX *AInv = malloc(sizeof *AInv); // keep

	AInv->layout  = A->layout;
	AInv->format  = A->format;
	AInv->NRows   = A->NRows;
	AInv->NCols   = A->NCols;

	double *I_d = identity_d(A->NRows); // free
	AInv->values = inverse_d(A->NRows,A->NCols,A->values,I_d); // keep
	free(I_d);

	return AInv;
}

struct S_MATRIX *identity_mat (unsigned int const N)
{
	struct S_MATRIX *A = malloc(sizeof *A); // keep

	A->layout = 'R';
	A->format = 'D';
	A->NRows  = N;
	A->NCols  = N;
	A->values = identity_d(N); // keep

	return A;
}

struct S_MATRIX *mat_constructor_move (char const layout, char const format, size_t const NRows, size_t const NCols,
                                       double *const values, struct S_OpCSR *const A_CSR)
{
	struct S_MATRIX *A = malloc(sizeof *A); // keep

	A->layout = layout;
	A->format = format;
	A->NRows  = NRows;
	A->NCols  = NCols;

	if (format == 'D') {
		A->values = values;
	} else if (format == 'S') {
		A->values   = A_CSR->values;
		A->rowIndex = A_CSR->rowIndex;
		A->columns  = A_CSR->columns;
	} else {
		EXIT_UNSUPPORTED;
	}

	return A;
}

struct S_MATRIX **mat_constructor2_move (char const layout, char const format, size_t const NRows, size_t const NCols,
                                         size_t const dim1, double *const *const values, struct S_OpCSR *const *const A_CSR)
{
	struct S_MATRIX **A = malloc(dim1 * sizeof *A); // keep

	for (size_t i = 0; i < dim1; i++) {
		A[i] = malloc(sizeof *A[i]); // keep
		if (format == 'D')
			mat_constructor_move(layout,format,NRows,NCols,values[i],NULL);
		else if (format == 'S')
			mat_constructor_move(layout,format,NRows,NCols,NULL,A_CSR[i]);
		else
			EXIT_UNSUPPORTED;
	}

	return A;
}


struct S_VECTOR *vec_constructor_move (size_t const NRows, double *const values)
{
	struct S_VECTOR *a = malloc(sizeof *a); // keep

	a->NRows  = NRows;
	a->values = values;

	return a;
}
