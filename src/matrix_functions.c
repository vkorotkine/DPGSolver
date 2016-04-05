#include <stdlib.h>
#include <stdio.h>

#include "mkl.h"

/*
 *	Purpose:
 *		Provide matrix functions:
 *			double *identity_d(const int N);
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

double *identity_d(const int N)
{
	int i, j;
	double *I;

	I = malloc(N*N * sizeof(double)); // keep (requires external free)
	for (i = 0; i < N; i++) {
	for (j = 0; j < N; j++) {
		if (i == j) I[i*N+j] = 1.0;
		else        I[i*N+j] = 0.0;
	}}
	return I;
}

double *inverse_d(int N, int NRHS, double *A, double *b)
{
	int i, iMax;
	double *x;

	lapack_int N_LA, NRHS_LA, ipiv[N], info;

	N_LA    = (lapack_int) N;
	NRHS_LA = (lapack_int) NRHS;

	x = malloc(N*NRHS * sizeof(double)); // keep (requires external free)
	for (i = 0, iMax = N*NRHS; i < iMax; i++)
		x[i] = b[i];

	info = LAPACKE_dgesv(LAPACK_ROW_MAJOR,N_LA,NRHS_LA,A,N_LA,ipiv,x,NRHS_LA);
	if (info > 0) {
		printf("The diagonal element of the triangular factor of A,\n");
		printf("U(%i,%i) is zero, so that A is singular;\n", info, info);
		printf("the solution could not be computed.\n");
		exit(1);
	}
	return x;
}

double *mm_Alloc_d(const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const int m, const int n, const int k,
	               const double alpha, const double *A, const double *B)
{
	/*
	 *	Purpose:
	 *		Returns: C = alpha*op(A)*op(B)
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

	int C_NRows, C_NCols;
	double *C;
	MKL_INT m_MKL, n_MKL, k_MKL,
	        ldA, ldB, ldC;

	m_MKL = (MKL_INT) m;
	n_MKL = (MKL_INT) n;
	k_MKL = (MKL_INT) k;

	if (transa == CblasNoTrans) {
		ldA = k_MKL;
		C_NRows = m;
	} else {
		ldA = m_MKL;
		C_NRows = k;
	}

	if (transb == CblasNoTrans) {
		ldB = n_MKL;
		C_NCols = k;
	} else {
		ldB = k_MKL;
		C_NCols = n;
	}

	ldC = n_MKL;

	C = malloc(C_NRows*C_NCols * sizeof(double)); // keep (requires external free)

	cblas_dgemm(CblasRowMajor,transa,transb,m_MKL,n_MKL,k_MKL,alpha,A,ldA,B,ldB,0.0,C,ldC);

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

	MKL_INT m_MKL, n_MKL, k_MKL,
	        ldA, ldB, ldC;

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
