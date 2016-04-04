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

void mm_d(const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const int m, const int n, const int k,
	      const double alpha, const double *A, const double *B, const double *C)
{
	/*
	 *	Purpose:
	 *		Returns: C = alpha*op(A)*op(B) with memory already allocated in calling function.
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

	cblas_dgemm(CblasRowMajor,transa,transb,m_MKL,n_MKL,k_MKL,alpha,A,ldA,B,ldB,0.0,C,ldC);
}



#ifdef DEBUG

int main(int argc, char **argv)
{
	printf("Entered\n");

	double *A, *B, *C;
	int m, n, k;

	m = 2;
	n = 4;
	k = 3;

	A = malloc(m*k * sizeof *A);
	B = malloc(k*n * sizeof *B);
	C = malloc(m*n * sizeof *C);

	A[0] = 8; A[1] = 1; A[2] = 6;
	A[3] = 3; A[4] = 5; A[5] = 7;

	B[0] = 8.1; B[1] = 5.1; B[2]  = 2.8; B[3]  = 0.4;
	B[4] = 9.1; B[5] = 6.3; B[6]  = 5.4; B[7]  = 1.7;
	B[8] = 1.3; B[9] = 9.8; B[10] = 9.6; B[11] = 2.9;
	mm_d(CblasNoTrans,CblasNoTrans,m,n,k,1.0,A,B,C);

	B[0] = 8.1; B[3] = 5.1; B[6] = 2.8; B[9]  = 0.4;
	B[1] = 9.1; B[4] = 6.3; B[7] = 5.4; B[10] = 1.7;
	B[2] = 1.3; B[5] = 9.8; B[8] = 9.6; B[11] = 2.9;
//	mm_d(CblasNoTrans,CblasTrans,m,n,k,1.0,A,B,C);

	int i, j;

	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++)
			printf("% .3e ",C[i*n+j]);
		printf("\n");
	}
	printf("\n");

	return 0;
}

#endif
