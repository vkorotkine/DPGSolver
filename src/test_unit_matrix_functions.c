// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "test_unit_matrix_functions.h"

#include <stdlib.h>
#include <stdio.h>

#include "mkl.h"

#include "Parameters.h"
#include "Test.h"
#include "S_OpCSR.h"

#include "test_support.h"
#include "array_norm.h"
#include "matrix_functions.h"
#include "array_free.h"
#include "array_print.h"

#undef I // No complex variables used here

/*
 *	Purpose:
 *		Test correctness of implementation of matrix functions.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

static void test_unit_matrix_diag(void)
{
	unsigned int pass;

	/*
	 *	diag_d:
	 *
	 *		Input:
	 *
	 *			N = 4
	 *			x = [1.0 2.0 3.0 4.0]
	 *
	 *		Expected Output:
	 *
	 *			X = [1.0 0.0 0.0 0.0
	 *			     0.0 2.0 0.0 0.0
	 *			     0.0 0.0 3.0 0.0
	 *			     0.0 0.0 0.0 4.0]
	 */

	unsigned int N = 4;
	double *X,
	       x[4] = { 1.0, 2.0, 3.0, 4.0 },
	       X4[16] = {1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 4.0};

	X = diag_d(x,N); // free

	pass = 0;
	if (array_norm_diff_d(N*N,X,X4,"Inf") < EPS)
		pass = 1;

	test_print2(pass,"matrix_diag_d:");

	free(X);
}

static void test_unit_matrix_identity(void)
{
	unsigned int pass;

	/*
	 *	identity_d:
	 *
	 *		Input:
	 *
	 *			N = 4
	 *
	 *		Expected Output:
	 *
	 *			I = [1.0 0.0 0.0 0.0
	 *			     0.0 1.0 0.0 0.0
	 *			     0.0 0.0 1.0 0.0
	 *			     0.0 0.0 0.0 1.0]
	 */

	unsigned int N = 4;
	double *I,
	       I4[16] = {1.0, 0.0, 0.0 , 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};

	I = identity_d(N); // free

	pass = 0;
	if (array_norm_diff_d(16,I,I4,"Inf") < EPS)
		pass = 1;

	test_print2(pass,"matrix_identity_d:");

	free(I);
}

static void test_unit_matrix_inverse(void)
{
	unsigned int pass;

	/*
	 *	inverse_d:
	 *
	 *		Input:
	 *
	 *			A = [0.81 0.91 0.28
	 *			     0.91 0.63 0.55
	 *			     0.13 0.10 0.96]
	 *
	 *		Expected Output:
	 *
	 *			AInv = [-1.949472564488963  2.998315752149631 -1.149188901693112
	 *			         2.844074106905415 -2.628135803563513  0.676181189610850
	 *			        -0.032266643028100 -0.132257778565730  1.126850456519812]
	 */

	unsigned int N = 3;
	double *I3, *Aroundtrip_c, *AInv_c,
	       A[9]    = { 0.81, 0.91, 0.28, 0.91, 0.63, 0.55, 0.13, 0.10, 0.96 },
	       AInv[9] = { -1.949472564488963,  2.998315752149631, -1.149188901693112,
	                    2.844074106905415, -2.628135803563513,  0.676181189610850,
	                   -0.032266643028100, -0.132257778565730,  1.126850456519812 };

	I3 = identity_d(N); // free

	AInv_c       = inverse_d(3,3,A,I3);      // free
	Aroundtrip_c = inverse_d(N,N,AInv_c,I3); // free

	pass = 0;
	if (array_norm_diff_d(9,AInv,AInv_c,"Inf")    < EPS &&
	    array_norm_diff_d(9,A,Aroundtrip_c,"Inf") < EPS)
		pass = 1;

	test_print2(pass,"matrix_inverse_d:");

	free(I3);
	free(AInv_c);
	free(Aroundtrip_c);
}

static void test_unit_matrix_mm(void)
{
	unsigned int pass;

	/*
	 *	mm_d (permutations):
	 *
	 *		Input:
	 *
	 *			A = [0.05 0.69 0.03 0.77
	 *			     0.10 0.32 0.44 0.80
	 *			     0.82 0.95 0.38 0.19]
	 *
	 *			B = [0.49 0.75
	 *			     0.45 0.28
	 *			     0.65 0.68
	 *			     0.71 0.66]
	 *
	 *			C = [0.9012 0.7593
	 *			     1.0470 0.9918
	 *			     1.2112 1.2648]
	 *
	 *			AT = [0.05 0.10 0.82
	 *			      0.69 0.32 0.95
	 *			      0.03 0.44 0.38
	 *			      0.77 0.80 0.19]
	 *
	 *			BT = [0.49 0.45 0.65 0.71
	 *			      0.75 0.28 0.68 0.66]
	 *
	 *			CT = [0.9012 1.0470 1.2112
	 *			      0.7593 0.9918 1.2648]
	 *
	 *
	 *		Expected Output:
	 *
	 *			Input: CblasRowMajor (i.e. all results interpreted as having a row-major layout)
	 *
	 *				C = A  *B
	 *				C = AT'*B
	 *				C = A  *BT'
	 *				C = AT'*BT'
	 *
	 *			Input: CblasRowMajor (i.e. all results interpreted as having a col-major layout)
	 *
	 *				CT = AT *BT
	 *				CT = A' *BT
	 *				CT = AT *B'
	 *				CT = A'*B'
	 */

	double A[12]  = { 0.05, 0.69, 0.03, 0.77, 0.10, 0.32, 0.44, 0.8, 0.82, 0.95, 0.38, 0.19 },
	       B[8]   = { 0.49, 0.75, 0.45, 0.28, 0.65, 0.68, 0.71, 0.66 },
	       C[6]   = { 0.9012, 0.7593, 1.0470, 0.9918, 1.2112, 1.2648 },
	       AT[12] = { 0.05, 0.10, 0.82, 0.69, 0.32, 0.95, 0.03, 0.44, 0.38, 0.77, 0.80, 0.19 },
	       BT[8]  = { 0.49, 0.45, 0.65, 0.71, 0.75, 0.28, 0.68, 0.66 },
	       CT[6]  = { 0.9012, 1.0470, 1.2112, 0.7593, 0.9918, 1.2648 };
	double *C_c, *CT_c;

	C_c  = malloc(6 * sizeof *C_c);  // free
	CT_c = malloc(6 * sizeof *CT_c); // free

	// Row-Major

	// NoTrans, NoTrans
	mm_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,3,2,4,1.0,0.0,A,B,C_c);

	pass = 0;
	if (array_norm_diff_d(6,C,C_c,"Inf") < EPS)
		pass = 1;

	test_print2(pass,"matrix_mm_d (RowMajor, NoTrans, NoTrans):");


	// Trans, NoTrans
	mm_d(CblasRowMajor,CblasTrans,CblasNoTrans,3,2,4,1.0,0.0,AT,B,C_c);

	pass = 0;
	if (array_norm_diff_d(6,C,C_c,"Inf") < EPS)
		pass = 1;

	test_print2(pass,"            (RowMajor,   Trans, NoTrans):");


	// NoTrans, Trans
	mm_d(CblasRowMajor,CblasNoTrans,CblasTrans,3,2,4,1.0,0.0,A,BT,C_c);

	pass = 0;
	if (array_norm_diff_d(6,C,C_c,"Inf") < EPS)
		pass = 1;

	test_print2(pass,"            (RowMajor, NoTrans,   Trans):");


	// Trans, Trans
	mm_d(CblasRowMajor,CblasTrans,CblasTrans,3,2,4,1.0,0.0,AT,BT,C_c);

	pass = 0;
	if (array_norm_diff_d(6,C,C_c,"Inf") < EPS)
		pass = 1;

	test_print2(pass,"            (RowMajor,   Trans,   Trans):");



	// Col-Major

	// NoTrans, NoTrans
	mm_d(CblasColMajor,CblasNoTrans,CblasNoTrans,3,2,4,1.0,0.0,AT,BT,CT_c);

	pass = 0;
	if (array_norm_diff_d(6,CT,CT_c,"Inf") < EPS)
		pass = 1;

	test_print2(pass,"            (ColMajor, NoTrans, NoTrans):");


	// Trans, NoTrans
	mm_d(CblasColMajor,CblasTrans,CblasNoTrans,3,2,4,1.0,0.0,A,BT,CT_c);

	pass = 0;
	if (array_norm_diff_d(6,CT,CT_c,"Inf") < EPS)
		pass = 1;

	test_print2(pass,"            (ColMajor,   Trans, NoTrans):");


	// NoTrans, Trans
	mm_d(CblasColMajor,CblasNoTrans,CblasTrans,3,2,4,1.0,0.0,AT,B,CT_c);

	pass = 0;
	if (array_norm_diff_d(6,CT,CT_c,"Inf") < EPS)
		pass = 1;

	test_print2(pass,"            (ColMajor, NoTrans,   Trans):");


	// Trans, Trans
	mm_d(CblasColMajor,CblasTrans,CblasTrans,3,2,4,1.0,0.0,A,B,CT_c);

	pass = 0;
	if (array_norm_diff_d(6,CT,CT_c,"Inf") < EPS)
		pass = 1;

	test_print2(pass,"            (ColMajor,   Trans,   Trans):");

	free(C_c);
	free(CT_c);

	/*
	 *	mm_Alloc_d:
	 *
	 *		Same example as above with alloc in function call.
	 */

	double *C_cNN,  *C_cTN,  *C_cNT,  *C_cTT,
	       *CT_cNN, *CT_cTN, *CT_cNT, *CT_cTT;

	// Row-Major
	C_cNN = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,3,2,4,1.0,A,B); // free
	C_cTN = mm_Alloc_d(CblasRowMajor,CblasTrans,CblasNoTrans,3,2,4,1.0,AT,B);  // free
	C_cNT = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasTrans,3,2,4,1.0,A,BT);  // free
	C_cTT = mm_Alloc_d(CblasRowMajor,CblasTrans,CblasTrans,3,2,4,1.0,AT,BT);   // free

	pass = 0;
	if (array_norm_diff_d(6,C,C_cNN,"Inf") < EPS &&
	    array_norm_diff_d(6,C,C_cTN,"Inf") < EPS &&
	    array_norm_diff_d(6,C,C_cNT,"Inf") < EPS &&
	    array_norm_diff_d(6,C,C_cTT,"Inf") < EPS)
			pass = 1;

	test_print2(pass,"matrix_mm_Alloc_d (RowMajor):");

	// Col-Major
	CT_cNN = mm_Alloc_d(CblasColMajor,CblasNoTrans,CblasNoTrans,3,2,4,1.0,AT,BT); // free
	CT_cTN = mm_Alloc_d(CblasColMajor,CblasTrans,CblasNoTrans,3,2,4,1.0,A,BT);    // free
	CT_cNT = mm_Alloc_d(CblasColMajor,CblasNoTrans,CblasTrans,3,2,4,1.0,AT,B);    // free
	CT_cTT = mm_Alloc_d(CblasColMajor,CblasTrans,CblasTrans,3,2,4,1.0,A,B);       // free

	pass = 0;
	if (array_norm_diff_d(6,CT,CT_cNN,"Inf") < EPS &&
	    array_norm_diff_d(6,CT,CT_cTN,"Inf") < EPS &&
	    array_norm_diff_d(6,CT,CT_cNT,"Inf") < EPS &&
	    array_norm_diff_d(6,CT,CT_cTT,"Inf") < EPS)
			pass = 1;

	test_print2(pass,"                  (ColMajor):");


	free(C_cNN);
	free(C_cTN);
	free(C_cNT);
	free(C_cTT);

	free(CT_cNN);
	free(CT_cTN);
	free(CT_cNT);
	free(CT_cTT);

	/*
	 *	mm_CTN_d (useBlas = 0, 1):
	 *
	 *		Input:
	 *
	 *			mv (useBLAS = 0): m = 3, n = 1, k = 4
	 *			mv (useBLAS = 1): m = 8, n = 1, k = 9
	 *			mm (useBLAS = 0): m = 3, n = 5, k = 4
	 *			mm (useBLAS = 1): m = 8, n = 5, k = 9
	 *
	 *			A = rand(k,m)
	 *			B = rand(k,n)
	 *			C = mm_Alloc_d(CblasColMajor,CblasTrans,CblasNoTrans,m,n,k,1.0,A,B)
	 *
	 *		Expected Output:
	 *
	 *			C_c = A'*B
	 */

	unsigned int i, iMax, m, n, k;
	double *A_CTN, *B_CTN, *C_CTN, *C_CTN_c, RAND_MAXd;

	RAND_MAXd = (double)RAND_MAX;

	printf("\nWarning: Ensure that all options are tested for mm_CTN (mv/mm, useBLAS = 0/1) "
	       "based on input dimensions (m, n, k)\n\n");
	TestDB.Nwarnings++;


	// mv, NoBLAS
	m = 3, n = 1, k = 4;

	A_CTN = malloc(k*m * sizeof *A_CTN); // free
	B_CTN = malloc(k*n * sizeof *B_CTN); // free

	for (i = 0, iMax = k*m; iMax--; i++)
		A_CTN[i] = (double)rand()/RAND_MAXd;

	for (i = 0, iMax = k*n; iMax--; i++)
		B_CTN[i] = (double)rand()/RAND_MAXd;

	C_CTN = mm_Alloc_d(CblasColMajor,CblasTrans,CblasNoTrans,m,n,k,1.0,A_CTN,B_CTN); // free

	C_CTN_c = malloc(m*n * sizeof *C_CTN_c); // free
	mm_CTN_d(m,n,k,A_CTN,B_CTN,C_CTN_c);

	pass = 0;
	if (array_norm_diff_d(m*n,C_CTN,C_CTN_c,"Inf") < EPS)
		pass = 1;

	test_print2(pass,"matrix_mm_CTN_d (mv, NoBLAS):");

	free(A_CTN);
	free(B_CTN);
	free(C_CTN);
	free(C_CTN_c);


	// mv, BLAS
	m = 8, n = 1, k = 9;

	A_CTN = malloc(k*m * sizeof *A_CTN); // free
	B_CTN = malloc(k*n * sizeof *B_CTN); // free

	for (i = 0, iMax = k*m; iMax--; i++)
		A_CTN[i] = (double)rand()/RAND_MAXd;

	for (i = 0, iMax = k*n; iMax--; i++)
		B_CTN[i] = (double)rand()/RAND_MAXd;

	C_CTN = mm_Alloc_d(CblasColMajor,CblasTrans,CblasNoTrans,m,n,k,1.0,A_CTN,B_CTN); // free

	C_CTN_c = malloc(m*n * sizeof *C_CTN_c); // free
	mm_CTN_d(m,n,k,A_CTN,B_CTN,C_CTN_c);

	pass = 0;
	if (array_norm_diff_d(m*n,C_CTN,C_CTN_c,"Inf") < EPS)
		pass = 1;

	test_print2(pass,"                (mv,   BLAS):");

	free(A_CTN);
	free(B_CTN);
	free(C_CTN);
	free(C_CTN_c);


	// mm, NoBLAS
	m = 3, n = 5, k = 4;

	A_CTN = malloc(k*m * sizeof *A_CTN); // free
	B_CTN = malloc(k*n * sizeof *B_CTN); // free

	for (i = 0, iMax = k*m; iMax--; i++)
		A_CTN[i] = (double)rand()/RAND_MAXd;

	for (i = 0, iMax = k*n; iMax--; i++)
		B_CTN[i] = (double)rand()/RAND_MAXd;

	C_CTN = mm_Alloc_d(CblasColMajor,CblasTrans,CblasNoTrans,m,n,k,1.0,A_CTN,B_CTN); // free

	C_CTN_c = malloc(m*n * sizeof *C_CTN_c); // free
	mm_CTN_d(m,n,k,A_CTN,B_CTN,C_CTN_c);

	pass = 0;
	if (array_norm_diff_d(m*n,C_CTN,C_CTN_c,"Inf") < EPS)
		pass = 1;

	test_print2(pass,"                (mm, NoBLAS):");

	free(A_CTN);
	free(B_CTN);
	free(C_CTN);
	free(C_CTN_c);


	// mm, BLAS
	m = 8, n = 5, k = 9;

	A_CTN = malloc(k*m * sizeof *A_CTN); // free
	B_CTN = malloc(k*n * sizeof *B_CTN); // free

	for (i = 0, iMax = k*m; iMax--; i++)
		A_CTN[i] = (double)rand()/RAND_MAXd;

	for (i = 0, iMax = k*n; iMax--; i++)
		B_CTN[i] = (double)rand()/RAND_MAXd;

	C_CTN = mm_Alloc_d(CblasColMajor,CblasTrans,CblasNoTrans,m,n,k,1.0,A_CTN,B_CTN); // free

	C_CTN_c = malloc(m*n * sizeof *C_CTN_c); // free
	mm_CTN_d(m,n,k,A_CTN,B_CTN,C_CTN_c);

	pass = 0;
	if (array_norm_diff_d(m*n,C_CTN,C_CTN_c,"Inf") < EPS)
		pass = 1;

	test_print2(pass,"                (mm,   BLAS):");

	free(A_CTN);
	free(B_CTN);
	free(C_CTN);
	free(C_CTN_c);
}

static void test_unit_matrix_mm_diag(void)
{
	unsigned int pass;

	/*
	 *	diag_d (permutations):
	 *
	 *		Input:
	 *
	 *			A = [0.05 0.69 0.03 0.77
	 *			     0.10 0.32 0.44 0.80
	 *			     0.82 0.95 0.38 0.19]
	 *
	 *		Expected Output:
	 *
	 *			Same output when compared to mm_d.
	 */

	unsigned int NRows = 3, NCols = 4;
	double *DR, *DC, *Odiag, *Omm,
	       A[12]  = { 0.05, 0.69, 0.03, 0.77, 0.10, 0.32, 0.44, 0.80, 0.82, 0.95, 0.38, 0.19 },
	       dR[3]  = { 0.49, 0.45, 0.65 },
	       dC[4]  = { 0.49, 0.45, 0.65, 0.71 };

	DR = diag_d(dR,NRows); // free
	DC = diag_d(dC,NCols); // free

	Odiag = calloc(NRows*NCols , sizeof *Odiag); // free
	Omm   = calloc(NRows*NCols , sizeof *Omm);   // free

	// (L)eft, (R)ow Major
	mm_diag_d(NRows,NCols,dR,A,Odiag,2.0,0.0,'L','R');
	mm_diag_d(NRows,NCols,dR,A,Odiag,2.0,1.0,'L','R');

	mm_d(CBRM,CBNT,CBNT,NRows,NCols,NRows,2.0,0.0,DR,A,Omm);
	mm_d(CBRM,CBNT,CBNT,NRows,NCols,NRows,2.0,1.0,DR,A,Omm);

	pass = 0;
	if (array_norm_diff_d(NRows*NCols,Odiag,Omm,"Inf") < EPS)
		pass = 1;

	test_print2(pass,"matrix_mm_diag_d ('L','R'):");

	// (R)ight, (R)ow Major
	mm_diag_d(NRows,NCols,dC,A,Odiag,0.5,2.0,'R','R');
	mm_d(CBRM,CBNT,CBNT,NRows,NCols,NCols,0.5,2.0,A,DC,Omm);

	pass = 0;
	if (array_norm_diff_d(NRows*NCols,Odiag,Omm,"Inf") < EPS)
		pass = 1;

	test_print2(pass,"                 ('R','R'):");

	// (L)eft, (C)ol Major
	mm_diag_d(NRows,NCols,dR,A,Odiag,2.0,1.5,'L','C');
	mm_d(CBCM,CBNT,CBNT,NRows,NCols,NRows,2.0,1.5,DR,A,Omm);

	pass = 0;
	if (array_norm_diff_d(NRows*NCols,Odiag,Omm,"Inf") < EPS)
		pass = 1;

	test_print2(pass,"                 ('L','C'):");

	// (R)ight, (C)ol Major
	mm_diag_d(NRows,NCols,dC,A,Odiag,0.5,0.2,'R','C');
	mm_d(CBCM,CBNT,CBNT,NRows,NCols,NCols,0.5,0.2,A,DC,Omm);

	pass = 0;
	if (array_norm_diff_d(NRows*NCols,Odiag,Omm,"Inf") < EPS)
		pass = 1;

	test_print2(pass,"                 ('R','C'):");

	free(DR);
	free(DC);
	free(Odiag);
	free(Omm);
}

static void test_unit_convert_to_CSR(void)
{
	unsigned int pass;
	unsigned int NRows = 5, NCols = 5,
	             rowIndex55[6] = {0, 3, 5, 8, 11, 13},
	             columns55[13] = {0, 1, 3, 0, 1, 2, 3, 4, 0, 2, 3, 1, 4};
	double       values55[13]  = {1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 2.0, 7.0, 8.0, -5.0},
	             A[25]  = { 1.0, -1.0,  0.0, -3.0,  0.0,
	                       -2.0,  5.0,  0.0,  0.0,  0.0,
	                        0.0,  0.0,  4.0,  6.0,  4.0,
	                       -4.0,  0.0,  2.0,  7.0,  0.0,
	                        0.0,  8.0,  0.0,  0.0, -5.0};

	struct S_OpCSR *A_sp;

	convert_to_CSR_d(NRows,NCols,A,&A_sp);

	pass = 0;
	if (array_norm_diff_ui(6,A_sp->rowIndex,rowIndex55,"Inf") < EPS &&
	    array_norm_diff_ui(13,A_sp->columns,columns55,"Inf")  < EPS &&
	    array_norm_diff_d(13,A_sp->values,values55,"Inf")     < EPS)
			pass = 1;

	test_print2(pass,"convert_to_CSR_d:");

	array_free1_CSR_d(A_sp);
}

void test_unit_matrix_functions(void)
{
	test_unit_matrix_diag();
	test_unit_matrix_identity();
	test_unit_matrix_inverse();
	test_unit_matrix_mm();
	test_unit_matrix_mm_diag();
	test_unit_convert_to_CSR();
}
