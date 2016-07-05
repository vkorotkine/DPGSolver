// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>

#include "test.h"
#include "parameters.h"
#include "functions.h"

#include "mkl.h"

/*
 *	Purpose:
 *		Test correctness of implementation of mm.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void test_imp_matrix_mm(void)
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
	mm_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,3,2,4,1.0,A,B,C_c);

	pass = 0;
	if (array_norm_diff_d(6,C,C_c,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("matrix_mm_d (RowMajor, NoTrans, NoTrans):        ");
	test_print(pass);


	// Trans, NoTrans
	mm_d(CblasRowMajor,CblasTrans,CblasNoTrans,3,2,4,1.0,AT,B,C_c);

	pass = 0;
	if (array_norm_diff_d(6,C,C_c,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("            (RowMajor,   Trans, NoTrans):        ");
	test_print(pass);


	// NoTrans, Trans
	mm_d(CblasRowMajor,CblasNoTrans,CblasTrans,3,2,4,1.0,A,BT,C_c);

	pass = 0;
	if (array_norm_diff_d(6,C,C_c,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("            (RowMajor, NoTrans,   Trans):        ");
	test_print(pass);


	// Trans, Trans
	mm_d(CblasRowMajor,CblasTrans,CblasTrans,3,2,4,1.0,AT,BT,C_c);

	pass = 0;
	if (array_norm_diff_d(6,C,C_c,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("            (RowMajor,   Trans,   Trans):        ");
	test_print(pass);



	// Col-Major

	// NoTrans, NoTrans
	mm_d(CblasColMajor,CblasNoTrans,CblasNoTrans,3,2,4,1.0,AT,BT,CT_c);

	pass = 0;
	if (array_norm_diff_d(6,CT,CT_c,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("            (ColMajor, NoTrans, NoTrans):        ");
	test_print(pass);


	// Trans, NoTrans
	mm_d(CblasColMajor,CblasTrans,CblasNoTrans,3,2,4,1.0,A,BT,CT_c);

	pass = 0;
	if (array_norm_diff_d(6,CT,CT_c,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("            (ColMajor,   Trans, NoTrans):        ");
	test_print(pass);


	// NoTrans, Trans
	mm_d(CblasColMajor,CblasNoTrans,CblasTrans,3,2,4,1.0,AT,B,CT_c);

	pass = 0;
	if (array_norm_diff_d(6,CT,CT_c,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("            (ColMajor, NoTrans,   Trans):        ");
	test_print(pass);


	// Trans, Trans
	mm_d(CblasColMajor,CblasTrans,CblasTrans,3,2,4,1.0,A,B,CT_c);

	pass = 0;
	if (array_norm_diff_d(6,CT,CT_c,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("            (ColMajor,   Trans,   Trans):        ");
	test_print(pass);

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
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("matrix_mm_Alloc_d (RowMajor):                    ");
	test_print(pass);

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
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("                  (ColMajor):                    ");
	test_print(pass);



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
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("matrix_mm_CTN_d (mv, NoBLAS):                    ");
	test_print(pass);

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
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("                (mv,   BLAS):                    ");
	test_print(pass);

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
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("                (mm, NoBLAS):                    ");
	test_print(pass);

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
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("                (mm,   BLAS):                    ");
	test_print(pass);

	free(A_CTN);
	free(B_CTN);
	free(C_CTN);
	free(C_CTN_c);
}
