// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "test_unit_matrix_functions.h"

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

void test_unit_matrix_diag(void)
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
	       X4[16] = {1.0, 0.0, 0.0 , 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 4.0};

	X = diag_d(x,N); // free

	pass = 0;
	if (array_norm_diff_d(16,X,X4,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("matrix_diag_d:                                   ");
	test_print(pass);

	free(X);
}

void test_unit_matrix_identity(void)
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
	double *Imat, // ToBeModified (name)
	       I4[16] = {1.0, 0.0, 0.0 , 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};

	Imat = identity_d(N); // free

	pass = 0;
	if (array_norm_diff_d(16,Imat,I4,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("matrix_identity_d:                               ");
	test_print(pass);

	free(Imat);
}

void test_unit_matrix_inverse(void)
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
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("matrix_inverse_d:                                ");
	test_print(pass);

	free(I3);
	free(AInv_c);
	free(Aroundtrip_c);
}

void test_unit_matrix_mm(void)
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

void test_unit_convert_to_CSR(void)
{
	/*
	 *	Purpose:
	 *		Test correctness of implementation of convert_to_CSR.
	 *
	 *	Comments:
	 *
	 *	Notation:
	 *
	 *	References:
	 */

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
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("convert_to_CSR_d:                                ");
	test_print(pass);

	array_free1_CSR_d(A_sp);
}
