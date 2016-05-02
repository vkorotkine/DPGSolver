#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "test.h"
#include "functions.h"
#include "parameters.h"

/*
 *	Purpose:
 *		Test correctness of implementation of basis_PYR.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

/*
static double *basis_TRI2(const double *rst, const unsigned int Nn)
{
	unsigned int n, d, Nbf;
	double       *ChiRef_rst, *a, *b, *c, a_n, b_n, con, con_b, con_i, con_j;

	d = 2;

	// Convert from rst to abc coordinates
	a = calloc(Nn , sizeof *a); // free
	b = calloc(Nn , sizeof *b); // free
	c = calloc(Nn , sizeof *c); // free

	rst_to_abc_SI(Nn,d,rst,a,b,c);

	Nbf = 6;

	ChiRef_rst = malloc(Nn*Nbf * sizeof *ChiRef_rst); // keep (requires external free)

	con = 2.0/pow(3.0,0.25);
	for (n = 0; n < Nn; n++) {
		a_n = a[n];
		b_n = b[n];

		get_scaling_basis_TRI(0,0,b_n,&con_i,&con_j,&con_b);
		ChiRef_rst[n*Nbf+0] = con*con_b
		                    * con_i*1.0
		                    * con_j*1.0;

		get_scaling_basis_TRI(0,1,b_n,&con_i,&con_j,&con_b);
		ChiRef_rst[n*Nbf+1] = con*con_b
		                    * con_i*1.0
		                    * con_j*1.0/2.0*(3.0*b_n+1.0);

		get_scaling_basis_TRI(0,2,b_n,&con_i,&con_j,&con_b);
		ChiRef_rst[n*Nbf+2] = con*con_b
		                    * con_i*1.0
		                    * con_j*(5.0/2.0*pow(b_n-1.0,2.0)+6.0*(b_n-1.0)+3.0);

		get_scaling_basis_TRI(1,0,b_n,&con_i,&con_j,&con_b);
		ChiRef_rst[n*Nbf+3] = con*con_b
		                    * con_i*a_n
		                    * con_j*1.0;

		get_scaling_basis_TRI(1,1,b_n,&con_i,&con_j,&con_b);
		ChiRef_rst[n*Nbf+4] = con*con_b
		                    * con_i*a_n
		                    * con_j*1.0/2.0*(5.0*b_n+3.0);

		get_scaling_basis_TRI(2,0,b_n,&con_i,&con_j,&con_b);
		ChiRef_rst[n*Nbf+5] = con*con_b
		                    * con_i*1.0/2.0*(3.0*pow(a_n,2.0)-1.0)
		                    * con_j*1.0;
	}

	free(a);
	free(b);
	free(c);

	return ChiRef_rst;
}

static double *basis_TRI3(const double *rst, const unsigned int Nn)
{
	unsigned int n, d, Nbf;
	double       *ChiRef_rst, *a, *b, *c, a_n, b_n, con, con_b, con_i, con_j;

	d = 2;

	// Convert from rst to abc coordinates
	a = calloc(Nn , sizeof *a); // free
	b = calloc(Nn , sizeof *b); // free
	c = calloc(Nn , sizeof *c); // free

	rst_to_abc_SI(Nn,d,rst,a,b,c);

	Nbf = 10;

	ChiRef_rst = malloc(Nn*Nbf * sizeof *ChiRef_rst); // keep (requires external free)

	con = 2.0/pow(3.0,0.25);
	for (n = 0; n < Nn; n++) {
		a_n = a[n];
		b_n = b[n];

		get_scaling_basis_TRI(0,0,b_n,&con_i,&con_j,&con_b);
		ChiRef_rst[n*Nbf+0] = con*con_b
		                    * con_i*1.0
		                    * con_j*1.0;

		get_scaling_basis_TRI(0,1,b_n,&con_i,&con_j,&con_b);
		ChiRef_rst[n*Nbf+1] = con*con_b
		                    * con_i*1.0
		                    * con_j*1.0/2.0*(3.0*b_n+1.0);

		get_scaling_basis_TRI(0,2,b_n,&con_i,&con_j,&con_b);
		ChiRef_rst[n*Nbf+2] = con*con_b
		                    * con_i*1.0
		                    * con_j*(5.0/2.0*pow(b_n-1.0,2.0)+6.0*(b_n-1.0)+3.0);

		get_scaling_basis_TRI(0,3,b_n,&con_i,&con_j,&con_b);
		ChiRef_rst[n*Nbf+3] = con*con_b
		                    * con_i*1.0
		                    * con_j*(35.0/8.0*pow(b_n-1.0,3.0)+15.0*pow(b_n-1.0,2.0)+15.0*(b_n-1)+4.0);

		get_scaling_basis_TRI(1,0,b_n,&con_i,&con_j,&con_b);
		ChiRef_rst[n*Nbf+4] = con*con_b
		                    * con_i*a_n
		                    * con_j*1.0;

		get_scaling_basis_TRI(1,1,b_n,&con_i,&con_j,&con_b);
		ChiRef_rst[n*Nbf+5] = con*con_b
		                    * con_i*a_n
		                    * con_j*1.0/2.0*(5.0*b_n+3.0);

		get_scaling_basis_TRI(1,2,b_n,&con_i,&con_j,&con_b);
		ChiRef_rst[n*Nbf+6] = con*con_b
		                    * con_i*a_n
		                    * con_j*(21.0/4.0*pow(b_n-1.0,2.0)+15.0*(b_n-1.0)+10.0);

		get_scaling_basis_TRI(2,0,b_n,&con_i,&con_j,&con_b);
		ChiRef_rst[n*Nbf+7] = con*con_b
		                    * con_i*1.0/2.0*(3.0*pow(a_n,2.0)-1.0)
		                    * con_j*1.0;

		get_scaling_basis_TRI(2,1,b_n,&con_i,&con_j,&con_b);
		ChiRef_rst[n*Nbf+8] = con*con_b
		                    * con_i*1.0/2.0*(3.0*pow(a_n,2.0)-1.0)
		                    * con_j*1.0/2.0*(7.0*b_n+5.0);

		get_scaling_basis_TRI(3,0,b_n,&con_i,&con_j,&con_b);
		ChiRef_rst[n*Nbf+9] = con*con_b
		                    * con_i*1.0/2.0*(5.0*pow(a_n,3.0)-3.0*a_n)
		                    * con_j*1.0;
	}

	free(a);
	free(b);
	free(c);

	return ChiRef_rst;
}

static double *basis_TET2(const double *rst, const unsigned int Nn)
{
	unsigned int n, d, Nbf;
	double       *ChiRef_rst, *a, *b, *c, a_n, b_n, c_n, con, con_b, con_c, con_i, con_j, con_k;

	d = 3;

	// Convert from rst to abc coordinates
	a = calloc(Nn , sizeof *a); // free
	b = calloc(Nn , sizeof *b); // free
	c = calloc(Nn , sizeof *c); // free

	rst_to_abc_SI(Nn,d,rst,a,b,c);

	Nbf = 10;

	ChiRef_rst = malloc(Nn*Nbf * sizeof *ChiRef_rst); // keep (requires external free)

	con = 4.0/pow(2.0,0.25);
	for (n = 0; n < Nn; n++) {
		a_n = a[n];
		b_n = b[n];
		c_n = c[n];

		get_scaling_basis_TET(0,0,0,b_n,c_n,&con_i,&con_j,&con_k,&con_b,&con_c);
		ChiRef_rst[n*Nbf+0] = con*con_b*con_c
		                    * con_i*1.0
		                    * con_j*1.0
		                    * con_k*1.0;

		get_scaling_basis_TET(0,0,1,b_n,c_n,&con_i,&con_j,&con_k,&con_b,&con_c);
		ChiRef_rst[n*Nbf+1] = con*con_b*con_c
		                    * con_i*1.0
		                    * con_j*1.0
		                    * con_k*1.0/2.0*(4.0*c_n+2.0);

		get_scaling_basis_TET(0,0,2,b_n,c_n,&con_i,&con_j,&con_k,&con_b,&con_c);
		ChiRef_rst[n*Nbf+2] = con*con_b*con_c
		                    * con_i*1.0
		                    * con_j*1.0
		                    * con_k*(15.0/4.0*pow(c_n-1.0,2.0)+10.0*(c_n-1.0)+6.0);

		get_scaling_basis_TET(0,1,0,b_n,c_n,&con_i,&con_j,&con_k,&con_b,&con_c);
		ChiRef_rst[n*Nbf+3] = con*con_b*con_c
		                    * con_i*1.0
		                    * con_j*1.0/2.0*(3.0*b_n+1.0)
		                    * con_k*1.0;

		get_scaling_basis_TET(0,1,1,b_n,c_n,&con_i,&con_j,&con_k,&con_b,&con_c);
		ChiRef_rst[n*Nbf+4] = con*con_b*con_c
		                    * con_i*1.0
		                    * con_j*1.0/2.0*(3.0*b_n+1.0)
		                    * con_k*1.0/2.0*(6.0*c_n+4.0);

		get_scaling_basis_TET(0,2,0,b_n,c_n,&con_i,&con_j,&con_k,&con_b,&con_c);
		ChiRef_rst[n*Nbf+5] = con*con_b*con_c
		                    * con_i*1.0
		                    * con_j*(5.0/2.0*pow(b_n-1.0,2.0)+6.0*(b_n-1.0)+3.0)
		                    * con_k*1.0;

		get_scaling_basis_TET(1,0,0,b_n,c_n,&con_i,&con_j,&con_k,&con_b,&con_c);
		ChiRef_rst[n*Nbf+6] = con*con_b*con_c
		                    * con_i*a_n
		                    * con_j*1.0
		                    * con_k*1.0;

		get_scaling_basis_TET(1,0,1,b_n,c_n,&con_i,&con_j,&con_k,&con_b,&con_c);
		ChiRef_rst[n*Nbf+7] = con*con_b*con_c
		                    * con_i*a_n
		                    * con_j*1.0
		                    * con_k*1.0/2.0*(6.0*c_n+4.0);

		get_scaling_basis_TET(1,1,0,b_n,c_n,&con_i,&con_j,&con_k,&con_b,&con_c);
		ChiRef_rst[n*Nbf+8] = con*con_b*con_c
		                    * con_i*a_n
		                    * con_j*1.0/2.0*(5.0*b_n+3.0)
		                    * con_k*1.0;

		get_scaling_basis_TET(2,0,0,b_n,c_n,&con_i,&con_j,&con_k,&con_b,&con_c);
		ChiRef_rst[n*Nbf+9] = con*con_b*con_c
		                    * con_i*1.0/2.0*(3.0*pow(a_n,2.0)-1.0)
		                    * con_j*1.0
		                    * con_k*1.0;
	}

	free(a);
	free(b);
	free(c);

	return ChiRef_rst;
}
*/

void test_imp_basis_PYR(void)
{
	unsigned int pass;

	/*
	 *	basis_SI (d = 2):
	 *
	 *		Input:
	 *
	 *			P, rst, Nn.
	 *
	 *		Expected output:
	 *
	 *			P = 2:
	 *				ChiRef_rst = @(r,s) [ See basis_TRI2 ]
	 *
	 *			P = 3:
	 *				ChiRef_rst = @(r,s) [ See basis_TRI3 ]
	 */
	unsigned int d, Nn, Ns, Nbf, P, Prst;
	unsigned int *symms;
	double *rst, *w;
/*

	d = 2;

	double *ChiRef22_code, *ChiRef22_test;

	P = 2;
	Prst = 4;

	cubature_TRI(&rst,&w,&symms,&Nn,&Ns,0,Prst,d,"AO"); // free

	ChiRef22_code = basis_SI(P,rst,Nn,&Nbf,d); // free
	ChiRef22_test = basis_TRI2(rst,Nn);        // free

	pass = 0;
	if (array_norm_diff_d(Nn*Nbf,ChiRef22_code,ChiRef22_test,"Inf") < EPS*10)
		pass = 1, TestDB.Npass++;


	//     0         10        20        30        40        50
	printf("basis_SI (d2, P2):                               ");
	test_print(pass);

	free(rst);
	free(symms);
	free(ChiRef22_code);
	free(ChiRef22_test);

	double *ChiRef23_code, *ChiRef23_test;

	P = 3;
	Prst = 4;

	cubature_TRI(&rst,&w,&symms,&Nn,&Ns,0,Prst,d,"AO"); // free

	ChiRef23_code = basis_SI(P,rst,Nn,&Nbf,d); // free
	ChiRef23_test = basis_TRI3(rst,Nn);        // free

	pass = 0;
	if (array_norm_diff_d(Nn*Nbf,ChiRef23_code,ChiRef23_test,"Inf") < EPS*100)
		pass = 1, TestDB.Npass++;


	//     0         10        20        30        40        50
	printf("         (d2, P3):                               ");
	test_print(pass);

	free(rst);
	free(symms);
	free(ChiRef23_code);
	free(ChiRef23_test);
*/

	/*
	 *	basis_SI (d = 3):
	 *
	 *		Input:
	 *
	 *			P, rst, Nn.
	 *
	 *		Expected output:
	 *
	 *			P = 2:
	 *				ChiRef_rst = @(r,s,t) [ See basis_TET2 ]
	 */
/*

	d = 3;

	double *ChiRef32_code, *ChiRef32_test;

	P = 2;
	Prst = 4;

	cubature_TET(&rst,&w,&symms,&Nn,&Ns,0,Prst,d,"AO"); // free

	ChiRef32_code = basis_SI(P,rst,Nn,&Nbf,d); // free
	ChiRef32_test = basis_TET2(rst,Nn);        // free

	pass = 0;
	if (array_norm_diff_d(Nn*Nbf,ChiRef32_code,ChiRef32_test,"Inf") < EPS*10)
		pass = 1, TestDB.Npass++;


	//     0         10        20        30        40        50
	printf("         (d3, P2):                               ");
	test_print(pass);

	free(rst);
	free(symms);
	free(ChiRef32_code);
	free(ChiRef32_test);
*/

	/*
	 *	basis_SI orthogonality:
	 *
	 *		Input:
	 *			rst_WS (d = 2), rst_SH (d = 3), rst_WV, P, d
	 *
	 *		Expected Output:
	 *
	 *			d = 2, P = 2, rst_WS (P = 2), rst_WV (P = 4):
	 *				M = ChiRef_rst'*W*ChiRef_rst = I
	 *
	 *			d = 2, P = 3, rst_WS (P = 3), rst_WV (P = 6):
	 *				M = ChiRef_rst'*W*ChiRef_rst = I with errro in highest order entries (WS)
	 *				M = ChiRef_rst'*W*ChiRef_rst = I (WV)
	 *
	 *			d = 3, P = 1, rst_SH (P = 1), rst_WV (P = 2):
	 *				M = ChiRef_rst'*W*ChiRef_rst = I
	 *
	 *			d = 2, P = 2, rst_SH (P = 2), rst_WV (P = 4):
	 *				M = ChiRef_rst'*W*ChiRef_rst = I with errro in highest order entries (SH)
	 *				M = ChiRef_rst'*W*ChiRef_rst = I (WV)
	 */

	double *ChiRef_rst, *W, *WChiRef_rst, *M, *I;

	d = 3;

	// P = 2
	P = 2;
	Nbf = 14;
	I = identity_d(Nbf); // free

	// WV
	Prst = 8;
//	Prst = 10;
	cubature_PYR(&rst,&w,&symms,&Nn,&Ns,1,Prst,d,"WV"); // free

	W = diag_d(w,Nn); // free
	free(w);

	ChiRef_rst = basis_PYR(P,rst,Nn,&Nbf,d); // free

	WChiRef_rst = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,Nbf,Nn,1.0,W,ChiRef_rst);          // free
	M           = mm_Alloc_d(CblasRowMajor,CblasTrans,CblasNoTrans,Nbf,Nbf,Nn,1.0,ChiRef_rst,WChiRef_rst); // free

	pass = 0;
	if (array_norm_diff_d(pow(Nbf,2),M,I,"Inf") < EPS*1e2)
		pass = 1, TestDB.Npass++;

array_print_d(Nbf,Nbf,M,'R');
//array_print_d(Nbf,Nbf,I,'R');
for (int i = 0; i < Nbf; i++) {
for (int j = 0; j < Nbf; j++) {
	if (i == j)
		printf("%d % .12e\n",i,M[i*Nbf+j]);
}}
printf("%d\n",Nn);

	//     0         10        20        30        40        50
	printf("                         WV (d2, P4):            ");
		test_print(pass);

	free(rst);
	free(symms);
	free(W);
	free(ChiRef_rst);
	free(WChiRef_rst);
	free(M);

//	free(I);

	Prst = 3;
	cubature_TP(&rst,&w,&symms,&Nn,&Ns,1,Prst,d,"GL"); // free

	double rst_tmp[3], a, b, c;

	for (int i = 0; i < Nn; i++) {
		a = rst[0*Nn+i];
		b = rst[1*Nn+i];
		c = rst[2*Nn+i];

		rst_tmp[0] = 0.5*(1.0-c)*a;
		rst_tmp[1] = 0.5*(1.0-c)*b;
		rst_tmp[2] = sqrt(2.0)/2.0*(0.6+c);

		for (int j = 0; j < d; j++)
			rst[j*Nn+i] = rst_tmp[j];
		// Add integral Jacobian and scaling to weights
		w[i] = w[i]*pow(2.0,-2.5)*pow(1.0-c,2.0);
	}

	W = diag_d(w,Nn); // free
	free(w);

	ChiRef_rst = basis_PYR(P,rst,Nn,&Nbf,d); // free

	WChiRef_rst = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,Nbf,Nn,1.0,W,ChiRef_rst);          // free
	M           = mm_Alloc_d(CblasRowMajor,CblasTrans,CblasNoTrans,Nbf,Nbf,Nn,1.0,ChiRef_rst,WChiRef_rst); // free

	pass = 0;
	if (array_norm_diff_d(pow(Nbf,2),M,I,"Inf") < EPS*1e2)
		pass = 1, TestDB.Npass++;

array_print_d(Nbf,Nbf,M,'R');
//array_print_d(Nbf,Nbf,I,'R');
for (int i = 0; i < Nbf; i++) {
for (int j = 0; j < Nbf; j++) {
	if (i == j)
		printf("%d % .12e\n",i,M[i*Nbf+j]);
}}
printf("%d\n",Nn);

	//     0         10        20        30        40        50
	printf("                         TP (d3, P?):            ");
		test_print(pass);
exit(1);

	free(rst);
	free(symms);
	free(W);
	free(ChiRef_rst);
	free(WChiRef_rst);
	free(M);

	free(I);




/*	// P = 3
	P = 3;

	// WV
	Prst = 2*P;
	cubature_TRI(&rst,&w,&symms,&Nn,&Ns,1,Prst,d,"WV"); // free

	W = diag_d(w,Nn); // free
	free(w);

	ChiRef_rst = basis_SI(P,rst,Nn,&Nbf,d); // free

	WChiRef_rst = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,Nbf,Nn,1.0,W,ChiRef_rst);          // free
	M           = mm_Alloc_d(CblasRowMajor,CblasTrans,CblasNoTrans,Nbf,Nbf,Nn,1.0,ChiRef_rst,WChiRef_rst); // free

	pass = 0;
	if (array_norm_diff_d(pow(Nbf,2),M,I,"Inf") < EPS*1e2)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("                         WV (d2, P6):            ");
		test_print(pass);

	free(rst);
	free(symms);
	free(W);
	free(ChiRef_rst);
	free(WChiRef_rst);
	free(M);

	free(I);
*/
}
