// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

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
 *		basis_PYR2 tests for the modal basis in Chan(2015) (See basis_PYR's "References" section).
 *		The WV and WVHToP nodes are not currently tested for their integration strength as their use was giving
 *		unexpected results. See comments in cubature_PYR for details.
 *
 *	Notation:
 *
 *	References:
 */

static double *basis_PYR2(const double *rst, const unsigned int Nn)
{
	unsigned int n, d, Nbf;
	double       *ChiRef_rst, *a, *b, *c, a_n, b_n, c_n, con, con_c, con_i, con_j, con_k;

	d = 3;

	// Convert from rst to abc coordinates
	a = calloc(Nn , sizeof *a); // free
	b = calloc(Nn , sizeof *b); // free
	c = calloc(Nn , sizeof *c); // free

	rst_to_abc_PYR(Nn,d,rst,a,b,c);

	Nbf = 14;

	ChiRef_rst = malloc(Nn*Nbf * sizeof *ChiRef_rst); // keep (requires external free)

	con = pow(2.0,1.25);
	for (n = 0; n < Nn; n++) {
		a_n = a[n];
		b_n = b[n];
		c_n = c[n];

		get_scaling_basis_PYR(0,0,0,c_n,&con_i,&con_j,&con_k,&con_c);
		ChiRef_rst[n*Nbf+0] = con*con_c
		                    * con_i*(1.0)
		                    * con_j*(1.0)
		                    * con_k*(1.0);

		get_scaling_basis_PYR(0,0,1,c_n,&con_i,&con_j,&con_k,&con_c);
		ChiRef_rst[n*Nbf+1] = con*con_c
		                    * con_i*(1.0)
		                    * con_j*(1.0)
		                    * con_k*(1.0/2.0*(4.0*c_n+2.0));

		get_scaling_basis_PYR(0,0,2,c_n,&con_i,&con_j,&con_k,&con_c);
		ChiRef_rst[n*Nbf+2] = con*con_c
		                    * con_i*(1.0)
		                    * con_j*(1.0)
		                    * con_k*(15.0/4.0*pow(c_n-1.0,2.0)+10.0*(c_n-1.0)+6.0);

		get_scaling_basis_PYR(0,1,0,c_n,&con_i,&con_j,&con_k,&con_c);
		ChiRef_rst[n*Nbf+3] = con*con_c
		                    * con_i*(1.0)
		                    * con_j*(b_n)
		                    * con_k*(1.0);

		get_scaling_basis_PYR(0,1,1,c_n,&con_i,&con_j,&con_k,&con_c);
		ChiRef_rst[n*Nbf+4] = con*con_c
		                    * con_i*(1.0)
		                    * con_j*(b_n)
		                    * con_k*(1.0/2.0*(6.0*c_n+4.0));

		get_scaling_basis_PYR(0,2,0,c_n,&con_i,&con_j,&con_k,&con_c);
		ChiRef_rst[n*Nbf+5] = con*con_c
		                    * con_i*(1.0)
		                    * con_j*(1.0/2.0*(3.0*pow(b_n,2.0)-1.0))
		                    * con_k*(1.0);

		get_scaling_basis_PYR(1,0,0,c_n,&con_i,&con_j,&con_k,&con_c);
		ChiRef_rst[n*Nbf+6] = con*con_c
		                    * con_i*(a_n)
		                    * con_j*(1.0)
		                    * con_k*(1.0);

		get_scaling_basis_PYR(1,0,1,c_n,&con_i,&con_j,&con_k,&con_c);
		ChiRef_rst[n*Nbf+7] = con*con_c
		                    * con_i*(a_n)
		                    * con_j*(1.0)
		                    * con_k*(1.0/2.0*(6.0*c_n+4.0));

		get_scaling_basis_PYR(1,1,0,c_n,&con_i,&con_j,&con_k,&con_c);
		ChiRef_rst[n*Nbf+8] = con*con_c
		                    * con_i*(a_n)
		                    * con_j*(b_n)
		                    * con_k*(1.0);

		get_scaling_basis_PYR(1,1,1,c_n,&con_i,&con_j,&con_k,&con_c);
		ChiRef_rst[n*Nbf+9] = con*con_c
		                    * con_i*(a_n)
		                    * con_j*(b_n)
		                    * con_k*(1.0/2.0*(6.0*c_n+4.0));

		get_scaling_basis_PYR(1,2,0,c_n,&con_i,&con_j,&con_k,&con_c);
		ChiRef_rst[n*Nbf+10] = con*con_c
		                     * con_i*(a_n)
		                     * con_j*(1.0/2.0*(3.0*pow(b_n,2.0)-1.0))
		                     * con_k*(1.0);

		get_scaling_basis_PYR(2,0,0,c_n,&con_i,&con_j,&con_k,&con_c);
		ChiRef_rst[n*Nbf+11] = con*con_c
		                     * con_i*(1.0/2.0*(3.0*pow(a_n,2.0)-1.0))
		                     * con_j*(1.0)
		                     * con_k*(1.0);

		get_scaling_basis_PYR(2,1,0,c_n,&con_i,&con_j,&con_k,&con_c);
		ChiRef_rst[n*Nbf+12] = con*con_c
		                     * con_i*(1.0/2.0*(3.0*pow(a_n,2.0)-1.0))
		                     * con_j*(b_n)
		                     * con_k*(1.0);

		get_scaling_basis_PYR(2,2,0,c_n,&con_i,&con_j,&con_k,&con_c);
		ChiRef_rst[n*Nbf+13] = con*con_c
		                     * con_i*(1.0/2.0*(3.0*pow(a_n,2.0)-1.0))
		                     * con_j*(1.0/2.0*(3.0*pow(b_n,2.0)-1.0))
		                     * con_k*(1.0);
	}

	free(a);
	free(b);
	free(c);

	return ChiRef_rst;
}

void test_imp_basis_PYR(void)
{
	unsigned int pass;

	/*
	 *	basis_PYR
	 *
	 *		Input:
	 *
	 *			P, rst, Nn.
	 *
	 *		Expected output:
	 *
	 *			P = 2:
	 *				ChiRef_rst = @(r,s,t) [ See basis_PYR2 ]
	 */
	unsigned int d, Nn, Ns, Nbf, P, Prst;
	unsigned int *symms;
	double *rst, *w;

	d = 3;

	double *ChiRef32_code, *ChiRef32_test;

	P = 2;
	Prst = 4;

	cubature_PYR(&rst,&w,&symms,&Nn,&Ns,0,Prst,d,"WV"); // free

	ChiRef32_code = basis_PYR(P,rst,Nn,&Nbf,d); // free
	ChiRef32_test = basis_PYR2(rst,Nn);         // free

	pass = 0;
	if (array_norm_diff_d(Nn*Nbf,ChiRef32_code,ChiRef32_test,"Inf") < EPS*10)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("basis_PYR (P2):                                  ");
	test_print(pass);

	free(rst);
	free(symms);
	free(ChiRef32_code);
	free(ChiRef32_test);

	/*
	 *	basis_PYR orthogonality:
	 *
	 *		Input:
	 *			rst_GL/GLL, P, Prst
	 *
	 *		Expected Output:
	 *
	 *			P = 3, rst_GL/GLL (P = P+1/P+2):
	 *				M = ChiRef_rst'*W*ChiRef_rst = I
	 */

	double *ChiRef_rst, *W, *WChiRef_rst, *M, *Imat; // ToBeModified (name)

	d = 3;

	// P = 3
	P = 3;
	Nbf = 30;
	Imat = identity_d(Nbf); // free

	// GL (HEX To PYR)
	Prst = P+1;
	cubature_PYR(&rst,&w,&symms,&Nn,&Ns,1,Prst,d,"GLW"); // free

	W = diag_d(w,Nn); // free
	free(w);

	ChiRef_rst = basis_PYR(P,rst,Nn,&Nbf,d); // free

	WChiRef_rst = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,Nbf,Nn,1.0,W,ChiRef_rst);          // free
	M           = mm_Alloc_d(CblasRowMajor,CblasTrans,CblasNoTrans,Nbf,Nbf,Nn,1.0,ChiRef_rst,WChiRef_rst); // free

	pass = 0;
	if (array_norm_diff_d(pow(Nbf,2),M,Imat,"Inf") < EPS*1e3)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("          orthogonality - GL  (P3):              ");
	test_print(pass);

	free(rst);
	free(symms);
	free(W);
	free(ChiRef_rst);
	free(WChiRef_rst);
	free(M);

	// GLL (HEX To PYR)
	Prst = P+2;
	cubature_PYR(&rst,&w,&symms,&Nn,&Ns,1,Prst,d,"GLW"); // free

	W = diag_d(w,Nn); // free
	free(w);

	ChiRef_rst = basis_PYR(P,rst,Nn,&Nbf,d); // free

	WChiRef_rst = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,Nbf,Nn,1.0,W,ChiRef_rst);          // free
	M           = mm_Alloc_d(CblasRowMajor,CblasTrans,CblasNoTrans,Nbf,Nbf,Nn,1.0,ChiRef_rst,WChiRef_rst); // free

	pass = 0;
	if (array_norm_diff_d(pow(Nbf,2),M,Imat,"Inf") < EPS*1e2)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("                          GLL (P3):              ");
	test_print(pass);

	free(rst);
	free(symms);
	free(W);
	free(ChiRef_rst);
	free(WChiRef_rst);
	free(M);

	free(Imat);
}
