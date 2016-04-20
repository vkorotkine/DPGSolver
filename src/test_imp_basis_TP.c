#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "test.h"
#include "functions.h"
#include "parameters.h"

/*
 *	Purpose:
 *		Test correctness of implementation of basis_TP.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

static double *basis_TP13(const double *rst, const unsigned int Nn)
{
	// ChiRef_rst = @(r) [sqrt(1/2)*1 sqrt(3/2)*r sqrt(5/2)*1/2*(3*r^2-1) sqrt(7/2)*1/2*(5*r^3-3*r)]

	unsigned int i, N, Nbf;
	double *ChiRef_rst, *r, r_i;

	N = 3+1;
	Nbf = pow(N,1);

	r = malloc(Nn * sizeof *r); // free
	for (i = 0; i < Nn; i++)
		r[i] = rst[0*Nn+i];

	ChiRef_rst = malloc(Nn*Nbf * sizeof *ChiRef_rst); // keep (requires external free)

	for (i = 0; i < Nn; i++) {
		r_i = r[i];

		ChiRef_rst[i*Nbf+0] = sqrt(1.0/2.0)*1.0;
		ChiRef_rst[i*Nbf+1] = sqrt(3.0/2.0)*r_i;
		ChiRef_rst[i*Nbf+2] = sqrt(5.0/2.0)*1.0/2.0*(3.0*pow(r_i,2.0)-1);
		ChiRef_rst[i*Nbf+3] = sqrt(7.0/2.0)*1.0/2.0*(5.0*pow(r_i,3.0)-3.0*r_i);
	}

	free(r);

	return ChiRef_rst;
}

static double *basis_TP22(const double *rst, const unsigned int Nn)
{
	/*	ChiRef_rst = @(r,s) [sqrt(1/2)*1*sqrt(1/2)*1
	 *	                     sqrt(3/2)*r*sqrt(1/2)*1
	 *	                     sqrt(5/2)*1/2*(3*r^2-1)*sqrt(1/2)*1
	 *	                     sqrt(1/2)*1*sqrt(3/2)*s
	 *	                     sqrt(3/2)*r*sqrt(3/2)*s
	 *	                     sqrt(5/2)*1/2*(3*r^2-1)*sqrt(3/2)*s
	 *	                     sqrt(1/2)*1*sqrt(5/2)*1/2*(3*s^2-1)
	 *	                     sqrt(3/2)*r*sqrt(5/2)*1/2*(3*s^2-1)
	 *	                     sqrt(5/2)*1/2*(3*r^2-1)*sqrt(5/2)*1/2*(3*s^2-1)]
	 */

	unsigned int i, N, Nbf;
	double *ChiRef_rst, *r, *s, r_i, s_i;

	N = 2+1;
	Nbf = pow(N,2);

	r = malloc(Nn * sizeof *r); // free
	s = malloc(Nn * sizeof *s); // free
	for (i = 0; i < Nn; i++) {
		r[i] = rst[0*Nn+i];
		s[i] = rst[1*Nn+i];
	}

	ChiRef_rst = malloc(Nn*Nbf * sizeof *ChiRef_rst); // keep (requires external free)

	for (i = 0; i < Nn; i++) {
		r_i = r[i];
		s_i = s[i];

		ChiRef_rst[i*Nbf+0] = sqrt(1.0/2.0)*1.0*sqrt(1.0/2.0)*1.0;
		ChiRef_rst[i*Nbf+1] = sqrt(3.0/2.0)*r_i*sqrt(1.0/2.0)*1.0;
		ChiRef_rst[i*Nbf+2] = sqrt(5.0/2.0)*1.0/2.0*(3.0*pow(r_i,2.0)-1.0)*sqrt(1.0/2.0)*1.0;
		ChiRef_rst[i*Nbf+3] = sqrt(1.0/2.0)*1.0*sqrt(3.0/2.0)*s_i;
		ChiRef_rst[i*Nbf+4] = sqrt(3.0/2.0)*r_i*sqrt(3.0/2.0)*s_i;
		ChiRef_rst[i*Nbf+5] = sqrt(5.0/2.0)*1.0/2.0*(3.0*pow(r_i,2.0)-1.0)*sqrt(3.0/2.0)*s_i;
		ChiRef_rst[i*Nbf+6] = sqrt(1.0/2.0)*1.0*sqrt(5.0/2.0)*1.0/2.0*(3.0*pow(s_i,2.0)-1.0);
		ChiRef_rst[i*Nbf+7] = sqrt(3.0/2.0)*r_i*sqrt(5.0/2.0)*1.0/2.0*(3.0*pow(s_i,2.0)-1.0);
		ChiRef_rst[i*Nbf+8] = sqrt(5.0/2.0)*1.0/2.0*(3.0*pow(r_i,2.0)-1.0)*sqrt(5.0/2.0)*1.0/2.0*(3.0*pow(s_i,2.0)-1.0);
	}

	free(r);
	free(s);

	return ChiRef_rst;
}

static double *basis_TP31(const double *rst, const unsigned int Nn)
{
	/*	ChiRef_rst = @(r,s,t) [sqrt(1/2)*1*sqrt(1/2)*1*sqrt(1/2)*1
	 *	                       sqrt(3/2)*r*sqrt(1/2)*1*sqrt(1/2)*1
	 *	                       sqrt(1/2)*1*sqrt(3/2)*s*sqrt(1/2)*1
	 *	                       sqrt(3/2)*r*sqrt(3/2)*s*sqrt(1/2)*1
	 *	                       sqrt(1/2)*1*sqrt(1/2)*1*sqrt(3/2)*t
	 *	                       sqrt(3/2)*r*sqrt(1/2)*1*sqrt(3/2)*t
	 *	                       sqrt(1/2)*1*sqrt(3/2)*s*sqrt(3/2)*t
	 *	                       sqrt(3/2)*r*sqrt(3/2)*s*sqrt(3/2)*t]
	 */

	unsigned int i, N, Nbf;
	double *ChiRef_rst, *r, *s, *t, r_i, s_i, t_i;

	N = 1+1;
	Nbf = pow(N,3);

	r = malloc(Nn * sizeof *r); // free
	s = malloc(Nn * sizeof *s); // free
	t = malloc(Nn * sizeof *t); // free
	for (i = 0; i < Nn; i++) {
		r[i] = rst[0*Nn+i];
		s[i] = rst[1*Nn+i];
		t[i] = rst[2*Nn+i];
	}

	ChiRef_rst = malloc(Nn*Nbf * sizeof *ChiRef_rst); // keep (requires external free)

	for (i = 0; i < Nn; i++) {
		r_i = r[i];
		s_i = s[i];
		t_i = t[i];

		ChiRef_rst[i*Nbf+0] = sqrt(1.0/2.0)*1.0*sqrt(1.0/2.0)*1.0*sqrt(1.0/2.0)*1.0;
		ChiRef_rst[i*Nbf+1] = sqrt(3.0/2.0)*r_i*sqrt(1.0/2.0)*1.0*sqrt(1.0/2.0)*1.0;
		ChiRef_rst[i*Nbf+2] = sqrt(1.0/2.0)*1.0*sqrt(3.0/2.0)*s_i*sqrt(1.0/2.0)*1.0;
		ChiRef_rst[i*Nbf+3] = sqrt(3.0/2.0)*r_i*sqrt(3.0/2.0)*s_i*sqrt(1.0/2.0)*1.0;
		ChiRef_rst[i*Nbf+4] = sqrt(1.0/2.0)*1.0*sqrt(1.0/2.0)*1.0*sqrt(3.0/2.0)*t_i;
		ChiRef_rst[i*Nbf+5] = sqrt(3.0/2.0)*r_i*sqrt(1.0/2.0)*1.0*sqrt(3.0/2.0)*t_i;
		ChiRef_rst[i*Nbf+6] = sqrt(1.0/2.0)*1.0*sqrt(3.0/2.0)*s_i*sqrt(3.0/2.0)*t_i;
		ChiRef_rst[i*Nbf+7] = sqrt(3.0/2.0)*r_i*sqrt(3.0/2.0)*s_i*sqrt(3.0/2.0)*t_i;
	}

	free(r);
	free(s);
	free(t);

	return ChiRef_rst;
}

void test_imp_basis_TP(void)
{
	unsigned int pass;

	/*
	 *	basis_TP (dE = 1):
	 *
	 *		Input:
	 *
	 *			P, rst, Nn, dE.
	 *
	 *		Expected output:
	 *
	 *			P = 3:
	 *				ChiRef_rst = @(r) [sqrt(1/2)*1 sqrt(3/2)*r sqrt(5/2)*1/2*(3*r^2-1) sqrt(7/2)*1/2*(5*r^3-3*r)]
	 */

	unsigned int dE, Nn, Nbf, P, Prst;
	double *rst, *w;

	dE = 1;

	double *ChiRef13_code, *ChiRef13_test;

	P = 3;
	Prst = 4;

	cubature_TP(&rst,&w,&Nn,0,Prst,dE,"GL"); // free

	ChiRef13_code = basis_TP(P,rst,Nn,&Nbf,dE); // free
	ChiRef13_test = basis_TP13(rst,Nn); // free

	pass = 0;
	if (array_norm_diff_d(Nn*(P+1),ChiRef13_code,ChiRef13_test,"Inf") < EPS*10)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("basis_TP (d1, P3):                               ");
	test_print(pass);

	free(rst);
	free(ChiRef13_code);
	free(ChiRef13_test);

	/*
	 *	basis_TP (dE = 2):
	 *
	 *		Input:
	 *
	 *			P, rst, Nn, dE.
	 *
	 *		Expected output:
	 *
	 *			P = 2:
	 *				ChiRef_rst = @(r,s) [sqrt(1/2)*1*sqrt(1/2)*1
	 *				                     sqrt(3/2)*r*sqrt(1/2)*1
	 *				                     sqrt(5/2)*1/2*(3*r^2-1)*sqrt(1/2)*1
	 *				                     sqrt(1/2)*1*sqrt(3/2)*s
	 *				                     sqrt(3/2)*r*sqrt(3/2)*s
	 *				                     sqrt(5/2)*1/2*(3*r^2-1)*sqrt(3/2)*s
	 *				                     sqrt(1/2)*1*sqrt(5/2)*1/2*(3*s^2-1)
	 *				                     sqrt(3/2)*r*sqrt(5/2)*1/2*(3*s^2-1)
	 *				                     sqrt(5/2)*1/2*(3*r^2-1)*sqrt(5/2)*1/2*(3*s^2-1)]
	 */

	dE = 2;

	double *ChiRef22_code, *ChiRef22_test;

	P = 2;
	Prst = 4;

	cubature_TP(&rst,&w,&Nn,0,Prst,dE,"GL"); // free

	ChiRef22_code = basis_TP(P,rst,Nn,&Nbf,dE); // free
	ChiRef22_test = basis_TP22(rst,Nn); // free

	pass = 0;
	if (array_norm_diff_d(Nn*pow(P+1,dE),ChiRef22_code,ChiRef22_test,"Inf") < EPS*10)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("         (d2, P2):                               ");
	test_print(pass);

	free(rst);
	free(ChiRef22_code);
	free(ChiRef22_test);

	/*
	 *	basis_TP (dE = 3):
	 *
	 *		Input:
	 *
	 *			P, rst, Nn, dE.
	 *
	 *		Expected output:
	 *
	 *			P = 1:
	 *				ChiRef_rst = @(r,s,t) [sqrt(1/2)*1*sqrt(1/2)*1*sqrt(1/2)*1
	 *				                       sqrt(3/2)*r*sqrt(1/2)*1*sqrt(1/2)*1
	 *				                       sqrt(1/2)*1*sqrt(3/2)*s*sqrt(1/2)*1
	 *				                       sqrt(3/2)*r*sqrt(3/2)*s*sqrt(1/2)*1
	 *				                       sqrt(1/2)*1*sqrt(1/2)*1*sqrt(3/2)*t
	 *				                       sqrt(3/2)*r*sqrt(1/2)*1*sqrt(3/2)*t
	 *				                       sqrt(1/2)*1*sqrt(3/2)*s*sqrt(3/2)*t
	 *				                       sqrt(3/2)*r*sqrt(3/2)*s*sqrt(3/2)*t]
	 */

	dE = 3;

	double *ChiRef31_code, *ChiRef31_test;

	P = 1;
	Prst = 4;

	cubature_TP(&rst,&w,&Nn,0,Prst,dE,"GL"); // free

	ChiRef31_code = basis_TP(P,rst,Nn,&Nbf,dE); // free
	ChiRef31_test = basis_TP31(rst,Nn); // free

	pass = 0;
	if (array_norm_diff_d(Nn*pow(P+1,dE),ChiRef31_code,ChiRef31_test,"Inf") < EPS*10)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("         (d3, P1):                               ");
	test_print(pass);

	free(rst);
	free(ChiRef31_code);
	free(ChiRef31_test);

	/*
	 *	basis_TP orthogonality:
	 *
	 *		Input:
	 *
	 *			rst_GL, rst_GLL
	 *
	 *		Expected Output:
	 *
	 *			M = ChiRef_rst'*W*ChiRef_rst = I (GL)
	 *			M = ChiRef_rst'*W*ChiRef_rst = I with error in highest order entries (GLL)
	 */

	double *ChiRef_rst, *W, *WChiRef_rst, *M, *I;

	// dE = 1
	dE = 1;

	P = 4;
	Nbf = pow(P+1,dE);
	I = identity_d(Nbf); // free

	// GL

	cubature_TP(&rst,&w,&Nn,1,P,dE,"GL"); // free

	W = diag_d(w,Nn); // free
	free(w);

	ChiRef_rst = basis_TP(P,rst,Nn,&Nbf,dE); // free

	WChiRef_rst = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,Nbf,Nn,1.0,W,ChiRef_rst); // free
	M = mm_Alloc_d(CblasRowMajor,CblasTrans,CblasNoTrans,Nbf,Nbf,Nn,1.0,ChiRef_rst,WChiRef_rst); // free

	pass = 0;
	if (array_norm_diff_d(pow(Nbf,2),M,I,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("         orthogonality - GL (d1):                ");
	test_print(pass);

	free(rst);
	free(W);
	free(ChiRef_rst);
	free(WChiRef_rst);
	free(M);

	// GLL

	cubature_TP(&rst,&w,&Nn,1,P,dE,"GLL"); // free

	W = diag_d(w,Nn); // free
	free(w);

	ChiRef_rst = basis_TP(P,rst,Nn,&Nbf,dE); // free

	WChiRef_rst = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,Nbf,Nn,1.0,W,ChiRef_rst); // free
	M = mm_Alloc_d(CblasRowMajor,CblasTrans,CblasNoTrans,Nbf,Nbf,Nn,1.0,ChiRef_rst,WChiRef_rst); // free

	pass = 0;
	if (array_norm_diff_d(pow(Nbf,2)-1,M,I,"Inf") < EPS*10)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("         orthogonality - GLL (d1):               ");
	test_print(pass);

	free(rst);
	free(W);
	free(ChiRef_rst);
	free(WChiRef_rst);
	free(M);

	free(I);

	// dE = 2
	dE = 2;

	P = 4;
	Nbf = pow(P+1,dE);
	I = identity_d(Nbf); // free

	// GL

	cubature_TP(&rst,&w,&Nn,1,P,dE,"GL"); // free

	W = diag_d(w,Nn); // free
	free(w);

	ChiRef_rst = basis_TP(P,rst,Nn,&Nbf,dE); // free

	WChiRef_rst = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,Nbf,Nn,1.0,W,ChiRef_rst); // free
	M = mm_Alloc_d(CblasRowMajor,CblasTrans,CblasNoTrans,Nbf,Nbf,Nn,1.0,ChiRef_rst,WChiRef_rst); // free

	pass = 0;
	if (array_norm_diff_d(pow(Nbf,2),M,I,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("         orthogonality - GL (d2):                ");
	test_print(pass);

	free(rst);
	free(W);
	free(ChiRef_rst);
	free(WChiRef_rst);
	free(M);

	free(I);

	// dE = 3
	dE = 3;

	P = 4;
	Nbf = pow(P+1,dE);
	I = identity_d(Nbf); // free

	// GL

	cubature_TP(&rst,&w,&Nn,1,P,dE,"GL"); // free

	W = diag_d(w,Nn); // free
	free(w);

	ChiRef_rst = basis_TP(P,rst,Nn,&Nbf,dE); // free

	WChiRef_rst = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,Nbf,Nn,1.0,W,ChiRef_rst); // free
	M = mm_Alloc_d(CblasRowMajor,CblasTrans,CblasNoTrans,Nbf,Nbf,Nn,1.0,ChiRef_rst,WChiRef_rst); // free

	pass = 0;
	if (array_norm_diff_d(pow(Nbf,2),M,I,"Inf") < EPS*10)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("         orthogonality - GL (d3):                ");
	test_print(pass);

	free(rst);
	free(W);
	free(ChiRef_rst);
	free(WChiRef_rst);
	free(M);

	free(I);
}
