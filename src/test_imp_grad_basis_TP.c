#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "test.h"
#include "functions.h"
#include "parameters.h"

/*
 *	Purpose:
 *		Test correctness of implementation of grad_basis_TP.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

static double **grad_basis_TP13(const double *rst, const unsigned int Nn)
{
	// GradChiRef_rst = @(r) [0 sqrt(3/2) sqrt(5/2)*3*r sqrt(7/2)*1/2*(15*r^2-3)]

	unsigned int i, N, Nbf, dim;
	double **GradChiRef_rst, *r, r_i;

	N = 3+1;
	Nbf = pow(N,1);

	r = malloc(Nn * sizeof *r); // free
	for (i = 0; i < Nn; i++)
		r[i] = rst[0*Nn+i];

	GradChiRef_rst = malloc(1 * sizeof *GradChiRef_rst); // keep (requires external free)
	for (dim = 0; dim < 1; dim++) {
		GradChiRef_rst[dim] = malloc(Nn*Nbf * sizeof **GradChiRef_rst); // keep (requires external free)
	}

	for (i = 0; i < Nn; i++) {
		r_i = r[i];

		GradChiRef_rst[0][i*Nbf+0] = 0.0;
		GradChiRef_rst[0][i*Nbf+1] = sqrt(3.0/2.0);
		GradChiRef_rst[0][i*Nbf+2] = sqrt(5.0/2.0)*3.0*r_i;
		GradChiRef_rst[0][i*Nbf+3] = sqrt(7.0/2.0)*1.0/2.0*(15.0*pow(r_i,2.0)-3.0);
	}

	free(r);

	return GradChiRef_rst;
}

static double **grad_basis_TP22(const double *rst, const unsigned int Nn)
{
	/*
	 *	GradChiRef_rst[0] = @(r,s) [0
	 *	                            sqrt(3/2)*sqrt(1/2)*1
	 *	                            sqrt(5/2)*3*r*sqrt(1/2)*1
	 *	                            0
	 *	                            sqrt(3/2)*sqrt(3/2)*s
	 *	                            sqrt(5/2)*3*r*sqrt(3/2)*s
	 *	                            0
	 *	                            sqrt(3/2)*sqrt(5/2)*1/2*(3*s^2-1)
	 *	                            sqrt(5/2)*3*r*sqrt(5/2)*1/2*(3*s^2-1)]
	 *
	 *	GradChiRef_rst[1] = @(r,s) [0
	 *	                            0
	 *	                            0
	 *	                            sqrt(1/2)*1*sqrt(3/2)
	 *	                            sqrt(3/2)*r*sqrt(3/2)
	 *	                            sqrt(5/2)*1/2*(3*r^2-1)*sqrt(3/2)
	 *	                            sqrt(1/2)*1*sqrt(5/2)*3*s
	 *	                            sqrt(3/2)*r*sqrt(5/2)*3*s
	 *	                            sqrt(5/2)*1/2*(3*r^2-1)*sqrt(5/2)*3*s]
	 */

	unsigned int i, N, Nbf, dim;
	double **GradChiRef_rst, *r, *s, r_i, s_i;

	N = 2+1;
	Nbf = pow(N,2);

	GradChiRef_rst = malloc(2 * sizeof *GradChiRef_rst); // keep (requires external free)
	for (dim = 0; dim < 2; dim++) {
		GradChiRef_rst[dim] = malloc(Nn*Nbf * sizeof **GradChiRef_rst); // keep (requires external free)
	}

	r = malloc(Nn * sizeof *r); // free
	s = malloc(Nn * sizeof *s); // free
	for (i = 0; i < Nn; i++) {
		r[i] = rst[0*Nn+i];
		s[i] = rst[1*Nn+i];
	}

	for (i = 0; i < Nn; i++) {
		r_i = r[i];
		s_i = s[i];

		GradChiRef_rst[0][i*Nbf+0] = 0.0;
		GradChiRef_rst[0][i*Nbf+1] = sqrt(3.0/2.0)*sqrt(1.0/2.0)*1.0;
		GradChiRef_rst[0][i*Nbf+2] = sqrt(5.0/2.0)*3.0*r_i*sqrt(1.0/2.0)*1.0;
		GradChiRef_rst[0][i*Nbf+3] = 0.0;
		GradChiRef_rst[0][i*Nbf+4] = sqrt(3.0/2.0)*sqrt(3.0/2.0)*s_i;
		GradChiRef_rst[0][i*Nbf+5] = sqrt(5.0/2.0)*3.0*r_i*sqrt(3.0/2.0)*s_i;
		GradChiRef_rst[0][i*Nbf+6] = 0.0;
		GradChiRef_rst[0][i*Nbf+7] = sqrt(3.0/2.0)*sqrt(5.0/2.0)*1.0/2.0*(3.0*pow(s_i,2.0)-1.0);
		GradChiRef_rst[0][i*Nbf+8] = sqrt(5.0/2.0)*3.0*r_i*sqrt(5.0/2.0)*1.0/2.0*(3.0*pow(s_i,2.0)-1.0);

		GradChiRef_rst[1][i*Nbf+0] = 0.0;
		GradChiRef_rst[1][i*Nbf+1] = 0.0;
		GradChiRef_rst[1][i*Nbf+2] = 0.0;
		GradChiRef_rst[1][i*Nbf+3] = sqrt(1.0/2.0)*1.0*sqrt(3.0/2.0);
		GradChiRef_rst[1][i*Nbf+4] = sqrt(3.0/2.0)*r_i*sqrt(3.0/2.0);
		GradChiRef_rst[1][i*Nbf+5] = sqrt(5.0/2.0)*1.0/2.0*(3.0*pow(r_i,2.0)-1.0)*sqrt(3.0/2.0);
		GradChiRef_rst[1][i*Nbf+6] = sqrt(1.0/2.0)*1.0*sqrt(5.0/2.0)*3.0*s_i;
		GradChiRef_rst[1][i*Nbf+7] = sqrt(3.0/2.0)*r_i*sqrt(5.0/2.0)*3.0*s_i;
		GradChiRef_rst[1][i*Nbf+8] = sqrt(5.0/2.0)*1.0/2.0*(3.0*pow(r_i,2.0)-1.0)*sqrt(5.0/2.0)*3.0*s_i;
	}

	free(r);
	free(s);

	return GradChiRef_rst;
}

static double **grad_basis_TP31(const double *rst, const unsigned int Nn)
{
	/*
	 *				GradChiRef_rst[0] = @(r,s,t) [0
	 *				                              sqrt(3/2)*sqrt(1/2)*1*sqrt(1/2)*1
	 *				                              0
	 *				                              sqrt(3/2)*sqrt(3/2)*1*sqrt(1/2)*1
	 *				                              0
	 *				                              sqrt(3/2)*sqrt(1/2)*1*sqrt(3/2)*t
	 *				                              0
	 *				                              sqrt(3/2)*sqrt(3/2)*s*sqrt(3/2)*t]
	 *
	 *				GradChiRef_rst[1] = @(r,s,t) [0
	 *				                              0
	 *				                              sqrt(1/2)*1*sqrt(3/2)*sqrt(1/2)*1
	 *				                              sqrt(3/2)*r*sqrt(3/2)*sqrt(1/2)*1
	 *				                              0
	 *				                              0
	 *				                              sqrt(1/2)*1*sqrt(3/2)*sqrt(3/2)*t
	 *				                              sqrt(3/2)*r*sqrt(3/2)*sqrt(3/2)*t
	 *
	 *				GradChiRef_rst[2] = @(r,s,t) [0
	 *				                              0
	 *				                              0
	 *				                              0
	 *				                              sqrt(1/2)*1*sqrt(1/2)*1*sqrt(3/2)
	 *				                              sqrt(3/2)*r*sqrt(1/2)*1*sqrt(3/2)
	 *				                              sqrt(1/2)*1*sqrt(3/2)*s*sqrt(3/2)
	 *				                              sqrt(3/2)*r*sqrt(3/2)*s*sqrt(3/2)]
	 */

	unsigned int i, N, Nbf, dim;
	double **GradChiRef_rst, *r, *s, *t, r_i, s_i, t_i;

	N = 1+1;
	Nbf = pow(N,3);

	GradChiRef_rst = malloc(3 * sizeof *GradChiRef_rst); // keep (requires external free)
	for (dim = 0; dim < 3; dim++) {
		GradChiRef_rst[dim] = malloc(Nn*Nbf * sizeof **GradChiRef_rst); // keep (requires external free)
	}

	r = malloc(Nn * sizeof *r); // free
	s = malloc(Nn * sizeof *s); // free
	t = malloc(Nn * sizeof *t); // free
	for (i = 0; i < Nn; i++) {
		r[i] = rst[0*Nn+i];
		s[i] = rst[1*Nn+i];
		t[i] = rst[2*Nn+i];
	}

	for (i = 0; i < Nn; i++) {
		r_i = r[i];
		s_i = s[i];
		t_i = t[i];

		GradChiRef_rst[0][i*Nbf+0] = sqrt(1.0/2.0)*0.0*sqrt(1.0/2.0)*1.0*sqrt(1.0/2.0)*1.0;
		GradChiRef_rst[0][i*Nbf+1] = sqrt(3.0/2.0)*1.0*sqrt(1.0/2.0)*1.0*sqrt(1.0/2.0)*1.0;
		GradChiRef_rst[0][i*Nbf+2] = sqrt(1.0/2.0)*0.0*sqrt(3.0/2.0)*s_i*sqrt(1.0/2.0)*1.0;
		GradChiRef_rst[0][i*Nbf+3] = sqrt(3.0/2.0)*1.0*sqrt(3.0/2.0)*s_i*sqrt(1.0/2.0)*1.0;
		GradChiRef_rst[0][i*Nbf+4] = sqrt(1.0/2.0)*0.0*sqrt(1.0/2.0)*1.0*sqrt(3.0/2.0)*t_i;
		GradChiRef_rst[0][i*Nbf+5] = sqrt(3.0/2.0)*1.0*sqrt(1.0/2.0)*1.0*sqrt(3.0/2.0)*t_i;
		GradChiRef_rst[0][i*Nbf+6] = sqrt(1.0/2.0)*0.0*sqrt(3.0/2.0)*s_i*sqrt(3.0/2.0)*t_i;
		GradChiRef_rst[0][i*Nbf+7] = sqrt(3.0/2.0)*1.0*sqrt(3.0/2.0)*s_i*sqrt(3.0/2.0)*t_i;

		GradChiRef_rst[1][i*Nbf+0] = sqrt(1.0/2.0)*1.0*sqrt(1.0/2.0)*0.0*sqrt(1.0/2.0)*1.0;
		GradChiRef_rst[1][i*Nbf+1] = sqrt(3.0/2.0)*r_i*sqrt(1.0/2.0)*0.0*sqrt(1.0/2.0)*1.0;
		GradChiRef_rst[1][i*Nbf+2] = sqrt(1.0/2.0)*1.0*sqrt(3.0/2.0)*1.0*sqrt(1.0/2.0)*1.0;
		GradChiRef_rst[1][i*Nbf+3] = sqrt(3.0/2.0)*r_i*sqrt(3.0/2.0)*1.0*sqrt(1.0/2.0)*1.0;
		GradChiRef_rst[1][i*Nbf+4] = sqrt(1.0/2.0)*1.0*sqrt(1.0/2.0)*0.0*sqrt(3.0/2.0)*t_i;
		GradChiRef_rst[1][i*Nbf+5] = sqrt(3.0/2.0)*r_i*sqrt(1.0/2.0)*0.0*sqrt(3.0/2.0)*t_i;
		GradChiRef_rst[1][i*Nbf+6] = sqrt(1.0/2.0)*1.0*sqrt(3.0/2.0)*1.0*sqrt(3.0/2.0)*t_i;
		GradChiRef_rst[1][i*Nbf+7] = sqrt(3.0/2.0)*r_i*sqrt(3.0/2.0)*1.0*sqrt(3.0/2.0)*t_i;

		GradChiRef_rst[2][i*Nbf+0] = sqrt(1.0/2.0)*1.0*sqrt(1.0/2.0)*1.0*sqrt(1.0/2.0)*0.0;
		GradChiRef_rst[2][i*Nbf+1] = sqrt(3.0/2.0)*r_i*sqrt(1.0/2.0)*1.0*sqrt(1.0/2.0)*0.0;
		GradChiRef_rst[2][i*Nbf+2] = sqrt(1.0/2.0)*1.0*sqrt(3.0/2.0)*s_i*sqrt(1.0/2.0)*0.0;
		GradChiRef_rst[2][i*Nbf+3] = sqrt(3.0/2.0)*r_i*sqrt(3.0/2.0)*s_i*sqrt(1.0/2.0)*0.0;
		GradChiRef_rst[2][i*Nbf+4] = sqrt(1.0/2.0)*1.0*sqrt(1.0/2.0)*1.0*sqrt(3.0/2.0)*1.0;
		GradChiRef_rst[2][i*Nbf+5] = sqrt(3.0/2.0)*r_i*sqrt(1.0/2.0)*1.0*sqrt(3.0/2.0)*1.0;
		GradChiRef_rst[2][i*Nbf+6] = sqrt(1.0/2.0)*1.0*sqrt(3.0/2.0)*s_i*sqrt(3.0/2.0)*1.0;
		GradChiRef_rst[2][i*Nbf+7] = sqrt(3.0/2.0)*r_i*sqrt(3.0/2.0)*s_i*sqrt(3.0/2.0)*1.0;
	}

	free(r);
	free(s);
	free(t);

	return GradChiRef_rst;
}

static void poly2(const double *r, const double *s, const double *t, const unsigned int Nn, double **f, double **f_r,
                  double **f_s, double **f_t)
{
	unsigned int i;
	double *f_rst, *f_r_rst, *f_s_rst, *f_t_rst;
	double r_i, s_i, t_i;

	f_rst   = malloc(Nn * sizeof *f_rst); // keep (requires external free)
	f_r_rst = malloc(Nn * sizeof *f_r_rst); // keep (requires external free)
	f_s_rst = malloc(Nn * sizeof *f_s_rst); // keep (requires external free)
	f_t_rst = malloc(Nn * sizeof *f_t_rst); // keep (requires external free)

	for (i = 0; i < Nn; i++) {
		r_i = r[i];
		s_i = s[i];
		t_i = t[i];

		f_rst[i] = 1.0
		         + 2.0*r_i + 3.0*s_i + 4.0*t_i
		         + 5.0*r_i*r_i + 6.0*r_i*s_i + 7.0*r_i*t_i + 8.0*s_i*s_i + 9.0*s_i*t_i + 10.0*t_i*t_i;
		f_r_rst[i] = 2.0
		           + 10.0*r_i + 6.0*s_i + 7.0*t_i;
		f_s_rst[i] = 3.0
		           + 6.0*r_i + 16.0*s_i + 9.0*t_i;
		f_t_rst[i] = 4.0
		           + 7.0*r_i + 9.0*s_i + 20.0*t_i;
	}

	*f   = f_rst;
	*f_r = f_r_rst;
	*f_s = f_s_rst;
	*f_t = f_t_rst;
}

void test_imp_grad_basis_TP(void)
{
	unsigned int pass;

	/*
	 *	grad_basis_TP (dE = 1):
	 *
	 *		Input:
	 *
	 *			P, rst, Nn, dE.
	 *
	 *		Expected output:
	 *
	 *			P = 3:
	 *				GradChiRef_rst = @(r) [0 sqrt(3/2) sqrt(5/2)*3*r sqrt(7/2)*1/2*(15*r^2-3)]
	 */

	unsigned int dE, Nn, Nbf, P, Prst;
	double *rst, *w;

	dE = 1;

	double **GradChiRef13_code, **GradChiRef13_test;

	P = 3;
	Prst = 4;

	cubature_TP(&rst,&w,&Nn,0,Prst,dE,"GL"); // free

	GradChiRef13_code = grad_basis_TP(P,rst,Nn,&Nbf,dE); // free
	GradChiRef13_test = grad_basis_TP13(rst,Nn); // free

	pass = 0;
	if (array_norm_diff_d(Nn*pow(P+1,dE),GradChiRef13_code[0],GradChiRef13_test[0],"Inf") < EPS*10)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("grad_basis_TP (d1, P3):                          ");
	test_print(pass);

	free(rst);
	array_free2_d(dE,GradChiRef13_code);
	array_free2_d(dE,GradChiRef13_test);

	/*
	 *	grad_basis_TP (dE = 2):
	 *
	 *		Input:
	 *
	 *			P, rst, Nn, dE.
	 *
	 *		Expected output:
	 *
	 *			P = 2:
	 *				GradChiRef_rst[0] = @(r,s) [0
	 *				                            sqrt(3/2)*sqrt(1/2)*1
	 *				                            sqrt(5/2)*3*r*sqrt(1/2)*1
	 *				                            0
	 *				                            sqrt(3/2)*sqrt(3/2)*s
	 *				                            sqrt(5/2)*3*r*sqrt(3/2)*s
	 *				                            0
	 *				                            sqrt(3/2)*sqrt(5/2)*1/2*(3*s^2-1)
	 *				                            sqrt(5/2)*3*r*sqrt(5/2)*1/2*(3*s^2-1)]
	 *
	 *				GradChiRef_rst[1] = @(r,s) [0
	 *				                            0
	 *				                            0
	 *				                            sqrt(1/2)*1*sqrt(3/2)
	 *				                            sqrt(3/2)*r*sqrt(3/2)
	 *				                            sqrt(5/2)*1/2*(3*r^2-1)*sqrt(3/2)
	 *				                            sqrt(1/2)*1*sqrt(5/2)*3*s
	 *				                            sqrt(3/2)*r*sqrt(5/2)*3*s
	 *				                            sqrt(5/2)*1/2*(3*r^2-1)*sqrt(5/2)*3*s]
	 */

	dE = 2;

	double **GradChiRef22_code, **GradChiRef22_test;

	P = 2;
	Prst = 4;

	cubature_TP(&rst,&w,&Nn,0,Prst,dE,"GL"); // free

	GradChiRef22_code = grad_basis_TP(P,rst,Nn,&Nbf,dE); // free
	GradChiRef22_test = grad_basis_TP22(rst,Nn); // free

	pass = 0;
	if (array_norm_diff_d(Nn*pow(P+1,dE),GradChiRef22_code[0],GradChiRef22_test[0],"Inf") < EPS*10 &&
	    array_norm_diff_d(Nn*pow(P+1,dE),GradChiRef22_code[1],GradChiRef22_test[1],"Inf") < EPS*10)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("              (d2, P2):                          ");
	test_print(pass);

	free(rst);
	array_free2_d(dE,GradChiRef22_code);
	array_free2_d(dE,GradChiRef22_test);

	/*
	 *	grad_basis_TP (dE = 3):
	 *
	 *		Input:
	 *
	 *			P, rst, Nn, dE.
	 *
	 *		Expected output:
	 *
	 *			P = 1:
	 *				GradChiRef_rst[0] = @(r,s,t) [0
	 *				                              sqrt(3/2)*sqrt(1/2)*1*sqrt(1/2)*1
	 *				                              0
	 *				                              sqrt(3/2)*sqrt(3/2)*1*sqrt(1/2)*1
	 *				                              0
	 *				                              sqrt(3/2)*sqrt(1/2)*1*sqrt(3/2)*t
	 *				                              0
	 *				                              sqrt(3/2)*sqrt(3/2)*s*sqrt(3/2)*t]
	 *
	 *				GradChiRef_rst[1] = @(r,s,t) [0
	 *				                              0
	 *				                              sqrt(1/2)*1*sqrt(3/2)*sqrt(1/2)*1
	 *				                              sqrt(3/2)*r*sqrt(3/2)*sqrt(1/2)*1
	 *				                              0
	 *				                              0
	 *				                              sqrt(1/2)*1*sqrt(3/2)*sqrt(3/2)*t
	 *				                              sqrt(3/2)*r*sqrt(3/2)*sqrt(3/2)*t
	 *
	 *				GradChiRef_rst[2] = @(r,s,t) [0
	 *				                              0
	 *				                              0
	 *				                              0
	 *				                              sqrt(1/2)*1*sqrt(1/2)*1*sqrt(3/2)
	 *				                              sqrt(3/2)*r*sqrt(1/2)*1*sqrt(3/2)
	 *				                              sqrt(1/2)*1*sqrt(3/2)*s*sqrt(3/2)
	 *				                              sqrt(3/2)*r*sqrt(3/2)*s*sqrt(3/2)]
	 */

	dE = 3;

	double **GradChiRef31_code, **GradChiRef31_test;

	P = 1;
	Prst = 4;

	cubature_TP(&rst,&w,&Nn,0,Prst,dE,"GL"); // free

	GradChiRef31_code = grad_basis_TP(P,rst,Nn,&Nbf,dE); // free
	GradChiRef31_test = grad_basis_TP31(rst,Nn); // free

	pass = 0;
	if (array_norm_diff_d(Nn*pow(P+1,dE),GradChiRef31_code[0],GradChiRef31_test[0],"Inf") < EPS*10 &&
	    array_norm_diff_d(Nn*pow(P+1,dE),GradChiRef31_code[1],GradChiRef31_test[1],"Inf") < EPS*10 &&
	    array_norm_diff_d(Nn*pow(P+1,dE),GradChiRef31_code[2],GradChiRef31_test[2],"Inf") < EPS*10)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("              (d3, P1):                          ");
	test_print(pass);

	free(rst);
	array_free2_d(dE,GradChiRef31_code);
	array_free2_d(dE,GradChiRef31_test);

	/*
	 *	grad_basis_TP, derivatives:
	 *
	 *		Input:
	 *
	 *			f(r,s,t), f_r(r,s,t), f_s(r,s,t), f_t(r,s,t), rst
	 *
	 *		Expected Output:
	 *
	 *			GradChiRef_rst[0]*ChiRefInvP_vP*f(rst_vP) = f_r(rst_vP)
	 *			GradChiRef_rst[1]*ChiRefInvP_vP*f(rst_vP) = f_s(rst_vP)
	 *			GradChiRef_rst[2]*ChiRefInvP_vP*f(rst_vP) = f_t(rst_vP)
	 */

	unsigned int i;
	double *f, *f_r, *f_s, *f_t, *r, *s, *t, *f_hat, *f_rcomp, *f_scomp, *f_tcomp;
	double *I, *ChiRef_rst, *ChiRefInv_rst, **GradChiRef_rst;

	P = 4; // should work for P >= 2

	// dE = 1
	dE = 1;

	cubature_TP(&rst,&w,&Nn,0,P,dE,"GL"); // free

	r = malloc(Nn * sizeof *r); // free
	s = malloc(Nn * sizeof *s); // free
	t = malloc(Nn * sizeof *t); // free
	for (i = 0; i < Nn; i++) {
		r[i] = rst[0*Nn+i];
		s[i] = 0.0;
		t[i] = 0.0;
	}

	I = identity_d(Nn); // free
	ChiRef_rst = basis_TP(P,rst,Nn,&Nbf,dE); // free
	ChiRefInv_rst = inverse_d(Nn,Nn,ChiRef_rst,I); // free

	GradChiRef_rst = grad_basis_TP(P,rst,Nn,&Nbf,dE); // free

	poly2(r,s,t,Nn,&f,&f_r,&f_s,&f_t); // free

	f_hat = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,1,Nn,1.0,ChiRefInv_rst,f); // free
	f_rcomp = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,1,Nn,1.0,GradChiRef_rst[0],f_hat); // free
	f_scomp = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,1,Nn,1.0,GradChiRef_rst[0],f_hat); // free
	f_tcomp = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,1,Nn,1.0,GradChiRef_rst[0],f_hat); // free

	pass = 0;
	if (array_norm_diff_d(Nn,f_r,f_rcomp,"Inf") < EPS*100)
		pass = 1, TestDB.Npass++;
	else
		printf("%e\n",array_norm_diff_d(Nn,f_r,f_rcomp,"Inf"));


	//     0         10        20        30        40        50
	printf("              derivative (d1):                   ");
	test_print(pass);

	free(rst), free(r), free(s), free(t);
	free(I);
	free(ChiRef_rst), free(ChiRefInv_rst);
	array_free2_d(dE,GradChiRef_rst);
	free(f), free(f_hat);
	free(f_r), free(f_s), free(f_t);
	free(f_rcomp), free(f_scomp), free(f_tcomp);

	// dE = 2
	dE = 2;

	cubature_TP(&rst,&w,&Nn,0,P,dE,"GL"); // free

	r = malloc(Nn * sizeof *r); // free
	s = malloc(Nn * sizeof *s); // free
	t = malloc(Nn * sizeof *t); // free
	for (i = 0; i < Nn; i++) {
		r[i] = rst[0*Nn+i];
		s[i] = rst[1*Nn+i];
		t[i] = 0.0;
	}

	I = identity_d(Nn); // free
	ChiRef_rst = basis_TP(P,rst,Nn,&Nbf,dE); // free
	ChiRefInv_rst = inverse_d(Nn,Nn,ChiRef_rst,I); // free

	GradChiRef_rst = grad_basis_TP(P,rst,Nn,&Nbf,dE); // free

	poly2(r,s,t,Nn,&f,&f_r,&f_s,&f_t); // free

	f_hat = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,1,Nn,1.0,ChiRefInv_rst,f); // free
	f_rcomp = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,1,Nn,1.0,GradChiRef_rst[0],f_hat); // free
	f_scomp = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,1,Nn,1.0,GradChiRef_rst[1],f_hat); // free
	f_tcomp = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,1,Nn,1.0,GradChiRef_rst[0],f_hat); // free

	pass = 0;
	if (array_norm_diff_d(Nn,f_r,f_rcomp,"Inf") < EPS*100 &&
	    array_norm_diff_d(Nn,f_s,f_scomp,"Inf") < EPS*100) {
			pass = 1, TestDB.Npass++;
	} else {
		printf("%e\n",array_norm_diff_d(Nn,f_r,f_rcomp,"Inf"));
		printf("%e\n",array_norm_diff_d(Nn,f_s,f_scomp,"Inf"));
	}

	//     0         10        20        30        40        50
	printf("              derivatives (d2):                  ");
	test_print(pass);

	free(rst), free(r), free(s), free(t);
	free(I);
	free(ChiRef_rst), free(ChiRefInv_rst);
	array_free2_d(dE,GradChiRef_rst);
	free(f), free(f_hat);
	free(f_r), free(f_s), free(f_t);
	free(f_rcomp), free(f_scomp), free(f_tcomp);

	// dE = 3
	dE = 3;

	cubature_TP(&rst,&w,&Nn,0,P,dE,"GL"); // free

	r = malloc(Nn * sizeof *r); // free
	s = malloc(Nn * sizeof *s); // free
	t = malloc(Nn * sizeof *t); // free
	for (i = 0; i < Nn; i++) {
		r[i] = rst[0*Nn+i];
		s[i] = rst[1*Nn+i];
		t[i] = rst[2*Nn+i];
	}

	I = identity_d(Nn); // free
	ChiRef_rst = basis_TP(P,rst,Nn,&Nbf,dE); // free
	ChiRefInv_rst = inverse_d(Nn,Nn,ChiRef_rst,I); // free

	GradChiRef_rst = grad_basis_TP(P,rst,Nn,&Nbf,dE); // free

	poly2(r,s,t,Nn,&f,&f_r,&f_s,&f_t); // free

	f_hat = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,1,Nn,1.0,ChiRefInv_rst,f); // free
	f_rcomp = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,1,Nn,1.0,GradChiRef_rst[0],f_hat); // free
	f_scomp = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,1,Nn,1.0,GradChiRef_rst[1],f_hat); // free
	f_tcomp = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,1,Nn,1.0,GradChiRef_rst[2],f_hat); // free

	pass = 0;
	if (array_norm_diff_d(Nn,f_r,f_rcomp,"Inf") < EPS*1000 &&
	    array_norm_diff_d(Nn,f_s,f_scomp,"Inf") < EPS*1000 &&
	    array_norm_diff_d(Nn,f_t,f_tcomp,"Inf") < EPS*1000) {
			pass = 1, TestDB.Npass++;
	} else {
		printf("%e\n",array_norm_diff_d(Nn,f_r,f_rcomp,"Inf"));
		printf("%e\n",array_norm_diff_d(Nn,f_s,f_scomp,"Inf"));
		printf("%e\n",array_norm_diff_d(Nn,f_t,f_tcomp,"Inf"));
	}

	//     0         10        20        30        40        50
	printf("              derivatives (d3):                  ");
	test_print(pass);

	free(rst), free(r), free(s), free(t);
	free(I);
	free(ChiRef_rst), free(ChiRefInv_rst);
	array_free2_d(dE,GradChiRef_rst);
	free(f), free(f_hat);
	free(f_r), free(f_s), free(f_t);
	free(f_rcomp), free(f_scomp), free(f_tcomp);
}
