// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "test_unit_grad_bases.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "mkl.h"

#include "Parameters.h"
#include "Test.h"

#include "test_support.h"
#include "test_code_bases.h"
#include "array_norm.h"
#include "cubature.h"
#include "bases.h"
#include "array_free.h"
#include "matrix_functions.h"

#undef I // No complex variables used here

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

void test_unit_grad_basis_TP(void)
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
	 *				GradChiRef_rst = @(r) [ See grad_basis_TP13 ]
	 */

	unsigned int dE, Nn, Nbf, P, Prst;
	double *rst;

	struct S_CUBATURE *CUBDATA = malloc(sizeof *CUBDATA); // free

	dE = 1;

	double **GradChiRef13_code, **GradChiRef13_test;

	P = 3;
	Prst = 4;

	set_cubdata(CUBDATA,false,false,"GL",dE,Prst,cubature_TP); // free
	set_from_cubdata(CUBDATA,&Nn,NULL,&rst,NULL,NULL);

	GradChiRef13_code = grad_basis_TP(P,rst,Nn,&Nbf,dE); // free
	GradChiRef13_test = grad_basis_TP13(rst,Nn); // free

	pass = 0;
	if (array_norm_diff_d(Nn*pow(P+1,dE),GradChiRef13_code[0],GradChiRef13_test[0],"Inf") < EPS*10)
		pass = 1;

	test_print2(pass,"grad_basis_TP (d1, P3):");

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
	 *				GradChiRef_rst[0] = @(r,s) [ See grad_basis_TP22[0] ]
	 *				GradChiRef_rst[1] = @(r,s) [ See grad_basis_TP22[1] ]
	 */

	dE = 2;

	double **GradChiRef22_code, **GradChiRef22_test;

	P = 2;
	Prst = 4;

	set_cubdata(CUBDATA,false,false,"GL",dE,Prst,cubature_TP); // free
	set_from_cubdata(CUBDATA,&Nn,NULL,&rst,NULL,NULL);

	GradChiRef22_code = grad_basis_TP(P,rst,Nn,&Nbf,dE); // free
	GradChiRef22_test = grad_basis_TP22(rst,Nn); // free

	pass = 0;
	if (array_norm_diff_d(Nn*pow(P+1,dE),GradChiRef22_code[0],GradChiRef22_test[0],"Inf") < EPS*10 &&
	    array_norm_diff_d(Nn*pow(P+1,dE),GradChiRef22_code[1],GradChiRef22_test[1],"Inf") < EPS*10)
			pass = 1;

	test_print2(pass,"              (d2, P2):");

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
	 *				GradChiRef_rst[0] = @(r,s,t) [ See grad_basis_TP31[0] ]
	 *				GradChiRef_rst[1] = @(r,s,t) [ See grad_basis_TP31[1] ]
	 *				GradChiRef_rst[2] = @(r,s,t) [ See grad_basis_TP31[2] ]
	 */

	dE = 3;

	double **GradChiRef31_code, **GradChiRef31_test;

	P = 1;
	Prst = 4;

	set_cubdata(CUBDATA,false,false,"GL",dE,Prst,cubature_TP); // free
	set_from_cubdata(CUBDATA,&Nn,NULL,&rst,NULL,NULL);

	GradChiRef31_code = grad_basis_TP(P,rst,Nn,&Nbf,dE); // free
	GradChiRef31_test = grad_basis_TP31(rst,Nn); // free

	pass = 0;
	if (array_norm_diff_d(Nn*pow(P+1,dE),GradChiRef31_code[0],GradChiRef31_test[0],"Inf") < EPS*10 &&
	    array_norm_diff_d(Nn*pow(P+1,dE),GradChiRef31_code[1],GradChiRef31_test[1],"Inf") < EPS*10 &&
	    array_norm_diff_d(Nn*pow(P+1,dE),GradChiRef31_code[2],GradChiRef31_test[2],"Inf") < EPS*10)
			pass = 1;

	test_print2(pass,"              (d3, P1):");

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

	set_cubdata(CUBDATA,false,false,"GL",dE,P,cubature_TP); // free
	set_from_cubdata(CUBDATA,&Nn,NULL,&rst,NULL,NULL);

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
		pass = 1;
	else
		printf("%e\n",array_norm_diff_d(Nn,f_r,f_rcomp,"Inf"));


	test_print2(pass,"              derivative  (d1):");

	free(rst), free(r), free(s), free(t);
	free(I);
	free(ChiRef_rst), free(ChiRefInv_rst);
	array_free2_d(dE,GradChiRef_rst);
	free(f), free(f_hat);
	free(f_r), free(f_s), free(f_t);
	free(f_rcomp), free(f_scomp), free(f_tcomp);

	// dE = 2
	dE = 2;

	set_cubdata(CUBDATA,false,false,"GL",dE,P,cubature_TP); // free
	set_from_cubdata(CUBDATA,&Nn,NULL,&rst,NULL,NULL);

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
			pass = 1;
	} else {
		printf("%e\n",array_norm_diff_d(Nn,f_r,f_rcomp,"Inf"));
		printf("%e\n",array_norm_diff_d(Nn,f_s,f_scomp,"Inf"));
	}

	test_print2(pass,"              derivatives (d2):");

	free(rst), free(r), free(s), free(t);
	free(I);
	free(ChiRef_rst), free(ChiRefInv_rst);
	array_free2_d(dE,GradChiRef_rst);
	free(f), free(f_hat);
	free(f_r), free(f_s), free(f_t);
	free(f_rcomp), free(f_scomp), free(f_tcomp);

	// dE = 3
	dE = 3;

	set_cubdata(CUBDATA,false,false,"GL",dE,P,cubature_TP); // free
	set_from_cubdata(CUBDATA,&Nn,NULL,&rst,NULL,NULL);

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
			pass = 1;
	} else {
		printf("%e\n",array_norm_diff_d(Nn,f_r,f_rcomp,"Inf"));
		printf("%e\n",array_norm_diff_d(Nn,f_s,f_scomp,"Inf"));
		printf("%e\n",array_norm_diff_d(Nn,f_t,f_tcomp,"Inf"));
	}

	test_print2(pass,"              derivatives (d3):");

	free(rst), free(r), free(s), free(t);
	free(I);
	free(ChiRef_rst), free(ChiRefInv_rst);
	array_free2_d(dE,GradChiRef_rst);
	free(f), free(f_hat);
	free(f_r), free(f_s), free(f_t);
	free(f_rcomp), free(f_scomp), free(f_tcomp);

	free(CUBDATA);
}

/*
 *	Purpose:
 *		Test correctness of implementation of grad_basis_SI.
 *
 *	Comments:
 *		May be better to test grad_basis_SI of order P > 1 to be thorough. (ToBeDeleted)
 *
 *	Notation:
 *
 *	References:
 */

static double **grad_basis_SI22(const double *rst, const unsigned int Nn)
{
	unsigned int n, Nbf, dim, i, j;
	double **GradChiRef_rst, *a, *b, *c, a_n, b_n, con, con_i, con_j, con_b;

	Nbf = 6;

	GradChiRef_rst = malloc(2 * sizeof *GradChiRef_rst); // keep (requires external free)
	for (dim = 0; dim < 2; dim++)
		GradChiRef_rst[dim] = malloc(Nn*Nbf * sizeof **GradChiRef_rst); // keep (requires external free)

	// Convert from rst to abc coordinates
	a = malloc(Nn * sizeof *a); // free
	b = malloc(Nn * sizeof *b); // free
	c = malloc(Nn * sizeof *c); // free

	rst_to_abc_SI(Nn,2,rst,a,b,c);

	con = 2.0/pow(3.0,0.25);
	for (n = 0; n < Nn; n++) {
		a_n = a[n];
		b_n = b[n];

		i = 0; j = 0;
		get_scaling_basis_TRI(i,j,b_n,&con_i,&con_j,&con_b);
		GradChiRef_rst[0][n*Nbf+0] = 0.0;
		GradChiRef_rst[1][n*Nbf+0] = 0.0;

		i = 0; j = 1;
		get_scaling_basis_TRI(i,j,b_n,&con_i,&con_j,&con_b);
		GradChiRef_rst[0][n*Nbf+1] = 0.0;
		GradChiRef_rst[1][n*Nbf+1] = con*( 2.0/3.0*sqrt(3.0)*pow(1.0-b_n,(double) i)
		                                 * con_i*(1.0)
		                                 * con_j*(3.0/2.0) );

		i = 0; j = 2;
		get_scaling_basis_TRI(i,j,b_n,&con_i,&con_j,&con_b);
		GradChiRef_rst[0][n*Nbf+2] = 0.0;
		GradChiRef_rst[1][n*Nbf+2] = con*( 2.0/3.0*sqrt(3.0)*pow(1.0-b_n,(double) i)
		                                 * con_i*(1.0)
		                                 * con_j*(5.0*(b_n-1.0)+6.0) );

		i = 1; j = 0;
		get_scaling_basis_TRI(i,j,b_n,&con_i,&con_j,&con_b);
		GradChiRef_rst[0][n*Nbf+3] = con*2.0*pow(1.0-b_n,i-1.0)
		                           * con_i*(1.0)
		                           * con_j*(1.0);
		GradChiRef_rst[1][n*Nbf+3] = 0.0;

		i = 1; j = 1;
		get_scaling_basis_TRI(i,j,b_n,&con_i,&con_j,&con_b);
		GradChiRef_rst[0][n*Nbf+4] = con*2.0*pow(1.0-b_n,i-1.0)
		                           * con_i*(1.0)
		                           * con_j*(1.0/2.0*(5.0*b_n+3.0));
		GradChiRef_rst[1][n*Nbf+4] = con*( -2.0/3.0*sqrt(3.0)*i*pow(1.0-b_n,i-1.0)
		                                 * con_i*(a_n)
		                                 * con_j*(1.0/2.0*(5.0*b_n+3.0))
		                                 + 2.0/3.0*sqrt(3.0)*a_n*pow(1.0-b_n,i-1.0)
		                                 * con_i*(1.0)
		                                 * con_j*(1.0/2.0*(5.0*b_n+3.0))
		                                 + 2.0/3.0*sqrt(3.0)*pow(1.0-b_n,(double) i)
		                                 * con_i*(a_n)
		                                 * con_j*(5.0/2.0) );

		i = 2; j = 0;
		get_scaling_basis_TRI(i,j,b_n,&con_i,&con_j,&con_b);
		GradChiRef_rst[0][n*Nbf+5] = con*2.0*pow(1.0-b_n,i-1.0)
		                           * con_i*(3.0*a_n)
		                           * con_j*(1.0);
		GradChiRef_rst[1][n*Nbf+5] = con*( -2.0/3.0*sqrt(3.0)*i*pow(1.0-b_n,i-1.0)
		                                 * con_i*(1.0/2.0*(3.0*pow(a_n,2.0)-1.0))
		                                 * con_j*(1.0)
		                                 + 2.0/3.0*sqrt(3.0)*a_n*pow(1.0-b_n,i-1.0)
		                                 * con_i*(3.0*a_n)
		                                 * con_j*(1.0) );
	}

	free(a);
	free(b);
	free(c);

	return GradChiRef_rst;
}

static double **grad_basis_SI31(const double *rst, const unsigned int Nn)
{
	unsigned int n, Nbf, dim, i, j, k;
	double **GradChiRef_rst, *a, *b, *c, a_n, b_n, c_n, con, con_i, con_j, con_k, con_b, con_c;

	// silence
	a_n = 0, b_n = a_n+1;

	Nbf = 4;

	GradChiRef_rst = malloc(3 * sizeof *GradChiRef_rst); // keep (requires external free)
	for (dim = 0; dim < 3; dim++)
		GradChiRef_rst[dim] = malloc(Nn*Nbf * sizeof **GradChiRef_rst); // keep (requires external free)

	// Convert from rst to abc coordinates
	a = malloc(Nn * sizeof *a); // free
	b = malloc(Nn * sizeof *b); // free
	c = malloc(Nn * sizeof *c); // free

	rst_to_abc_SI(Nn,3,rst,a,b,c);

	con = 4.0/pow(2.0,0.25);
	for (n = 0; n < Nn; n++) {
		a_n = a[n];
		b_n = b[n];
		c_n = c[n];

		i = 0; j = 0; k = 0;
		get_scaling_basis_TET(i,j,k,b_n,c_n,&con_i,&con_j,&con_k,&con_b,&con_c);
		GradChiRef_rst[0][n*Nbf+0] = 0.0;
		GradChiRef_rst[1][n*Nbf+0] = 0.0;
		GradChiRef_rst[2][n*Nbf+0] = 0.0;

		i = 0; j = 0; k = 1;
		get_scaling_basis_TET(i,j,k,b_n,c_n,&con_i,&con_j,&con_k,&con_b,&con_c);
		GradChiRef_rst[0][n*Nbf+1] = 0.0;
		GradChiRef_rst[1][n*Nbf+1] = 0.0;
		GradChiRef_rst[2][n*Nbf+1] = con*( 1.0/2.0*sqrt(6.0)*pow(1.0-b_n,(double) i)*pow(1.0-c_n,(double) (i+j))
		                                 * con_i*(1.0)
		                                 * con_j*(1.0)
		                                 * con_k*(2.0) );

		i = 0; j = 1; k = 0;
		get_scaling_basis_TET(i,j,k,b_n,c_n,&con_i,&con_j,&con_k,&con_b,&con_c);
		GradChiRef_rst[0][n*Nbf+2] = 0.0;
		GradChiRef_rst[1][n*Nbf+2] = con*( 4.0/3.0*sqrt(3.0)*pow(1.0-b_n,(double) i)*pow(1.0-c_n,i+j-1.0)
		                                 * con_i*(1.0)
		                                 * con_j*(3.0/2.0)
		                                 * con_k*(1.0) );
		GradChiRef_rst[2][n*Nbf+2] = 0.0;

		i = 1; j = 0; k = 0;
		get_scaling_basis_TET(i,j,k,b_n,c_n,&con_i,&con_j,&con_k,&con_b,&con_c);
		GradChiRef_rst[0][n*Nbf+3] = con*( 4.0*pow(1.0-b_n,i-1.0)*pow(1.0-c_n,i+j-1.0)
		                                 * con_i*(1.0)
		                                 * con_j*(1.0)
		                                 * con_k*(1.0) );
		GradChiRef_rst[1][n*Nbf+3] = 0.0;
		GradChiRef_rst[2][n*Nbf+3] = 0.0;
	}

	free(a);
	free(b);
	free(c);

	return GradChiRef_rst;
}

void test_unit_grad_basis_SI(void)
{
	unsigned int pass;

	/*
	 *	grad_basis_SI (d = 2):
	 *
	 *		Input:
	 *
	 *			P, rst, Nn
	 *
	 *		Expected output:
	 *
	 *			P = 2:
	 *				GradChiRef_rst[0] = @(r,s) [ See grad_basis_SI22[0] ]
	 *				GradChiRef_rst[1] = @(r,s) [ See grad_basis_SI22[1] ]
	 */

	unsigned int d, Nn, Nbf, P, Prst;
	double *rst;

	struct S_CUBATURE *CUBDATA = malloc(sizeof *CUBDATA); // free

	d = 2;

	double **GradChiRef22_code, **GradChiRef22_test;

	P = 2;
	Prst = 4;

	set_cubdata(CUBDATA,false,false,"AO",d,Prst,cubature_TRI); // free
	set_from_cubdata(CUBDATA,&Nn,NULL,&rst,NULL,NULL);

	GradChiRef22_code = grad_basis_SI(P,rst,Nn,&Nbf,d); // free
	GradChiRef22_test = grad_basis_SI22(rst,Nn);        // free

	pass = 0;
	if (array_norm_diff_d(Nn*Nbf,GradChiRef22_code[0],GradChiRef22_test[0],"Inf") < EPS*10 &&
	    array_norm_diff_d(Nn*Nbf,GradChiRef22_code[1],GradChiRef22_test[1],"Inf") < EPS*10)
			pass = 1;

	test_print2(pass,"grad_basis_SI (d2, P2):");

	free(rst);
	array_free2_d(d,GradChiRef22_code);
	array_free2_d(d,GradChiRef22_test);

	/*
	 *	grad_basis_SI (d = 3):
	 *
	 *		Input:
	 *
	 *			P, rst, Nn
	 *
	 *		Expected output:
	 *
	 *			P = 1:
	 *				GradChiRef_rst[0] = @(r,s,t) [ See grad_basis_SI31[0] ]
	 *				GradChiRef_rst[1] = @(r,s,t) [ See grad_basis_SI31[1] ]
	 *				GradChiRef_rst[2] = @(r,s,t) [ See grad_basis_SI31[2] ]
	 */

	d = 3;

	double **GradChiRef31_code, **GradChiRef31_test;

	P = 1;
	Prst = 4;

	set_cubdata(CUBDATA,false,false,"AO",d,Prst,cubature_TET); // free
	set_from_cubdata(CUBDATA,&Nn,NULL,&rst,NULL,NULL);

	GradChiRef31_code = grad_basis_SI(P,rst,Nn,&Nbf,d); // free
	GradChiRef31_test = grad_basis_SI31(rst,Nn); // free

	pass = 0;
	if (array_norm_diff_d(Nn*Nbf,GradChiRef31_code[0],GradChiRef31_test[0],"Inf") < EPS*10 &&
	    array_norm_diff_d(Nn*Nbf,GradChiRef31_code[1],GradChiRef31_test[1],"Inf") < EPS*10 &&
	    array_norm_diff_d(Nn*Nbf,GradChiRef31_code[2],GradChiRef31_test[2],"Inf") < EPS*10)
			pass = 1;

	test_print2(pass,"              (d3, P1):");

	free(rst);
	array_free2_d(d,GradChiRef31_code);
	array_free2_d(d,GradChiRef31_test);

	/*
	 *	grad_basis_SI, derivatives:
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

	// d = 2
	d = 2;

	set_cubdata(CUBDATA,false,false,"AO",d,P,cubature_TRI); // free
	set_from_cubdata(CUBDATA,&Nn,NULL,&rst,NULL,NULL);

	r = malloc(Nn * sizeof *r); // free
	s = malloc(Nn * sizeof *s); // free
	t = malloc(Nn * sizeof *t); // free
	for (i = 0; i < Nn; i++) {
		r[i] = rst[0*Nn+i];
		s[i] = rst[1*Nn+i];
		t[i] = 0.0;
	}

	I = identity_d(Nn); // free
	ChiRef_rst = basis_SI(P,rst,Nn,&Nbf,d); // free
	ChiRefInv_rst = inverse_d(Nn,Nn,ChiRef_rst,I); // free

	GradChiRef_rst = grad_basis_SI(P,rst,Nn,&Nbf,d); // free

	poly2(r,s,t,Nn,&f,&f_r,&f_s,&f_t); // free

	f_hat = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,1,Nn,1.0,ChiRefInv_rst,f); // free
	f_rcomp = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,1,Nn,1.0,GradChiRef_rst[0],f_hat); // free
	f_scomp = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,1,Nn,1.0,GradChiRef_rst[1],f_hat); // free
	f_tcomp = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,1,Nn,1.0,GradChiRef_rst[0],f_hat); // free

	pass = 0;
	if (array_norm_diff_d(Nn,f_r,f_rcomp,"Inf") < EPS*100 &&
	    array_norm_diff_d(Nn,f_s,f_scomp,"Inf") < EPS*100) {
			pass = 1;
	} else {
		printf("%e\n",array_norm_diff_d(Nn,f_r,f_rcomp,"Inf"));
		printf("%e\n",array_norm_diff_d(Nn,f_s,f_scomp,"Inf"));
	}

	test_print2(pass,"              derivatives (d2):");

	free(rst), free(r), free(s), free(t);
	free(I);
	free(ChiRef_rst), free(ChiRefInv_rst);
	array_free2_d(d,GradChiRef_rst);
	free(f), free(f_hat);
	free(f_r), free(f_s), free(f_t);
	free(f_rcomp), free(f_scomp), free(f_tcomp);

	// d = 3
	d = 3;

	set_cubdata(CUBDATA,false,false,"AO",d,P,cubature_TET); // free
	set_from_cubdata(CUBDATA,&Nn,NULL,&rst,NULL,NULL);

	r = malloc(Nn * sizeof *r); // free
	s = malloc(Nn * sizeof *s); // free
	t = malloc(Nn * sizeof *t); // free
	for (i = 0; i < Nn; i++) {
		r[i] = rst[0*Nn+i];
		s[i] = rst[1*Nn+i];
		t[i] = rst[2*Nn+i];
	}

	I = identity_d(Nn); // free
	ChiRef_rst = basis_SI(P,rst,Nn,&Nbf,d); // free
	ChiRefInv_rst = inverse_d(Nn,Nn,ChiRef_rst,I); // free

	GradChiRef_rst = grad_basis_SI(P,rst,Nn,&Nbf,d); // free

	poly2(r,s,t,Nn,&f,&f_r,&f_s,&f_t); // free

	f_hat = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,1,Nn,1.0,ChiRefInv_rst,f); // free
	f_rcomp = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,1,Nn,1.0,GradChiRef_rst[0],f_hat); // free
	f_scomp = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,1,Nn,1.0,GradChiRef_rst[1],f_hat); // free
	f_tcomp = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,1,Nn,1.0,GradChiRef_rst[2],f_hat); // free

	pass = 0;
	if (array_norm_diff_d(Nn,f_r,f_rcomp,"Inf") < EPS*1e3 &&
	    array_norm_diff_d(Nn,f_s,f_scomp,"Inf") < EPS*1e3 &&
	    array_norm_diff_d(Nn,f_t,f_tcomp,"Inf") < EPS*1e3) {
			pass = 1;
	} else {
		printf("%e\n",array_norm_diff_d(Nn,f_r,f_rcomp,"Inf"));
		printf("%e\n",array_norm_diff_d(Nn,f_s,f_scomp,"Inf"));
		printf("%e\n",array_norm_diff_d(Nn,f_t,f_tcomp,"Inf"));
	}

	test_print2(pass,"              derivatives (d3):");

	free(rst), free(r), free(s), free(t);
	free(I);
	free(ChiRef_rst), free(ChiRefInv_rst);
	array_free2_d(d,GradChiRef_rst);
	free(f), free(f_hat);
	free(f_r), free(f_s), free(f_t);
	free(f_rcomp), free(f_scomp), free(f_tcomp);

	free(CUBDATA);
}

/*
 *	Purpose:
 *		Test correctness of implementation of grad_basis_PYR.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

static double **grad_basis_PYR32(const double *rst, const unsigned int Nn)
{
	unsigned int n, Nbf, dim;
	double **GradChiRef_rst, *a, *b, *c, a_n, b_n, c_n, con, con_i, con_j, con_k, con_c;

	Nbf = 14;

	GradChiRef_rst = malloc(3 * sizeof *GradChiRef_rst); // keep (requires external free)
	for (dim = 0; dim < 3; dim++)
		GradChiRef_rst[dim] = calloc(Nn*Nbf , sizeof **GradChiRef_rst); // keep (requires external free)
//		GradChiRef_rst[dim] = malloc(Nn*Nbf * sizeof **GradChiRef_rst); // keep (requires external free)

	// Convert from rst to abc coordinates
	a = malloc(Nn * sizeof *a); // free
	b = malloc(Nn * sizeof *b); // free
	c = malloc(Nn * sizeof *c); // free

	rst_to_abc_PYR(Nn,3,rst,a,b,c);

	con = pow(2.0,1.25);
	for (n = 0; n < Nn; n++) {
		a_n = a[n];
		b_n = b[n];
		c_n = c[n];

		get_scaling_basis_PYR(0,0,0,c_n,&con_i,&con_j,&con_k,&con_c);
		GradChiRef_rst[0][n*Nbf+0] = 0.0;
		GradChiRef_rst[1][n*Nbf+0] = 0.0;
		GradChiRef_rst[2][n*Nbf+0] = 0.0;

		get_scaling_basis_PYR(0,0,1,c_n,&con_i,&con_j,&con_k,&con_c);
		GradChiRef_rst[0][n*Nbf+1] = 0.0;
		GradChiRef_rst[1][n*Nbf+1] = 0.0;
		GradChiRef_rst[2][n*Nbf+1] = con*( sqrt(2.0)*pow(1.0-c_n,0.0)
		                                 * con_i*(1.0)
		                                 * con_j*(1.0)
		                                 * con_k*(1.0/2.0*(4.0)));

		get_scaling_basis_PYR(0,0,2,c_n,&con_i,&con_j,&con_k,&con_c);
		GradChiRef_rst[0][n*Nbf+2] = 0.0;
		GradChiRef_rst[1][n*Nbf+2] = 0.0;
		GradChiRef_rst[2][n*Nbf+2] = con*( sqrt(2.0)*pow(1.0-c_n,0.0)
		                                 * con_i*(1.0)
		                                 * con_j*(1.0)
		                                 * con_k*(15.0/4.0*2.0*pow(c_n-1.0,1.0)+10.0));

		get_scaling_basis_PYR(0,1,0,c_n,&con_i,&con_j,&con_k,&con_c);
		GradChiRef_rst[0][n*Nbf+3] = 0.0;
		GradChiRef_rst[1][n*Nbf+3] = con*( 2.0*pow(1.0-c_n,0.0)
		                                 * con_i*(1.0)
		                                 * con_j*(1.0)
		                                 * con_k*(1.0));
		GradChiRef_rst[2][n*Nbf+3] = con*( sqrt(2.0)*b_n*pow(1.0-c_n,0.0)
		                                 * con_i*(1.0)
		                                 * con_j*(1.0)
		                                 * con_k*(1.0)
		                                 - 1.0*sqrt(2.0)*pow(1.0-c_n,0.0)
		                                 * con_i*(1.0)
		                                 * con_j*(b_n)
		                                 * con_k*(1.0));

		get_scaling_basis_PYR(0,1,1,c_n,&con_i,&con_j,&con_k,&con_c);
		GradChiRef_rst[0][n*Nbf+4] = 0.0;
		GradChiRef_rst[1][n*Nbf+4] = con*( 2.0*pow(1.0-c_n,0.0)
		                                 * con_i*(1.0)
		                                 * con_j*(1.0)
		                                 * con_k*(1.0/2.0*(6.0*c_n+4.0)));
		GradChiRef_rst[2][n*Nbf+4] = con*( sqrt(2.0)*b_n*pow(1.0-c_n,0.0)
		                                 * con_i*(1.0)
		                                 * con_j*(1.0)
		                                 * con_k*(1.0/2.0*(6.0*c_n+4.0))
		                                 + sqrt(2.0)*pow(1.0-c_n,1.0)
		                                 * con_i*(1.0)
		                                 * con_j*(b_n)
		                                 * con_k*(1.0/2.0*(6.0))
		                                 - 1.0*sqrt(2.0)*pow(1.0-c_n,0.0)
		                                 * con_i*(1.0)
		                                 * con_j*(b_n)
		                                 * con_k*(1.0/2.0*(6.0*c_n+4.0)));

		get_scaling_basis_PYR(0,2,0,c_n,&con_i,&con_j,&con_k,&con_c);
		GradChiRef_rst[0][n*Nbf+5] = 0.0;
		GradChiRef_rst[1][n*Nbf+5] = con*( 2.0*pow(1.0-c_n,1.0)
		                                 * con_i*(1.0)
		                                 * con_j*(1.0/2.0*(3.0*2.0*b_n))
		                                 * con_k*(1.0));
		GradChiRef_rst[2][n*Nbf+5] = con*( sqrt(2.0)*b_n*pow(1.0-c_n,1.0)
		                                 * con_i*(1.0)
		                                 * con_j*(1.0/2.0*(3.0*2.0*b_n))
		                                 * con_k*(1.0)
		                                 - 2.0*sqrt(2.0)*pow(1.0-c_n,1.0)
		                                 * con_i*(1.0)
		                                 * con_j*(1.0/2.0*(3.0*pow(b_n,2.0)-1.0))
		                                 * con_k*(1.0));

		get_scaling_basis_PYR(1,0,0,c_n,&con_i,&con_j,&con_k,&con_c);
		GradChiRef_rst[0][n*Nbf+6] = con*( 2.0*pow(1.0-c_n,0.0)
		                                 * con_i*(1.0)
		                                 * con_j*(1.0)
		                                 * con_k*(1.0));
		GradChiRef_rst[1][n*Nbf+6] = 0.0;
		GradChiRef_rst[2][n*Nbf+6] = con*( sqrt(2.0)*a_n*pow(1.0-c_n,0.0)
		                                 * con_i*(1.0)
		                                 * con_j*(1.0)
		                                 * con_k*(1.0)
		                                 - 1.0*sqrt(2.0)*pow(1.0-c_n,0.0)
		                                 * con_i*(a_n)
		                                 * con_j*(1.0)
		                                 * con_k*(1.0));

		get_scaling_basis_PYR(1,0,1,c_n,&con_i,&con_j,&con_k,&con_c);
		GradChiRef_rst[0][n*Nbf+7] = con*( 2.0*pow(1.0-c_n,0.0)
		                                 * con_i*(1.0)
		                                 * con_j*(1.0)
		                                 * con_k*(1.0/2.0*(6.0*c_n+4.0)));
		GradChiRef_rst[1][n*Nbf+7] = 0.0;
		GradChiRef_rst[2][n*Nbf+7] = con*( sqrt(2.0)*a_n*pow(1.0-c_n,0.0)
		                                 * con_i*(1.0)
		                                 * con_j*(1.0)
		                                 * con_k*(1.0/2.0*(6.0*c_n+4.0))
		                                 + sqrt(2.0)*pow(1.0-c_n,1.0)
		                                 * con_i*(a_n)
		                                 * con_j*(1.0)
		                                 * con_k*(1.0/2.0*(6.0))
		                                 - 1.0*sqrt(2.0)*pow(1.0-c_n,0.0)
		                                 * con_i*(a_n)
		                                 * con_j*(1.0)
		                                 * con_k*(1.0/2.0*(6.0*c_n+4.0)));

		get_scaling_basis_PYR(1,1,0,c_n,&con_i,&con_j,&con_k,&con_c);
		GradChiRef_rst[0][n*Nbf+8] = con*( 2.0*pow(1.0-c_n,0.0)
		                                 * con_i*(1.0)
		                                 * con_j*(b_n)
		                                 * con_k*(1.0));
		GradChiRef_rst[1][n*Nbf+8] = con*( 2.0*pow(1.0-c_n,0.0)
		                                 * con_i*(a_n)
		                                 * con_j*(1.0)
		                                 * con_k*(1.0));
		GradChiRef_rst[2][n*Nbf+8] = con*( sqrt(2.0)*a_n*pow(1.0-c_n,0.0)
		                                 * con_i*(1.0)
		                                 * con_j*(b_n)
		                                 * con_k*(1.0)
		                                 + sqrt(2.0)*b_n*pow(1.0-c_n,0.0)
		                                 * con_i*(a_n)
		                                 * con_j*(1.0)
		                                 * con_k*(1.0)
		                                 - 1.0*sqrt(2.0)*pow(1.0-c_n,0.0)
		                                 * con_i*(a_n)
		                                 * con_j*(b_n)
		                                 * con_k*(1.0));

		get_scaling_basis_PYR(1,1,1,c_n,&con_i,&con_j,&con_k,&con_c);
		GradChiRef_rst[0][n*Nbf+9] = con*( 2.0*pow(1.0-c_n,0.0)
		                                 * con_i*(1.0)
		                                 * con_j*(b_n)
		                                 * con_k*(1.0/2.0*(6.0*c_n+4.0)));
		GradChiRef_rst[1][n*Nbf+9] = con*( 2.0*pow(1.0-c_n,0.0)
		                                 * con_i*(a_n)
		                                 * con_j*(1.0)
		                                 * con_k*(1.0/2.0*(6.0*c_n+4.0)));
		GradChiRef_rst[2][n*Nbf+9] = con*( sqrt(2.0)*a_n*pow(1.0-c_n,0.0)
		                                 * con_i*(1.0)
		                                 * con_j*(b_n)
		                                 * con_k*(1.0/2.0*(6.0*c_n+4.0))
		                                 + sqrt(2.0)*b_n*pow(1.0-c_n,0.0)
		                                 * con_i*(a_n)
		                                 * con_j*(1.0)
		                                 * con_k*(1.0/2.0*(6.0*c_n+4.0))
		                                 + sqrt(2.0)*pow(1.0-c_n,1.0)
		                                 * con_i*(a_n)
		                                 * con_j*(b_n)
		                                 * con_k*(1.0/2.0*(6.0))
		                                 - 1.0*sqrt(2.0)*pow(1.0-c_n,0.0)
		                                 * con_i*(a_n)
		                                 * con_j*(b_n)
		                                 * con_k*(1.0/2.0*(6.0*c_n+4.0)));

		get_scaling_basis_PYR(1,2,0,c_n,&con_i,&con_j,&con_k,&con_c);
		GradChiRef_rst[0][n*Nbf+10] = con*( 2.0*pow(1.0-c_n,1.0)
		                                  * con_i*(1.0)
		                                  * con_j*(1.0/2.0*(3.0*pow(b_n,2.0)-1.0))
		                                  * con_k*(1.0));
		GradChiRef_rst[1][n*Nbf+10] = con*( 2.0*pow(1.0-c_n,1.0)
		                                  * con_i*(a_n)
		                                  * con_j*(1.0/2.0*(3.0*2.0*b_n))
		                                  * con_k*(1.0));
		GradChiRef_rst[2][n*Nbf+10] = con*( sqrt(2.0)*a_n*pow(1.0-c_n,1.0)
		                                  * con_i*(1.0)
		                                  * con_j*(1.0/2.0*(3.0*pow(b_n,2.0)-1.0))
		                                  * con_k*(1.0)
		                                  + sqrt(2.0)*b_n*pow(1.0-c_n,1.0)
		                                  * con_i*(a_n)
		                                  * con_j*(1.0/2.0*(3.0*2.0*b_n))
		                                  * con_k*(1.0)
		                                  - 2.0*sqrt(2.0)*pow(1.0-c_n,1.0)
		                                  * con_i*(a_n)
		                                  * con_j*(1.0/2.0*(3.0*pow(b_n,2.0)-1.0))
		                                  * con_k*(1.0));

		get_scaling_basis_PYR(2,0,0,c_n,&con_i,&con_j,&con_k,&con_c);
		GradChiRef_rst[0][n*Nbf+11] = con*( 2.0*pow(1.0-c_n,1.0)
		                                  * con_i*(1.0/2.0*(3.0*2.0*a_n))
		                                  * con_j*(1.0)
		                                  * con_k*(1.0));
		GradChiRef_rst[1][n*Nbf+11] = 0.0;
		GradChiRef_rst[2][n*Nbf+11] = con*( sqrt(2.0)*a_n*pow(1.0-c_n,1.0)
		                                  * con_i*(1.0/2.0*(3.0*2.0*a_n))
		                                  * con_j*(1.0)
		                                  * con_k*(1.0)
		                                  - 2.0*sqrt(2.0)*pow(1.0-c_n,1.0)
		                                  * con_i*(1.0/2.0*(3.0*pow(a_n,2.0)-1.0))
		                                  * con_j*(1.0)
		                                  * con_k*(1.0));

		get_scaling_basis_PYR(2,1,0,c_n,&con_i,&con_j,&con_k,&con_c);
		GradChiRef_rst[0][n*Nbf+12] = con*( 2.0*pow(1.0-c_n,1.0)
		                                  * con_i*(1.0/2.0*(3.0*2.0*a_n))
		                                  * con_j*(b_n)
		                                  * con_k*(1.0));
		GradChiRef_rst[1][n*Nbf+12] = con*( 2.0*pow(1.0-c_n,1.0)
		                                  * con_i*(1.0/2.0*(3.0*pow(a_n,2.0)-1.0))
		                                  * con_j*(1.0)
		                                  * con_k*(1.0));
		GradChiRef_rst[2][n*Nbf+12] = con*( sqrt(2.0)*a_n*pow(1.0-c_n,1.0)
		                                  * con_i*(1.0/2.0*(3.0*2.0*a_n))
		                                  * con_j*(b_n)
		                                  * con_k*(1.0)
		                                  + sqrt(2.0)*b_n*pow(1.0-c_n,1.0)
		                                  * con_i*(1.0/2.0*(3.0*pow(a_n,2.0)-1.0))
		                                  * con_j*(1.0)
		                                  * con_k*(1.0)
		                                  - 2.0*sqrt(2.0)*pow(1.0-c_n,1.0)
		                                  * con_i*(1.0/2.0*(3.0*pow(a_n,2.0)-1.0))
		                                  * con_j*(b_n)
		                                  * con_k*(1.0));

		get_scaling_basis_PYR(2,2,0,c_n,&con_i,&con_j,&con_k,&con_c);
		GradChiRef_rst[0][n*Nbf+13] = con*( 2.0*pow(1.0-c_n,1.0)
		                                  * con_i*(1.0/2.0*(3.0*2.0*a_n))
		                                  * con_j*(1.0/2.0*(3.0*pow(b_n,2.0)-1.0))
		                                  * con_k*(1.0));
		GradChiRef_rst[1][n*Nbf+13] = con*( 2.0*pow(1.0-c_n,1.0)
		                                  * con_i*(1.0/2.0*(3.0*pow(a_n,2.0)-1.0))
		                                  * con_j*(1.0/2.0*(3.0*2.0*b_n))
		                                  * con_k*(1.0));
		GradChiRef_rst[2][n*Nbf+13] = con*( sqrt(2.0)*a_n*pow(1.0-c_n,1.0)
		                                  * con_i*(1.0/2.0*(3.0*2.0*a_n))
		                                  * con_j*(1.0/2.0*(3.0*pow(b_n,2.0)-1.0))
		                                  * con_k*(1.0)
		                                  + sqrt(2.0)*b_n*pow(1.0-c_n,1.0)
		                                  * con_i*(1.0/2.0*(3.0*pow(a_n,2.0)-1.0))
		                                  * con_j*(1.0/2.0*(3.0*2.0*b_n))
		                                  * con_k*(1.0)
		                                  - 2.0*sqrt(2.0)*pow(1.0-c_n,1.0)
		                                  * con_i*(1.0/2.0*(3.0*pow(a_n,2.0)-1.0))
		                                  * con_j*(1.0/2.0*(3.0*pow(b_n,2.0)-1.0))
		                                  * con_k*(1.0));
	}

	free(a);
	free(b);
	free(c);

	return GradChiRef_rst;
}

void test_unit_grad_basis_PYR(void)
{
	unsigned int pass;

	/*
	 *	grad_basis_PYR:
	 *
	 *		Input:
	 *
	 *			P, rst, Nn
	 *
	 *		Expected output:
	 *
	 *			P = 2:
	 *				GradChiRef_rst[0] = @(r,s,t) [ See grad_basis_PYR32[0] ]
	 *				GradChiRef_rst[1] = @(r,s,t) [ See grad_basis_PYR32[1] ]
	 *				GradChiRef_rst[2] = @(r,s,t) [ See grad_basis_PYR32[2] ]
	 */

	unsigned int d, Nn, Nbf, P, Prst;
	double *rst;

	struct S_CUBATURE *CUBDATA = malloc(sizeof *CUBDATA); // free

	d = 3;

	double **GradChiRef32_code, **GradChiRef32_test;

	P = 2;
	Prst = 4;

	set_cubdata(CUBDATA,false,false,"WV",d,Prst,cubature_PYR); // free
	set_from_cubdata(CUBDATA,&Nn,NULL,&rst,NULL,NULL);

	GradChiRef32_code = grad_basis_PYR(P,rst,Nn,&Nbf,d); // free
	GradChiRef32_test = grad_basis_PYR32(rst,Nn);        // free

/*
for (unsigned int dim = 2; dim < d; dim++) {
array_print_d(Nn,Nbf,GradChiRef32_code[dim],'R');
array_print_d(Nn,Nbf,GradChiRef32_test[dim],'R');
}
exit(1);
*/

	pass = 0;
	if (array_norm_diff_d(Nn*Nbf,GradChiRef32_code[0],GradChiRef32_test[0],"Inf") < EPS*10 &&
	    array_norm_diff_d(Nn*Nbf,GradChiRef32_code[1],GradChiRef32_test[1],"Inf") < EPS*10 &&
	    array_norm_diff_d(Nn*Nbf,GradChiRef32_code[2],GradChiRef32_test[2],"Inf") < EPS*10)
			pass = 1;

	test_print2(pass,"grad_basis_PYR (P2):");

	free(rst);
	array_free2_d(d,GradChiRef32_code);
	array_free2_d(d,GradChiRef32_test);

	/*
	 *	grad_basis_PYR, derivatives:
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

	d = 3;

	set_cubdata(CUBDATA,false,false,"GLL",d,P,cubature_PYR); // free
	set_from_cubdata(CUBDATA,&Nn,NULL,&rst,NULL,NULL);

	r = malloc(Nn * sizeof *r); // free
	s = malloc(Nn * sizeof *s); // free
	t = malloc(Nn * sizeof *t); // free
	for (i = 0; i < Nn; i++) {
		r[i] = rst[0*Nn+i];
		s[i] = rst[1*Nn+i];
		t[i] = rst[2*Nn+i];
	}

	I = identity_d(Nn); // free
	ChiRef_rst = basis_PYR(P,rst,Nn,&Nbf,d); // free
	ChiRefInv_rst = inverse_d(Nn,Nn,ChiRef_rst,I); // free

	GradChiRef_rst = grad_basis_PYR(P,rst,Nn,&Nbf,d); // free

	poly2(r,s,t,Nn,&f,&f_r,&f_s,&f_t); // free

	f_hat = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,1,Nn,1.0,ChiRefInv_rst,f); // free
	f_rcomp = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,1,Nn,1.0,GradChiRef_rst[0],f_hat); // free
	f_scomp = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,1,Nn,1.0,GradChiRef_rst[1],f_hat); // free
	f_tcomp = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,Nn,1,Nn,1.0,GradChiRef_rst[2],f_hat); // free

	pass = 0;
	if (array_norm_diff_d(Nn,f_r,f_rcomp,"Inf") < EPS*1e3 &&
	    array_norm_diff_d(Nn,f_s,f_scomp,"Inf") < EPS*1e3 &&
	    array_norm_diff_d(Nn,f_t,f_tcomp,"Inf") < EPS*1e3) {
			pass = 1;
	} else {
		printf("%e\n",array_norm_diff_d(Nn,f_r,f_rcomp,"Inf"));
		printf("%e\n",array_norm_diff_d(Nn,f_s,f_scomp,"Inf"));
		printf("%e\n",array_norm_diff_d(Nn,f_t,f_tcomp,"Inf"));
	}

	test_print2(pass,"               derivatives (d3):");

	free(rst), free(r), free(s), free(t);
	free(I);
	free(ChiRef_rst), free(ChiRefInv_rst);
	array_free2_d(d,GradChiRef_rst);
	free(f), free(f_hat);
	free(f_r), free(f_s), free(f_t);
	free(f_rcomp), free(f_scomp), free(f_tcomp);

	free(CUBDATA);
}
