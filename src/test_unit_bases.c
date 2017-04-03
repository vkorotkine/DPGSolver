// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "test_unit_bases.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#include "mkl.h"

#include "Parameters.h"
#include "Macros.h"
#include "Test.h"
#include "S_DB.h"

#include "test_support.h"
#include "test_code_bases.h"
#include "select_functions.h"
#include "bases.h"
#include "cubature.h"
#include "math_functions.h"
#include "matrix_functions.h"
#include "element_functions.h"
#include "memory_free.h"
#include "array_norm.h"
#include "array_print.h"

#undef I // No complex variables used here

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

struct S_orthogonality {
	unsigned int d, P, PIv, type;
	char         NodeType[STRLEN_MIN], PrntName[STRLEN_MAX];
};


static void test_unit_basis_TP_modal(void);
static void test_unit_basis_SI_modal(void);
static void test_unit_basis_PYR_modal(void);

static void test_unit_basis_TP_Bezier(void);
static void test_unit_basis_SI_Bezier(void);

void test_unit_basis_TP(void)
{
	test_unit_basis_TP_modal();
	test_unit_basis_TP_Bezier();
}

void test_unit_basis_SI(void)
{
	test_unit_basis_SI_modal();

	DB.NP = 8;
	initialize_ELEMENTs();
	test_unit_basis_SI_Bezier();
	memory_free_ELEMENTs();

	printf("\nWarning: Simplex Bezier basis for d = 3 is not being tested.\n\n"); TestDB.Nwarnings++;
}

void test_unit_basis_PYR(void)
{
	test_unit_basis_PYR_modal();
	printf("\nWarning: Pyramidal Bezier basis is not being tested.\n\n"); TestDB.Nwarnings++;
}

static void test_basis_orthogonality(struct S_orthogonality *data, const unsigned int no_last_entry)
{
	// Standard datatypes
	unsigned int pass = 0;

	char         *NodeType, *PrntName;
	unsigned int d, P, Nn, Ns, Nbf, *symms;
	double       *rst, *w, *ChiRef_rst, *WChiRef_rst, *M, *I;

	// Function pointers
	cubature_tdef cubature;
	basis_tdef    basis;

	select_functions_cubature(&cubature,data->type);
	select_functions_basis(&basis,data->type);

	d        = data->d;
	P        = data->P;
	NodeType = data->NodeType;
	PrntName = data->PrntName;

	cubature(&rst,&w,&symms,&Nn,&Ns,1,data->PIv,d,NodeType); // free

	ChiRef_rst = basis(P,rst,Nn,&Nbf,d); // free

	WChiRef_rst = calloc(Nn*Nbf , sizeof *WChiRef_rst); // free
	mm_diag_d(Nn,Nbf,w,ChiRef_rst,WChiRef_rst,1.0,0.0,'L','R');

	M = mm_Alloc_d(CBRM,CBT,CBNT,Nbf,Nbf,Nn,1.0,ChiRef_rst,WChiRef_rst); // free
	I = identity_d(Nbf); // free

	if (no_last_entry < 2) {
		if (array_norm_diff_d(pow(Nbf,2)-no_last_entry,M,I,"Inf") < EPS*1e3)
			pass = 1;

		char *PrintName = malloc(STRLEN_MAX * sizeof *PrintName); // free

		sprintf(PrintName,"          orthogonality - %s (d%d, P%d):",PrntName,d,P);
		test_print2(pass,PrintName);

		free(PrintName);
	} else if (no_last_entry < 4) {
		if (no_last_entry == 3) {
			array_print_d(Nbf,Nbf,M,'R');
			array_print_d(Nbf,Nbf,I,'R');
		}
	} else {
		EXIT_UNSUPPORTED;
	}

	free(rst);
	free(w);
	free(symms);
	free(ChiRef_rst);
	free(WChiRef_rst);
	free(M);
	free(I);
}

static double *basis_TP13(const double *rst, const unsigned int Nn)
{
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

static void test_unit_basis_TP_modal(void)
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
	 *				ChiRef_rst = @(r) [ See basis_TP13 ]
	 */

	unsigned int dE, Nn, Ns, Nbf, P, Prst;
	unsigned int *symms;
	double *rst, *w;

	dE = 1;

	double *ChiRef13_code, *ChiRef13_test;

	P = 3;
	Prst = 4;

	cubature_TP(&rst,&w,&symms,&Nn,&Ns,0,Prst,dE,"GL"); // free

	ChiRef13_code = basis_TP(P,rst,Nn,&Nbf,dE); // free
	ChiRef13_test = basis_TP13(rst,Nn); // free

	pass = 0;
	if (array_norm_diff_d(Nn*(P+1),ChiRef13_code,ChiRef13_test,"Inf") < EPS*10)
		pass = 1;

	test_print2(pass,"basis_TP (d1, P3):");

	free(rst);
	free(symms);
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
	 *				ChiRef_rst = @(r,s) [ See basis_TP22 ]
	 */

	dE = 2;

	double *ChiRef22_code, *ChiRef22_test;

	P = 2;
	Prst = 4;

	cubature_TP(&rst,&w,&symms,&Nn,&Ns,0,Prst,dE,"GL"); // free

	ChiRef22_code = basis_TP(P,rst,Nn,&Nbf,dE); // free
	ChiRef22_test = basis_TP22(rst,Nn); // free

	pass = 0;
	if (array_norm_diff_d(Nn*pow(P+1,dE),ChiRef22_code,ChiRef22_test,"Inf") < EPS*10)
		pass = 1;

	test_print2(pass,"         (d2, P2):");

	free(rst);
	free(symms);
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
	 *				ChiRef_rst = @(r,s,t) [ See basis_TP31 ]
	 */

	dE = 3;

	double *ChiRef31_code, *ChiRef31_test;

	P = 1;
	Prst = 4;

	cubature_TP(&rst,&w,&symms,&Nn,&Ns,0,Prst,dE,"GL"); // free

	ChiRef31_code = basis_TP(P,rst,Nn,&Nbf,dE); // free
	ChiRef31_test = basis_TP31(rst,Nn); // free

	pass = 0;
	if (array_norm_diff_d(Nn*pow(P+1,dE),ChiRef31_code,ChiRef31_test,"Inf") < EPS*10)
		pass = 1;

	test_print2(pass,"         (d3, P1):");

	free(rst);
	free(symms);
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

	struct S_orthogonality *data;

	data = malloc(sizeof *data); // free

	// dE = 1
	dE = 1;
	P  = 4;
	data->d = dE;
	data->P = P;
	data->type = LINE;

	// GL
	data->PIv = P;
	strcpy(data->NodeType,"GL");
	strcpy(data->PrntName,"GL ");
	test_basis_orthogonality(data,0);

	// GLL
	data->PIv = P;
	strcpy(data->NodeType,"GLL");
	strcpy(data->PrntName,"GLL");
	test_basis_orthogonality(data,1);


	// dE = 2
	dE = 2;
	P  = 4;
	data->d = dE;
	data->P = P;
	data->type = QUAD;

	// GL
	data->PIv = P;
	strcpy(data->NodeType,"GL");
	strcpy(data->PrntName,"GL ");
	test_basis_orthogonality(data,0);


	// dE = 3
	dE = 3;
	P  = 4;
	data->d = dE;
	data->P = P;
	data->type = HEX;

	// GL
	data->PIv = P;
	strcpy(data->NodeType,"GL");
	strcpy(data->PrntName,"GL ");
	test_basis_orthogonality(data,0);

	free(data);
}

/*
 *	Purpose:
 *		Test correctness of implementation of basis_SI.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

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

static void test_unit_basis_SI_modal(void)
{
	printf("\nWarning: Ensure that discrepancy between M and I for WSH nodes is as expected\n\n");
	TestDB.Nwarnings++;

	unsigned int pass;
	unsigned int printMI = 0;

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

	d = 2;

	double *ChiRef22_code, *ChiRef22_test;

	P = 2;
	Prst = 4;

	cubature_TRI(&rst,&w,&symms,&Nn,&Ns,0,Prst,d,"AO"); // free

	ChiRef22_code = basis_SI(P,rst,Nn,&Nbf,d); // free
	ChiRef22_test = basis_TRI2(rst,Nn);        // free

	pass = 0;
	if (array_norm_diff_d(Nn*Nbf,ChiRef22_code,ChiRef22_test,"Inf") < EPS*10)
		pass = 1;


	test_print2(pass,"basis_SI (d2, P2):");

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
		pass = 1;

	test_print2(pass,"         (d2, P3):");

	free(rst);
	free(symms);
	free(ChiRef23_code);
	free(ChiRef23_test);

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

	d = 3;

	double *ChiRef32_code, *ChiRef32_test;

	P = 2;
	Prst = 4;

	cubature_TET(&rst,&w,&symms,&Nn,&Ns,0,Prst,d,"AO"); // free

	ChiRef32_code = basis_SI(P,rst,Nn,&Nbf,d); // free
	ChiRef32_test = basis_TET2(rst,Nn);        // free

	pass = 0;
	if (array_norm_diff_d(Nn*Nbf,ChiRef32_code,ChiRef32_test,"Inf") < EPS*10)
		pass = 1;


	test_print2(pass,"         (d3, P2):");

	free(rst);
	free(symms);
	free(ChiRef32_code);
	free(ChiRef32_test);

	/*
	 *	basis_SI orthogonality:
	 *
	 *		Input:
	 *			rst_WSH (d = 2), rst_WSH (d = 3), rst_WV, P, d
	 *
	 *		Expected Output:
	 *
	 *			d = 2, P = 2, rst_WSH (P = 2), rst_WV (P = 4):
	 *				M = ChiRef_rst'*W*ChiRef_rst = I
	 *
	 *			d = 2, P = 3, rst_WSH (P = 3), rst_WV (P = 6):
	 *				M = ChiRef_rst'*W*ChiRef_rst = I with error in highest order entries (WSH)
	 *				M = ChiRef_rst'*W*ChiRef_rst = I (WV)
	 *
	 *			d = 3, P = 1, rst_WSH (P = 1), rst_WV (P = 2):
	 *				M = ChiRef_rst'*W*ChiRef_rst = I
	 *
	 *			d = 2, P = 2, rst_WSH (P = 2), rst_WV (P = 4):
	 *				M = ChiRef_rst'*W*ChiRef_rst = I with error in highest order entries (WSH)
	 *				M = ChiRef_rst'*W*ChiRef_rst = I (WV)
	 */

	struct S_orthogonality *data;

	data = malloc(sizeof *data); // free

	// d = 2
	d = 2;
	data->d = d;
	data->type = TRI;

	// WSH
	for (P = 2; P <= 3; P++) {
		data->P   = P;
		data->PIv = P;
		strcpy(data->NodeType,"WSH");
		strcpy(data->PrntName,"WSH");
		if (P == 2) {
			test_basis_orthogonality(data,0);
		} else if (P == 3) {
			if (!printMI)
				test_basis_orthogonality(data,2);
			else
				test_basis_orthogonality(data,3);
		}
	}

	// WV
	for (P = 2; P <= 9; P++) {
		data->P   = P;
		data->PIv = 2*P;
		strcpy(data->NodeType,"WV");
		strcpy(data->PrntName,"WV ");
		test_basis_orthogonality(data,0);
	}

	// d = 3
	d = 3;
	data->d = d;
	data->type = TET;

	for (P = 1; P <= 2; P++) {
		data->P = P;

		// WSH
		data->PIv = P;
		strcpy(data->NodeType,"WSH");
		strcpy(data->PrntName,"WSH");
		if (P == 1) {
			test_basis_orthogonality(data,0);
		} else if (P == 2) {
			if (!printMI)
				test_basis_orthogonality(data,2);
			else
				test_basis_orthogonality(data,3);
		}

		// WV
		data->PIv = 2*P;
		strcpy(data->NodeType,"WV");
		strcpy(data->PrntName,"WV ");
		test_basis_orthogonality(data,0);
	}

	free(data);
}

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

static void test_unit_basis_PYR_modal(void)
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
		pass = 1;

	test_print2(pass,"basis_PYR (P2):");

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

	struct S_orthogonality *data;

	data = malloc(sizeof *data); // free

	d = 3;
	P = 3;
	data->d = d;
	data->P = P;
	data->type = PYR;

	// GL (HEX To PYR)
	data->PIv = P+1;
	strcpy(data->NodeType,"GLW");
	strcpy(data->PrntName,"GL ");
	test_basis_orthogonality(data,0);

	// GLL (HEX To PYR)
	data->PIv = P+2;
	strcpy(data->NodeType,"GLLW");
	strcpy(data->PrntName,"GLL");
	test_basis_orthogonality(data,0);

	// GJ (HEX To PYR)
	data->PIv = P;
	strcpy(data->NodeType,"GJW");
	strcpy(data->PrntName,"GJ ");
	test_basis_orthogonality(data,0);

	free(data);
}

static double *basis_TP13_Bezier(const double *rst, const unsigned int Nn)
{
	unsigned int n, N, Nbf;
	double       *ChiBez_rst;
	const double *r_ptr;

	N   = 3+1;
	Nbf = pow(N,1);

	r_ptr = &rst[0*Nn];

	ChiBez_rst = malloc(Nn*Nbf * sizeof *ChiBez_rst); // keep (requires external free)

	for (n = 0; n < Nn; n++) {
		double r, b0, b1;

		r  = r_ptr[n];
		b0 = (1.0-r)/2.0;
		b1 = (1.0+r)/2.0;

		ChiBez_rst[n*Nbf+0] = pow(b0,3.0)*pow(b1,0.0);
		ChiBez_rst[n*Nbf+1] = 3.0*pow(b0,2.0)*pow(b1,1.0);
		ChiBez_rst[n*Nbf+2] = 3.0*pow(b0,1.0)*pow(b1,2.0);
		ChiBez_rst[n*Nbf+3] = pow(b0,0.0)*pow(b1,3.0);
	}

	return ChiBez_rst;
}

static double *basis_TP22_Bezier(const double *rst, const unsigned int Nn)
{
	unsigned int n, N, Nbf;
	double       *ChiBez_rst;
	const double *r_ptr, *s_ptr;

	N   = 2+1;
	Nbf = pow(N,2);

	r_ptr = &rst[0*Nn];
	s_ptr = &rst[1*Nn];

	ChiBez_rst = malloc(Nn*Nbf * sizeof *ChiBez_rst); // keep (requires external free)

	for (n = 0; n < Nn; n++) {
		double r, s, b0r, b1r, b0s, b1s;

		r   = r_ptr[n];
		b0r = (1.0-r)/2.0;
		b1r = (1.0+r)/2.0;

		s   = s_ptr[n];
		b0s = (1.0-s)/2.0;
		b1s = (1.0+s)/2.0;

		ChiBez_rst[n*Nbf+0] = pow(b0r,2.0)*pow(b1r,0.0)     * pow(b0s,2.0)*pow(b1s,0.0)    ;
		ChiBez_rst[n*Nbf+1] = 2.0*pow(b0r,1.0)*pow(b1r,1.0) * pow(b0s,2.0)*pow(b1s,0.0)    ;
		ChiBez_rst[n*Nbf+2] = pow(b0r,0.0)*pow(b1r,2.0)     * pow(b0s,2.0)*pow(b1s,0.0)    ;
		ChiBez_rst[n*Nbf+3] = pow(b0r,2.0)*pow(b1r,0.0)     * 2.0*pow(b0s,1.0)*pow(b1s,1.0);
		ChiBez_rst[n*Nbf+4] = 2.0*pow(b0r,1.0)*pow(b1r,1.0) * 2.0*pow(b0s,1.0)*pow(b1s,1.0);
		ChiBez_rst[n*Nbf+5] = pow(b0r,0.0)*pow(b1r,2.0)     * 2.0*pow(b0s,1.0)*pow(b1s,1.0);
		ChiBez_rst[n*Nbf+6] = pow(b0r,2.0)*pow(b1r,0.0)     * pow(b0s,0.0)*pow(b1s,2.0)    ;
		ChiBez_rst[n*Nbf+7] = 2.0*pow(b0r,1.0)*pow(b1r,1.0) * pow(b0s,0.0)*pow(b1s,2.0)    ;
		ChiBez_rst[n*Nbf+8] = pow(b0r,0.0)*pow(b1r,2.0)     * pow(b0s,0.0)*pow(b1s,2.0)    ;
	}

	return ChiBez_rst;
}

static double *basis_TP32_Bezier(const double *rst, const unsigned int Nn)
{
	unsigned int n, N, Nbf;
	double       *ChiBez_rst;
	const double *r_ptr, *s_ptr, *t_ptr;

	N   = 2+1;
	Nbf = pow(N,3);

	r_ptr = &rst[0*Nn];
	s_ptr = &rst[1*Nn];
	t_ptr = &rst[2*Nn];

	ChiBez_rst = malloc(Nn*Nbf * sizeof *ChiBez_rst); // keep (requires external free)

	for (n = 0; n < Nn; n++) {
		double r, s, t, b0r, b1r, b0s, b1s, b0t, b1t;

		r   = r_ptr[n];
		b0r = (1.0-r)/2.0;
		b1r = (1.0+r)/2.0;

		s   = s_ptr[n];
		b0s = (1.0-s)/2.0;
		b1s = (1.0+s)/2.0;

		t   = t_ptr[n];
		b0t = (1.0-t)/2.0;
		b1t = (1.0+t)/2.0;

		ChiBez_rst[n*Nbf+0 ] =     pow(b0r,2.0)*pow(b1r,0.0) *     pow(b0s,2.0)*pow(b1s,0.0) *     pow(b0t,2.0)*pow(b1t,0.0);
		ChiBez_rst[n*Nbf+1 ] = 2.0*pow(b0r,1.0)*pow(b1r,1.0) *     pow(b0s,2.0)*pow(b1s,0.0) *     pow(b0t,2.0)*pow(b1t,0.0);
		ChiBez_rst[n*Nbf+2 ] =     pow(b0r,0.0)*pow(b1r,2.0) *     pow(b0s,2.0)*pow(b1s,0.0) *     pow(b0t,2.0)*pow(b1t,0.0);
		ChiBez_rst[n*Nbf+3 ] =     pow(b0r,2.0)*pow(b1r,0.0) * 2.0*pow(b0s,1.0)*pow(b1s,1.0) *     pow(b0t,2.0)*pow(b1t,0.0);
		ChiBez_rst[n*Nbf+4 ] = 2.0*pow(b0r,1.0)*pow(b1r,1.0) * 2.0*pow(b0s,1.0)*pow(b1s,1.0) *     pow(b0t,2.0)*pow(b1t,0.0);
		ChiBez_rst[n*Nbf+5 ] =     pow(b0r,0.0)*pow(b1r,2.0) * 2.0*pow(b0s,1.0)*pow(b1s,1.0) *     pow(b0t,2.0)*pow(b1t,0.0);
		ChiBez_rst[n*Nbf+6 ] =     pow(b0r,2.0)*pow(b1r,0.0) *     pow(b0s,0.0)*pow(b1s,2.0) *     pow(b0t,2.0)*pow(b1t,0.0);
		ChiBez_rst[n*Nbf+7 ] = 2.0*pow(b0r,1.0)*pow(b1r,1.0) *     pow(b0s,0.0)*pow(b1s,2.0) *     pow(b0t,2.0)*pow(b1t,0.0);
		ChiBez_rst[n*Nbf+8 ] =     pow(b0r,0.0)*pow(b1r,2.0) *     pow(b0s,0.0)*pow(b1s,2.0) *     pow(b0t,2.0)*pow(b1t,0.0);
		ChiBez_rst[n*Nbf+9 ] =     pow(b0r,2.0)*pow(b1r,0.0) *     pow(b0s,2.0)*pow(b1s,0.0) * 2.0*pow(b0t,1.0)*pow(b1t,1.0);
		ChiBez_rst[n*Nbf+10] = 2.0*pow(b0r,1.0)*pow(b1r,1.0) *     pow(b0s,2.0)*pow(b1s,0.0) * 2.0*pow(b0t,1.0)*pow(b1t,1.0);
		ChiBez_rst[n*Nbf+11] =     pow(b0r,0.0)*pow(b1r,2.0) *     pow(b0s,2.0)*pow(b1s,0.0) * 2.0*pow(b0t,1.0)*pow(b1t,1.0);
		ChiBez_rst[n*Nbf+12] =     pow(b0r,2.0)*pow(b1r,0.0) * 2.0*pow(b0s,1.0)*pow(b1s,1.0) * 2.0*pow(b0t,1.0)*pow(b1t,1.0);
		ChiBez_rst[n*Nbf+13] = 2.0*pow(b0r,1.0)*pow(b1r,1.0) * 2.0*pow(b0s,1.0)*pow(b1s,1.0) * 2.0*pow(b0t,1.0)*pow(b1t,1.0);
		ChiBez_rst[n*Nbf+14] =     pow(b0r,0.0)*pow(b1r,2.0) * 2.0*pow(b0s,1.0)*pow(b1s,1.0) * 2.0*pow(b0t,1.0)*pow(b1t,1.0);
		ChiBez_rst[n*Nbf+15] =     pow(b0r,2.0)*pow(b1r,0.0) *     pow(b0s,0.0)*pow(b1s,2.0) * 2.0*pow(b0t,1.0)*pow(b1t,1.0);
		ChiBez_rst[n*Nbf+16] = 2.0*pow(b0r,1.0)*pow(b1r,1.0) *     pow(b0s,0.0)*pow(b1s,2.0) * 2.0*pow(b0t,1.0)*pow(b1t,1.0);
		ChiBez_rst[n*Nbf+17] =     pow(b0r,0.0)*pow(b1r,2.0) *     pow(b0s,0.0)*pow(b1s,2.0) * 2.0*pow(b0t,1.0)*pow(b1t,1.0);
		ChiBez_rst[n*Nbf+18] =     pow(b0r,2.0)*pow(b1r,0.0) *     pow(b0s,2.0)*pow(b1s,0.0) *     pow(b0t,0.0)*pow(b1t,2.0);
		ChiBez_rst[n*Nbf+19] = 2.0*pow(b0r,1.0)*pow(b1r,1.0) *     pow(b0s,2.0)*pow(b1s,0.0) *     pow(b0t,0.0)*pow(b1t,2.0);
		ChiBez_rst[n*Nbf+20] =     pow(b0r,0.0)*pow(b1r,2.0) *     pow(b0s,2.0)*pow(b1s,0.0) *     pow(b0t,0.0)*pow(b1t,2.0);
		ChiBez_rst[n*Nbf+21] =     pow(b0r,2.0)*pow(b1r,0.0) * 2.0*pow(b0s,1.0)*pow(b1s,1.0) *     pow(b0t,0.0)*pow(b1t,2.0);
		ChiBez_rst[n*Nbf+22] = 2.0*pow(b0r,1.0)*pow(b1r,1.0) * 2.0*pow(b0s,1.0)*pow(b1s,1.0) *     pow(b0t,0.0)*pow(b1t,2.0);
		ChiBez_rst[n*Nbf+23] =     pow(b0r,0.0)*pow(b1r,2.0) * 2.0*pow(b0s,1.0)*pow(b1s,1.0) *     pow(b0t,0.0)*pow(b1t,2.0);
		ChiBez_rst[n*Nbf+24] =     pow(b0r,2.0)*pow(b1r,0.0) *     pow(b0s,0.0)*pow(b1s,2.0) *     pow(b0t,0.0)*pow(b1t,2.0);
		ChiBez_rst[n*Nbf+25] = 2.0*pow(b0r,1.0)*pow(b1r,1.0) *     pow(b0s,0.0)*pow(b1s,2.0) *     pow(b0t,0.0)*pow(b1t,2.0);
		ChiBez_rst[n*Nbf+26] =     pow(b0r,0.0)*pow(b1r,2.0) *     pow(b0s,0.0)*pow(b1s,2.0) *     pow(b0t,0.0)*pow(b1t,2.0);
	}

	return ChiBez_rst;
}

typedef double *(*basis_test_tdef) (const double *rst, const unsigned int Nn);

/*
 *	Purpose:
 *		Test correctness of implementation of basis_TP_Bezier.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

static void test_unit_basis_TP_Bezier(void)
{
	unsigned int pass;

	char *PrintName = malloc(STRLEN_MAX * sizeof *PrintName); // free

	/*
	 *	basis_TP_Bezier
	 *
	 *		Input:
	 *
	 *			P, rst, Nn, dE.
	 *
	 *		Expected output:
	 *
	 *			ChiBez_rst = @(r)     [ See basis_TP13_Bezier ]
	 *			ChiBez_rst = @(r,s)   [ See basis_TP22_Bezier ]
	 *			ChiBez_rst = @(r,s,t) [ See basis_TP32_Bezier ]
	 */

	for (unsigned int dE = 1; dE <= DMAX; dE++) {
		unsigned int P, Prst, Nn, Ns, Nbf, *symms;
		double       *rst, *w, *ChiBez_code, *ChiBez_test;

		basis_test_tdef basis_test;

		Prst = 4;
		if (dE == 1) {
			P = 3;
			basis_test = basis_TP13_Bezier;
		} else if (dE == 2) {
			P = 2;
			basis_test = basis_TP22_Bezier;
		} else if (dE == 3) {
			P = 2;
			basis_test = basis_TP32_Bezier;
		} else {
			EXIT_UNSUPPORTED;
		}
		cubature_TP(&rst,&w,&symms,&Nn,&Ns,0,Prst,dE,"GLL"); // free
		ChiBez_code = basis_TP_Bezier(P,rst,Nn,&Nbf,dE);
		ChiBez_test = basis_test(rst,Nn);

		pass = 0;
		if (array_norm_diff_d(Nn*pow(P+1,dE),ChiBez_code,ChiBez_test,"Inf") < EPS*10)
			pass = 1;

		if (dE == 1)
			sprintf(PrintName,"basis_TP_Bezier (d%d, P%d):",dE,P);
		else
			sprintf(PrintName,"                (d%d, P%d):",dE,P);

		test_print2(pass,PrintName);

		free(rst);
		free(symms);
		free(ChiBez_test);
		free(ChiBez_code);
	}

	/*
	 *	basis_TP_Bezier partition of unity/positivity:
	 *
	 *		Input:
	 *
	 *			rst
	 *
	 *		Expected Output:
	 *
	 *			Basis functions at each node sum to one.
	 *			Basis functions evaluated at each node have values greater than or equal to zero.
	 */

	for (unsigned int dE = 1; dE <= DMAX; dE++) {
		unsigned int P, Prst, Nn, Ns, Nbf, *symms;
		double       *rst, *w, *ChiBez;

		Prst = 4;
		P    = 3;

		cubature_TP(&rst,&w,&symms,&Nn,&Ns,0,Prst,dE,"GLL"); // free
		ChiBez = basis_TP_Bezier(P,rst,Nn,&Nbf,dE);

		double *RowSum, *ones;
		RowSum = calloc(Nn , sizeof *RowSum);
		ones   = calloc(Nn , sizeof *RowSum);
		for (size_t i = 0; i < Nn; i++) {
			ones[i] = 1.0;
			for (size_t j = 0; j < Nbf; j++) {
				RowSum[i] += ChiBez[j+i*Nbf];
			}
		}

		bool   Positive = 1;
		double minVal   = 1e15;
		for (size_t i = 0, iMax = Nn*Nbf; i < iMax; i++) {
			if (ChiBez[i] < minVal)
				minVal = ChiBez[i];
		}

		if (minVal < 0.0)
			Positive = 0;

		pass = 0;
		if (array_norm_diff_d(Nn,RowSum,ones,"Inf") < EPS*10 && Positive)
			pass = 1;

		if (dE == 1) {
			sprintf(PrintName,"basis_TP_Bezier partition of unity (d%d, P%d):",dE,P);
		} else {
			sprintf(PrintName,"                                   (d%d, P%d):",dE,P);
		}
		test_print2(pass,PrintName);

		free(rst);
		free(symms);
		free(ChiBez);
		free(RowSum);
		free(ones);
	}

	free(PrintName);
}

static void get_BCoord_Exponents_test(const unsigned int P, const unsigned int d, unsigned int *NExp,
                                      unsigned int **NpermsOut, unsigned int **ExponentsOut)
{
	unsigned int Nbf, *Nperms, *Exp;

	Nbf = (unsigned int) (factorial_ull(d+P)/(factorial_ull(d)*factorial_ull(P)));

	Nperms = malloc(Nbf       * sizeof *Nperms); // keep
	Exp   = malloc(Nbf*(d+1) * sizeof *Exp);   // keep

	size_t Ind = 0;
	if (d == 2) {
		if (P == 0) {
			Exp[Ind*(d+1)+0] = 0; Exp[Ind*(d+1)+1] = 0; Exp[Ind*(d+1)+2] = 0; Nperms[Ind++] = 1;
		} else if (P == 1) {
			Exp[Ind*(d+1)+0] = 1; Exp[Ind*(d+1)+1] = 0; Exp[Ind*(d+1)+2] = 0; Nperms[Ind++] = 3;
		} else if (P == 2) {
			Exp[Ind*(d+1)+0] = 1; Exp[Ind*(d+1)+1] = 1; Exp[Ind*(d+1)+2] = 0; Nperms[Ind++] = 3;
			Exp[Ind*(d+1)+0] = 2; Exp[Ind*(d+1)+1] = 0; Exp[Ind*(d+1)+2] = 0; Nperms[Ind++] = 3;
		} else if (P == 3) {
			Exp[Ind*(d+1)+0] = 1; Exp[Ind*(d+1)+1] = 1; Exp[Ind*(d+1)+2] = 1; Nperms[Ind++] = 1;
			Exp[Ind*(d+1)+0] = 2; Exp[Ind*(d+1)+1] = 1; Exp[Ind*(d+1)+2] = 0; Nperms[Ind++] = 6;
			Exp[Ind*(d+1)+0] = 3; Exp[Ind*(d+1)+1] = 0; Exp[Ind*(d+1)+2] = 0; Nperms[Ind++] = 3;
		} else if (P == 4) {
			Exp[Ind*(d+1)+0] = 2; Exp[Ind*(d+1)+1] = 1; Exp[Ind*(d+1)+2] = 1; Nperms[Ind++] = 3;
			Exp[Ind*(d+1)+0] = 2; Exp[Ind*(d+1)+1] = 2; Exp[Ind*(d+1)+2] = 0; Nperms[Ind++] = 3;
			Exp[Ind*(d+1)+0] = 3; Exp[Ind*(d+1)+1] = 1; Exp[Ind*(d+1)+2] = 0; Nperms[Ind++] = 6;
			Exp[Ind*(d+1)+0] = 4; Exp[Ind*(d+1)+1] = 0; Exp[Ind*(d+1)+2] = 0; Nperms[Ind++] = 3;
		} else if (P == 5) {
			Exp[Ind*(d+1)+0] = 2; Exp[Ind*(d+1)+1] = 2; Exp[Ind*(d+1)+2] = 1; Nperms[Ind++] = 3;
			Exp[Ind*(d+1)+0] = 3; Exp[Ind*(d+1)+1] = 1; Exp[Ind*(d+1)+2] = 1; Nperms[Ind++] = 3;
			Exp[Ind*(d+1)+0] = 3; Exp[Ind*(d+1)+1] = 2; Exp[Ind*(d+1)+2] = 0; Nperms[Ind++] = 6;
			Exp[Ind*(d+1)+0] = 4; Exp[Ind*(d+1)+1] = 1; Exp[Ind*(d+1)+2] = 0; Nperms[Ind++] = 6;
			Exp[Ind*(d+1)+0] = 5; Exp[Ind*(d+1)+1] = 0; Exp[Ind*(d+1)+2] = 0; Nperms[Ind++] = 3;
		} else if (P == 6) {
			Exp[Ind*(d+1)+0] = 2; Exp[Ind*(d+1)+1] = 2; Exp[Ind*(d+1)+2] = 2; Nperms[Ind++] = 1;
			Exp[Ind*(d+1)+0] = 3; Exp[Ind*(d+1)+1] = 2; Exp[Ind*(d+1)+2] = 1; Nperms[Ind++] = 6;
			Exp[Ind*(d+1)+0] = 3; Exp[Ind*(d+1)+1] = 3; Exp[Ind*(d+1)+2] = 0; Nperms[Ind++] = 3;
			Exp[Ind*(d+1)+0] = 4; Exp[Ind*(d+1)+1] = 1; Exp[Ind*(d+1)+2] = 1; Nperms[Ind++] = 3;
			Exp[Ind*(d+1)+0] = 4; Exp[Ind*(d+1)+1] = 2; Exp[Ind*(d+1)+2] = 0; Nperms[Ind++] = 6;
			Exp[Ind*(d+1)+0] = 5; Exp[Ind*(d+1)+1] = 1; Exp[Ind*(d+1)+2] = 0; Nperms[Ind++] = 6;
			Exp[Ind*(d+1)+0] = 6; Exp[Ind*(d+1)+1] = 0; Exp[Ind*(d+1)+2] = 0; Nperms[Ind++] = 3;
		} else {
			EXIT_UNSUPPORTED;
		}
	} else if (d == 3) {
		if (P == 0) {
			Exp[Ind*(d+1)+0] = 0; Exp[Ind*(d+1)+1] = 0; Exp[Ind*(d+1)+2] = 0; Exp[Ind*(d+1)+3] = 0; Nperms[Ind++] = 1;
		} else if (P == 1) {
			Exp[Ind*(d+1)+0] = 1; Exp[Ind*(d+1)+1] = 0; Exp[Ind*(d+1)+2] = 0; Exp[Ind*(d+1)+3] = 0; Nperms[Ind++] = 4;
		} else if (P == 2) {
			Exp[Ind*(d+1)+0] = 1; Exp[Ind*(d+1)+1] = 1; Exp[Ind*(d+1)+2] = 0; Exp[Ind*(d+1)+3] = 0; Nperms[Ind++] = 6;
			Exp[Ind*(d+1)+0] = 2; Exp[Ind*(d+1)+1] = 0; Exp[Ind*(d+1)+2] = 0; Exp[Ind*(d+1)+3] = 0; Nperms[Ind++] = 4;
		} else if (P == 3) {
			Exp[Ind*(d+1)+0] = 1; Exp[Ind*(d+1)+1] = 1; Exp[Ind*(d+1)+2] = 1; Exp[Ind*(d+1)+3] = 0; Nperms[Ind++] = 4;
			Exp[Ind*(d+1)+0] = 2; Exp[Ind*(d+1)+1] = 1; Exp[Ind*(d+1)+2] = 0; Exp[Ind*(d+1)+3] = 0; Nperms[Ind++] = 12;
			Exp[Ind*(d+1)+0] = 3; Exp[Ind*(d+1)+1] = 0; Exp[Ind*(d+1)+2] = 0; Exp[Ind*(d+1)+3] = 0; Nperms[Ind++] = 4;
		} else if (P == 4) {
			Exp[Ind*(d+1)+0] = 1; Exp[Ind*(d+1)+1] = 1; Exp[Ind*(d+1)+2] = 1; Exp[Ind*(d+1)+3] = 1; Nperms[Ind++] = 1;
			Exp[Ind*(d+1)+0] = 2; Exp[Ind*(d+1)+1] = 1; Exp[Ind*(d+1)+2] = 1; Exp[Ind*(d+1)+3] = 0; Nperms[Ind++] = 12;
			Exp[Ind*(d+1)+0] = 2; Exp[Ind*(d+1)+1] = 2; Exp[Ind*(d+1)+2] = 0; Exp[Ind*(d+1)+3] = 0; Nperms[Ind++] = 6;
			Exp[Ind*(d+1)+0] = 3; Exp[Ind*(d+1)+1] = 1; Exp[Ind*(d+1)+2] = 0; Exp[Ind*(d+1)+3] = 0; Nperms[Ind++] = 12;
			Exp[Ind*(d+1)+0] = 4; Exp[Ind*(d+1)+1] = 0; Exp[Ind*(d+1)+2] = 0; Exp[Ind*(d+1)+3] = 0; Nperms[Ind++] = 4;
		} else if (P == 5) {
			Exp[Ind*(d+1)+0] = 2; Exp[Ind*(d+1)+1] = 1; Exp[Ind*(d+1)+2] = 1; Exp[Ind*(d+1)+3] = 1; Nperms[Ind++] = 4;
			Exp[Ind*(d+1)+0] = 2; Exp[Ind*(d+1)+1] = 2; Exp[Ind*(d+1)+2] = 1; Exp[Ind*(d+1)+3] = 0; Nperms[Ind++] = 12;
			Exp[Ind*(d+1)+0] = 3; Exp[Ind*(d+1)+1] = 1; Exp[Ind*(d+1)+2] = 1; Exp[Ind*(d+1)+3] = 0; Nperms[Ind++] = 12;
			Exp[Ind*(d+1)+0] = 3; Exp[Ind*(d+1)+1] = 2; Exp[Ind*(d+1)+2] = 0; Exp[Ind*(d+1)+3] = 0; Nperms[Ind++] = 12;
			Exp[Ind*(d+1)+0] = 4; Exp[Ind*(d+1)+1] = 1; Exp[Ind*(d+1)+2] = 0; Exp[Ind*(d+1)+3] = 0; Nperms[Ind++] = 12;
			Exp[Ind*(d+1)+0] = 5; Exp[Ind*(d+1)+1] = 0; Exp[Ind*(d+1)+2] = 0; Exp[Ind*(d+1)+3] = 0; Nperms[Ind++] = 4;
		} else if (P == 6) {
			Exp[Ind*(d+1)+0] = 2; Exp[Ind*(d+1)+1] = 2; Exp[Ind*(d+1)+2] = 1; Exp[Ind*(d+1)+3] = 1; Nperms[Ind++] = 6;
			Exp[Ind*(d+1)+0] = 2; Exp[Ind*(d+1)+1] = 2; Exp[Ind*(d+1)+2] = 2; Exp[Ind*(d+1)+3] = 0; Nperms[Ind++] = 4;
			Exp[Ind*(d+1)+0] = 3; Exp[Ind*(d+1)+1] = 1; Exp[Ind*(d+1)+2] = 1; Exp[Ind*(d+1)+3] = 1; Nperms[Ind++] = 4;
			Exp[Ind*(d+1)+0] = 3; Exp[Ind*(d+1)+1] = 2; Exp[Ind*(d+1)+2] = 1; Exp[Ind*(d+1)+3] = 0; Nperms[Ind++] = 24;
			Exp[Ind*(d+1)+0] = 3; Exp[Ind*(d+1)+1] = 3; Exp[Ind*(d+1)+2] = 0; Exp[Ind*(d+1)+3] = 0; Nperms[Ind++] = 6;
			Exp[Ind*(d+1)+0] = 4; Exp[Ind*(d+1)+1] = 1; Exp[Ind*(d+1)+2] = 1; Exp[Ind*(d+1)+3] = 0; Nperms[Ind++] = 12;
			Exp[Ind*(d+1)+0] = 4; Exp[Ind*(d+1)+1] = 2; Exp[Ind*(d+1)+2] = 0; Exp[Ind*(d+1)+3] = 0; Nperms[Ind++] = 12;
			Exp[Ind*(d+1)+0] = 5; Exp[Ind*(d+1)+1] = 1; Exp[Ind*(d+1)+2] = 0; Exp[Ind*(d+1)+3] = 0; Nperms[Ind++] = 12;
			Exp[Ind*(d+1)+0] = 6; Exp[Ind*(d+1)+1] = 0; Exp[Ind*(d+1)+2] = 0; Exp[Ind*(d+1)+3] = 0; Nperms[Ind++] = 4;
		} else {
			EXIT_UNSUPPORTED;
		}
	} else {
		EXIT_UNSUPPORTED;
	}

	*NExp         = Ind;
	*NpermsOut    = Nperms;
	*ExponentsOut = Exp;
}

static double *basis_SI2_Bezier(const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut,
                                const unsigned int d)
{
	unsigned int n, Nbf;
	double       *ChiBez_rst, *BCoords;

	Nbf = (unsigned int) (factorial_ull(d+P)/(factorial_ull(d)*factorial_ull(P)));
	*NbfOut = Nbf;

	ChiBez_rst = malloc(Nn*Nbf * sizeof *ChiBez_rst); // keep

	BCoords = malloc(Nn*(d+1) * sizeof *BCoords); // free

	rst_to_barycentric_SI(Nn,d,rst,BCoords);
	for (n = 0; n < Nn; n++) {
		double b0, b1, b2;

		b0 = BCoords[0*Nn+n];
		b1 = BCoords[1*Nn+n];
		b2 = BCoords[2*Nn+n];

		if (P == 0) {
			ChiBez_rst[n*Nbf+0] = pow(b0,0.0)*pow(b1,0.0)*pow(b2,0.0);
		} else if (P == 1) {
			ChiBez_rst[n*Nbf+0] = pow(b0,1.0)*pow(b1,0.0)*pow(b2,0.0);
			ChiBez_rst[n*Nbf+1] = pow(b0,0.0)*pow(b1,1.0)*pow(b2,0.0);
			ChiBez_rst[n*Nbf+2] = pow(b0,0.0)*pow(b1,0.0)*pow(b2,1.0);
		} else if (P == 2) {
			ChiBez_rst[n*Nbf+0] =     pow(b0,2.0)*pow(b1,0.0)*pow(b2,0.0);
			ChiBez_rst[n*Nbf+1] =     pow(b0,0.0)*pow(b1,2.0)*pow(b2,0.0);
			ChiBez_rst[n*Nbf+2] =     pow(b0,0.0)*pow(b1,0.0)*pow(b2,2.0);
			ChiBez_rst[n*Nbf+3] = 2.0*pow(b0,1.0)*pow(b1,1.0)*pow(b2,0.0);
			ChiBez_rst[n*Nbf+4] = 2.0*pow(b0,0.0)*pow(b1,1.0)*pow(b2,1.0);
			ChiBez_rst[n*Nbf+5] = 2.0*pow(b0,1.0)*pow(b1,0.0)*pow(b2,1.0);
		} else if (0&&P == 3) {
			ChiBez_rst[n*Nbf+0] =     pow(b0,3.0)*pow(b1,0.0)*pow(b2,0.0);
			ChiBez_rst[n*Nbf+1] =     pow(b0,0.0)*pow(b1,3.0)*pow(b2,0.0);
			ChiBez_rst[n*Nbf+2] =     pow(b0,0.0)*pow(b1,0.0)*pow(b2,3.0);
			ChiBez_rst[n*Nbf+3] = 3.0*pow(b0,2.0)*pow(b1,1.0)*pow(b2,0.0);
			ChiBez_rst[n*Nbf+4] = 3.0*pow(b0,0.0)*pow(b1,2.0)*pow(b2,1.0);
			ChiBez_rst[n*Nbf+5] = 3.0*pow(b0,1.0)*pow(b1,0.0)*pow(b2,2.0);
			ChiBez_rst[n*Nbf+6] = 3.0*pow(b0,2.0)*pow(b1,0.0)*pow(b2,1.0);
			ChiBez_rst[n*Nbf+7] = 3.0*pow(b0,1.0)*pow(b1,2.0)*pow(b2,0.0);
			ChiBez_rst[n*Nbf+8] = 3.0*pow(b0,0.0)*pow(b1,1.0)*pow(b2,2.0);
			ChiBez_rst[n*Nbf+9] = 6.0*pow(b0,1.0)*pow(b1,1.0)*pow(b2,1.0);
		} else if (P == 4) {
			ChiBez_rst[n*Nbf+0 ] =      pow(b0,4.0)*pow(b1,0.0)*pow(b2,0.0);
			ChiBez_rst[n*Nbf+1 ] =      pow(b0,0.0)*pow(b1,4.0)*pow(b2,0.0);
			ChiBez_rst[n*Nbf+2 ] =      pow(b0,0.0)*pow(b1,0.0)*pow(b2,4.0);
			ChiBez_rst[n*Nbf+3 ] = 4.0* pow(b0,3.0)*pow(b1,1.0)*pow(b2,0.0);
			ChiBez_rst[n*Nbf+4 ] = 4.0* pow(b0,0.0)*pow(b1,3.0)*pow(b2,1.0);
			ChiBez_rst[n*Nbf+5 ] = 4.0* pow(b0,1.0)*pow(b1,0.0)*pow(b2,3.0);
			ChiBez_rst[n*Nbf+6 ] = 4.0* pow(b0,3.0)*pow(b1,0.0)*pow(b2,1.0);
			ChiBez_rst[n*Nbf+7 ] = 4.0* pow(b0,1.0)*pow(b1,3.0)*pow(b2,0.0);
			ChiBez_rst[n*Nbf+8 ] = 4.0* pow(b0,0.0)*pow(b1,1.0)*pow(b2,3.0);
			ChiBez_rst[n*Nbf+9 ] = 6.0* pow(b0,2.0)*pow(b1,2.0)*pow(b2,0.0);
			ChiBez_rst[n*Nbf+10] = 6.0* pow(b0,0.0)*pow(b1,2.0)*pow(b2,2.0);
			ChiBez_rst[n*Nbf+11] = 6.0* pow(b0,2.0)*pow(b1,0.0)*pow(b2,2.0);
			ChiBez_rst[n*Nbf+12] = 12.0*pow(b0,2.0)*pow(b1,1.0)*pow(b2,1.0);
			ChiBez_rst[n*Nbf+13] = 12.0*pow(b0,1.0)*pow(b1,2.0)*pow(b2,1.0);
			ChiBez_rst[n*Nbf+14] = 12.0*pow(b0,1.0)*pow(b1,1.0)*pow(b2,2.0);
		} else if (P < 9) {
			unsigned int permutationsTRI[18]  = { 0, 1, 2, 2, 0, 1, 1, 2, 0, 0, 2, 1, 1, 0, 2, 2, 1, 0};
			unsigned int NExp, *Nperms, *Exponents;

			get_BCoord_Exponents(P,d,&NExp,&Nperms,&Exponents); // free

			size_t IndB = Nbf-1;
			for (size_t e = 0; e < NExp; e++) {
				unsigned int *Exp_ptr = &Exponents[e*(d+1)];

				for (size_t perm = Nperms[e]; perm ; perm--) {
					unsigned int basis_exp[d+1], *perm_ptr;

					// Get exponents of each barycentric coordinate for the current basis function
					perm_ptr = &permutationsTRI[(perm-1)*(d+1)];
					for (size_t i = 0; i < d+1; i++)
						basis_exp[i] = Exp_ptr[perm_ptr[i]];

					// Compute scaling
					double Num, Den, Scaling;

					Num = factorial_d(P);
					Den = 1.0;
					for (size_t i = 0; i < d+1; i++)
						Den *= factorial_d(basis_exp[i]);

					Scaling = Num/Den;

					// Compute basis function
					for (size_t n = 0; n < Nn; n++)
						ChiBez_rst[n*Nbf+IndB] = Scaling;

					for (size_t i = 0; i < d+1; i++) {
						double *b = &BCoords[i*Nn];
						for (size_t n = 0; n < Nn; n++)
							ChiBez_rst[n*Nbf+IndB] *= pow(b[n],basis_exp[i]);
					}
					IndB--;
				}
			}

			free(Nperms);
			free(Exponents);
		} else {
			EXIT_UNSUPPORTED;
		}
	}

	free(BCoords);

	return ChiBez_rst;
}

static double *basis_SI3_Bezier(const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut,
                                const unsigned int d)
{
	if (0)
		printf("%d % .3e %d %d %d\n",P,rst[0],Nn,*NbfOut,d);
	return NULL;
}



/*
 *	Purpose:
 *		Test correctness of implementation of basis_SI_Bezier.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

static void test_unit_basis_SI_Bezier(void)
{
	unsigned int pass;

	char *PrintName = malloc(STRLEN_MAX * sizeof *PrintName); // free
	/*
	 *	barycentric_exponents_Bezier
	 *
	 *		Input:
	 *			P, dE.
	 *
	 *		Expected output:
	 *			Exponents, Number of permutations.
	 */

	for (unsigned int dE = 2; dE <= DMAX; dE++) {
		unsigned int PMin, PMax;

		if (dE == 2) {
			PMin = 0; PMax = 6; // First 6-symmetry occurs for P = 3
		} else if (dE == 3) {
			PMin = 0; PMax = 6; // First 24-symmetry occurs for P = 6
		} else {
			EXIT_UNSUPPORTED;
		}

		pass = 1;
		for (unsigned int P = PMin; P <= PMax; P++) {
			unsigned int NExp, *Nperms, *Nperms_test, *BCoord_Exponents, *BCoord_Exponents_test;

			get_BCoord_Exponents_test(P,dE,&NExp,&Nperms_test,&BCoord_Exponents_test); // free
			get_BCoord_Exponents(P,dE,&NExp,&Nperms,&BCoord_Exponents);                // free

			if (!(array_norm_diff_ui(NExp,Nperms,Nperms_test,"Inf") == 0 &&
			      array_norm_diff_ui(NExp*(dE+1),BCoord_Exponents,BCoord_Exponents_test,"Inf") == 0)) {
					pass = 0;

				printf("%d %d\n",dE,P);
				array_print_ui(1,NExp,Nperms,'R');
				array_print_ui(1,NExp,Nperms_test,'R');

				array_print_ui(NExp,dE+1,BCoord_Exponents,'R');
				array_print_ui(NExp,dE+1,BCoord_Exponents_test,'R');
			}

			free(Nperms);
			free(Nperms_test);
			free(BCoord_Exponents);
			free(BCoord_Exponents_test);
		}

		if (dE == 2)
			sprintf(PrintName,"BCoord_Exp_Bezier (d%d, P = [%d,%d]):",dE,PMin,PMax);
		else if (dE == 3)
			sprintf(PrintName,"                  (d%d, P = [%d,%d]):",dE,PMin,PMax);
		test_print2(pass,PrintName);
	}

	/*
	 *	basis_SI_Bezier
	 *
	 *		Input:
	 *
	 *			P, rst, Nn, dE.
	 *
	 *		Expected output:
	 *
	 *			ChiBez_rst = @(r,s)   [ See basis_SI2_Bezier ]
	 *			ChiBez_rst = @(r,s,t) [ See basis_SI3_Bezier ]
	 */

//	for (unsigned int dE = 2; dE <= DMAX; dE++) {
	for (unsigned int dE = 2; dE <= 2; dE++) {
		unsigned int PMin, PMax, Prst, Nn, Ns, Nbf, *symms;
		double       *rst, *w, *ChiBez_code, *ChiBez_test;

		cubature_tdef cubature;
		basis_tdef    basis_test;

		if (dE == 2) {
			PMin = 0; PMax = 8;
			basis_test = basis_SI2_Bezier;
			select_functions_cubature(&cubature,TRI);
		} else if (dE == 3) {
			PMin = 2; PMax = 4;
			basis_test = basis_SI3_Bezier;
			select_functions_cubature(&cubature,TET);
		} else {
			EXIT_UNSUPPORTED;
		}

		Prst = 4;
		cubature(&rst,&w,&symms,&Nn,&Ns,0,Prst,dE,"WV"); // free

		pass = 1;
		for (unsigned int P = PMin; P <= PMax; P++) {
			ChiBez_code = basis_SI_Bezier(P,rst,Nn,&Nbf,dE);
			ChiBez_test = basis_test(P,rst,Nn,&Nbf,dE);

			if (!(array_norm_diff_d(Nn*Nbf,ChiBez_code,ChiBez_test,"Inf") < EPS*10))
				pass = 0;

			free(ChiBez_test);
			free(ChiBez_code);
		}

		if (dE == 2) {
			sprintf(PrintName,"basis_SI_Bezier (d%d, P = [%d,%d]):",dE,PMin,PMax);
		} else {
			sprintf(PrintName,"                (d%d, P = [%d,%d]):",dE,PMin,PMax);
		}
		test_print2(pass,PrintName);

		free(rst);
		free(symms);
	}

	/*
	 *	basis_SI_Bezier partition of unity/positivity:
	 *
	 *		Input:
	 *
	 *			rst
	 *
	 *		Expected Output:
	 *
	 *			Basis functions at each node sum to one.
	 *			Basis functions evaluated at each node have values greater than or equal to zero.
	 */

//	for (unsigned int dE = 2; dE <= DMAX; dE++) {
	for (unsigned int dE = 2; dE <= 2; dE++) {
		unsigned int PMin, PMax, Prst, Nn, Ns, Nbf, *symms;
		double       *rst, *w, *ChiBez;

		cubature_tdef cubature;

		Prst = 4;

		if (dE == 2) {
			PMin = 0; PMax = 8;
			select_functions_cubature(&cubature,TRI);
		} else if (dE == 3) {
			PMin = 0; PMax = 8;
			select_functions_cubature(&cubature,TET);
		} else {
			EXIT_UNSUPPORTED;
		}

		cubature(&rst,&w,&symms,&Nn,&Ns,0,Prst,dE,"AO"); // free

		pass = 1;
		for (size_t P = PMin; P <= PMax; P++) {
			ChiBez = basis_SI_Bezier(P,rst,Nn,&Nbf,dE); // free

			double *RowSum, *ones;
			RowSum = calloc(Nn , sizeof *RowSum); // free
			ones   = calloc(Nn , sizeof *RowSum); // free
			for (size_t i = 0; i < Nn; i++) {
				ones[i] = 1.0;
				for (size_t j = 0; j < Nbf; j++) {
					RowSum[i] += ChiBez[j+i*Nbf];
				}
			}

			bool   Positive = 1;
			double minVal   = 1e15;
			for (size_t i = 0, iMax = Nn*Nbf; i < iMax; i++) {
				if (ChiBez[i] < minVal)
					minVal = ChiBez[i];
			}

			if (minVal < 0.0)
				Positive = 0;

			if (!(array_norm_diff_d(Nn,RowSum,ones,"Inf") < EPS*10 && Positive))
				pass = 0;

			free(ChiBez);
			free(RowSum);
			free(ones);
		}

		if (dE == 2) {
			sprintf(PrintName,"basis_SI_Bezier part. of unity (d%d, P = [%d,%d]):",dE,PMin,PMax);
		} else {
			sprintf(PrintName,"                               (d%d, P = [%d,%d]):",dE,PMin,PMax);
		}
		test_print2(pass,PrintName);

		free(rst);
		free(symms);
	}

	free(PrintName);
}
