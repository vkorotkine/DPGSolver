#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "test.h"
#include "functions.h"
#include "parameters.h"

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

static void poly2(const double *r, const double *s, const double *t, const unsigned int Nn, double **f, double **f_r,
                  double **f_s, double **f_t)
{
	unsigned int i;
	double *f_rst, *f_r_rst, *f_s_rst, *f_t_rst;
	double r_i, s_i, t_i;

	f_rst   = malloc(Nn * sizeof *f_rst);   // keep (requires external free)
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

void test_imp_grad_basis_SI(void)
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

	unsigned int d, Nn, Ns, Nbf, P, Prst;
	unsigned int *symms;
	double *rst, *w;

	d = 2;

	double **GradChiRef22_code, **GradChiRef22_test;

	P = 2;
	Prst = 4;

	cubature_TRI(&rst,&w,&symms,&Nn,&Ns,0,Prst,d,"AO"); // free

	GradChiRef22_code = grad_basis_SI(P,rst,Nn,&Nbf,d); // free
	GradChiRef22_test = grad_basis_SI22(rst,Nn);        // free

	pass = 0;
	if (array_norm_diff_d(Nn*Nbf,GradChiRef22_code[0],GradChiRef22_test[0],"Inf") < EPS*10 &&
	    array_norm_diff_d(Nn*Nbf,GradChiRef22_code[1],GradChiRef22_test[1],"Inf") < EPS*10)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("grad_basis_SI (d2, P2):                          ");
	test_print(pass);

	free(rst);
	free(symms);
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

	cubature_TET(&rst,&w,&symms,&Nn,&Ns,0,Prst,d,"AO"); // free

	GradChiRef31_code = grad_basis_SI(P,rst,Nn,&Nbf,d); // free
	GradChiRef31_test = grad_basis_SI31(rst,Nn); // free

	pass = 0;
	if (array_norm_diff_d(Nn*Nbf,GradChiRef31_code[0],GradChiRef31_test[0],"Inf") < EPS*10 &&
	    array_norm_diff_d(Nn*Nbf,GradChiRef31_code[1],GradChiRef31_test[1],"Inf") < EPS*10 &&
	    array_norm_diff_d(Nn*Nbf,GradChiRef31_code[2],GradChiRef31_test[2],"Inf") < EPS*10)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("              (d3, P1):                          ");
	test_print(pass);

	free(rst);
	free(symms);
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

	cubature_TRI(&rst,&w,&symms,&Nn,&Ns,0,P,d,"AO"); // free

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
			pass = 1, TestDB.Npass++;
	} else {
		printf("%e\n",array_norm_diff_d(Nn,f_r,f_rcomp,"Inf"));
		printf("%e\n",array_norm_diff_d(Nn,f_s,f_scomp,"Inf"));
	}

	//     0         10        20        30        40        50
	printf("              derivatives (d2):                  ");
	test_print(pass);

	free(rst), free(r), free(s), free(t);
	free(symms);
	free(I);
	free(ChiRef_rst), free(ChiRefInv_rst);
	array_free2_d(d,GradChiRef_rst);
	free(f), free(f_hat);
	free(f_r), free(f_s), free(f_t);
	free(f_rcomp), free(f_scomp), free(f_tcomp);

	// d = 3
	d = 3;

	cubature_TET(&rst,&w,&symms,&Nn,&Ns,0,P,d,"AO"); // free

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
	free(symms);
	free(I);
	free(ChiRef_rst), free(ChiRefInv_rst);
	array_free2_d(d,GradChiRef_rst);
	free(f), free(f_hat);
	free(f_r), free(f_s), free(f_t);
	free(f_rcomp), free(f_scomp), free(f_tcomp);
}
