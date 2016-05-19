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

void test_imp_grad_basis_PYR(void)
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

	unsigned int d, Nn, Ns, Nbf, P, Prst;
	unsigned int *symms;
	double *rst, *w;

	d = 3;

	double **GradChiRef32_code, **GradChiRef32_test;

	P = 2;
	Prst = 4;

	cubature_PYR(&rst,&w,&symms,&Nn,&Ns,0,Prst,d,"WV");  // free

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
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("grad_basis_PYR (P2):                             ");
	test_print(pass);

	free(rst);
	free(symms);
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

	cubature_PYR(&rst,&w,&symms,&Nn,&Ns,0,P,d,"GLL"); // free

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
			pass = 1, TestDB.Npass++;
	} else {
		printf("%e\n",array_norm_diff_d(Nn,f_r,f_rcomp,"Inf"));
		printf("%e\n",array_norm_diff_d(Nn,f_s,f_scomp,"Inf"));
		printf("%e\n",array_norm_diff_d(Nn,f_t,f_tcomp,"Inf"));
	}

	//     0         10        20        30        40        50
	printf("               derivatives (d3):                 ");
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
