// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "bases.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "mkl.h"

#include "Parameters.h"
#include "Macros.h"

#include "math_functions.h"

/*
 *	Purpose:
 *		Return the "matrix" ChiRef_rst representing the orthonormal basis functions evaluated at the provided quadrature
 *		nodes for polynomial order P.
 *
 *	Comments:
 *		The "matrix" is returned as a 1D array.
 *		Careful of overflow in gamma function in JacobiP for alpha > 20.0 (P ~ 20) for basis_SI.
 *
 *	Notation:
 *
 *	References:
 *		Hesthaven(2008)-Nodal_Discontinuous_Galerkin_Methods (Appendix A.1)
 *		Chan(2015)-Orthogonal_Bases_for_Vertex-Mapped_Pyramids
 *			Modal basis        : eq. 2.1
 *			Semi-nodal basis   : Following Lemma 2.2
 *		Witherden(2015)-On_the_development_and_Implementatin_of_High-Order_Flux_Reconstruction_Schemes_for_
 *		                Computational_Fluid_Dynamics
 *			pyfr (modal) basis : eq. 3.20
 */

double jacobiP      (const double x, const double alpha, const double beta, const int N);
double grad_jacobiP (const double x, const double alpha, const double beta, const int N);

void rst_to_abc_SI  (const unsigned int Nn, const unsigned int d, const double *rst, double *a, double *b, double *c);
void rst_to_abc_PYR (const unsigned int Nn, const unsigned int d, const double *rst, double *a, double *b, double *c);

double *basis_TP(const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut,
                 const unsigned int d)
{
	unsigned int i, j, k, iMax, jMax, kMax, Indbf, Indn, u1,
	             N, Nbf;
	int          sd, sN;
	double       *ChiRef_rst, *r, *s, *t;

	// silence
	s = t = NULL;

	N = P+1;

	u1 = 1;
	sd = d;
	sN = N;

	r = malloc(Nn * sizeof *r); // free
	if (d > 1) s = malloc(Nn * sizeof *s); // free
	if (d > 2) t = malloc(Nn * sizeof *t); // free

	for (i = 0; i < Nn; i++) {
		r[i] = rst[0*Nn+i];
		if (d > 1) s[i] = rst[1*Nn+i];
		if (d > 2) t[i] = rst[2*Nn+i];
	}

	Nbf = pow(N,d);

	ChiRef_rst = malloc(Nn*Nbf * sizeof *ChiRef_rst); // keep (requires external free)

	Indbf = 0;
	for (k = 0, kMax = (unsigned int) min(max((sd-2)*sN,1),sN); k < kMax; k++) {
	for (j = 0, jMax = min(max((d-1)*N,u1),N); j < jMax; j++) {
	for (i = 0, iMax = min(max((d)  *N,u1),N); i < iMax; i++) {
		for (Indn = 0; Indn < Nn; Indn++) {
			           ChiRef_rst[Indbf*Nn+Indn] =  jacobiP(r[Indn],0.0,0.0,(double) i);
			if (d > 1) ChiRef_rst[Indbf*Nn+Indn] *= jacobiP(s[Indn],0.0,0.0,(double) j);
			if (d > 2) ChiRef_rst[Indbf*Nn+Indn] *= jacobiP(t[Indn],0.0,0.0,(double) k);
		}
		Indbf++;
	}}}

	// Transpose ChiRef_rst
	mkl_dimatcopy('R','T',Nbf,Nn,1.,ChiRef_rst,Nn,Nbf);

// array_print_d(Nn,Nbf,ChiRef_rst,'R');

	free(r);
	if (d > 1) free(s);
	if (d > 2) free(t);

	*NbfOut = Nbf;
	return ChiRef_rst;
}

double *basis_SI(const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut,
                 const unsigned int d)
{
	unsigned int i, j, k, iMax, jMax, kMax, Indbf, Indn, Nbf;
	double       *ChiRef_rst, *a, *b, *c;

	if (d < 2 || d > 3)
		printf("Error: basis_SI only supports d = [2,3].\n"), EXIT_MSG;

	// Convert from rst to abc coordinates
	a = malloc(Nn * sizeof *a); // free
	b = malloc(Nn * sizeof *b); // free
	c = malloc(Nn * sizeof *c); // free

	rst_to_abc_SI(Nn,d,rst,a,b,c);

	Nbf = (unsigned int) (factorial_ull(d+P)/(factorial_ull(d)*factorial_ull(P)));

	ChiRef_rst = malloc(Nn*Nbf * sizeof *ChiRef_rst); // keep (requires external free)

	Indbf = 0;
	for (i = 0, iMax = P;   i <= iMax; i++) {
	for (j = 0, jMax = P-i; j <= jMax; j++) {
	for (k = 0, kMax = (d-2)*(P-i-j); k <= kMax; k++) {
		for (Indn = 0; Indn < Nn; Indn++) {
			if (d == 2)
				ChiRef_rst[Indbf*Nn+Indn] = 2.0/pow(3.0,0.25)*pow(1.0-b[Indn],i);
			else
				ChiRef_rst[Indbf*Nn+Indn] = 4.0/pow(2.0,0.25)*pow(1.0-b[Indn],i)*pow(1.0-c[Indn],i+j);

			ChiRef_rst[Indbf*Nn+Indn] *= jacobiP(a[Indn],0.0,0.0,(double) i);
			ChiRef_rst[Indbf*Nn+Indn] *= jacobiP(b[Indn],2.0*i+1.0,0.0,(double) j);
			if (d == 3)
				ChiRef_rst[Indbf*Nn+Indn] *= jacobiP(c[Indn],2.0*(i+j+1.0),0.0,(double) k);
		}
		Indbf++;
	}}}

	// Transpose ChiRef_rst
	mkl_dimatcopy('R','T',Nbf,Nn,1.,ChiRef_rst,Nn,Nbf);

//	array_print_d(Nn,Nbf,ChiRef_rst,'R');

	free(a);
	free(b);
	free(c);

	*NbfOut = Nbf;
	return ChiRef_rst;
}

double *basis_PYR(const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut,
                  const unsigned int d)
{
/*
 *	Comments:
 *		The modal PYR basis can be transformed to Chan's semi-nodal basis using a transformation matrix as the two bases
 *		span the same space.
 */

	unsigned int i, j, k, iMax, jMax, kMax, Indbf, Indn, Nbf, mu_ij, IndChi;
	double       *ChiRef_rst, *a, *b, *c, a_n, b_n, c_n;

	if (d != 3)
		printf("Error: basis_PYR only supports d = 3.\n"), EXIT_MSG;

	// Convert from rst to abc coordinates
	a = malloc(Nn * sizeof *a); // free
	b = malloc(Nn * sizeof *b); // free
	c = malloc(Nn * sizeof *c); // free

	rst_to_abc_PYR(Nn,d,rst,a,b,c);

	Nbf = (unsigned int) (1.0/6.0*((P+1)*(P+2)*(2*P+3)));

	ChiRef_rst = malloc(Nn*Nbf * sizeof *ChiRef_rst); // keep (requires external free)

	Indbf = 0;
	for (i = 0, iMax = P; i <= iMax; i++) {
	for (j = 0, jMax = P; j <= jMax; j++) {
		mu_ij = max(i,j);
		for (k = 0, kMax = P-mu_ij; k <= kMax; k++) {
//			printf("%d %d %d %d %d %d\n",Indbf,i,j,k,i+j+k,mu_ij+k);
			for (Indn = 0; Indn < Nn; Indn++) {
				IndChi = Indbf*Nn+Indn;
				a_n = a[Indn];
				b_n = b[Indn];
				c_n = c[Indn];

// Orthonormal (GL nodes as expected, WV nodes of much higher order than expected)
				ChiRef_rst[IndChi]  = pow(2.0,1.25)*pow(1.0-c_n,(double) mu_ij);
				ChiRef_rst[IndChi] *= jacobiP(c_n,2.0*(mu_ij+1.0),0.0,(double) k);

// Orthonormal (GL nodes as expected, WV nodes of slightly higher order than expected)
				// Basis from pyfr code (ToBeDeleted)
//				ChiRef_rst[IndChi]  = pow(2.0,1.25)*pow(1.0-c_n,(double) (i+j));
//				ChiRef_rst[IndChi] *= jacobiP(c_n,2.0*(i+j+1.0),0.0,(double) k);

// Remaining part of the basis
				// Uncomment to test integration capability of variation only in a anb b
//				ChiRef_rst[IndChi]  = pow(2.0,-0.25)*sqrt(3.0);
				ChiRef_rst[IndChi] *= jacobiP(a_n,0.0,0.0,(double) i);
				ChiRef_rst[IndChi] *= jacobiP(b_n,0.0,0.0,(double) j);
		}
		Indbf++;
	}}}

	// Transpose ChiRef_rst
	mkl_dimatcopy('R','T',Nbf,Nn,1.,ChiRef_rst,Nn,Nbf);

//	array_print_d(Nn,Nbf,ChiRef_rst,'R');

	free(a);
	free(b);
	free(c);

	*NbfOut = Nbf;
	return ChiRef_rst;
}


/*
 *	Purpose:
 *		Return the "matrices" GradChiRef_rst representing the gradients of the orthonormal basis functions evaluated at
 *		the provided quadrature nodes for polynomial order P.
 *
 *	Comments:
 *		The "matrices" are returned as d pointers to 1D arrays.
 *
 *	Notation:
 *
 *	References:
 *		ToBeModified: Add reference to thesis where it is explained how these derivatives are computed.
 */

double **grad_basis_TP(const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut,
                       const unsigned int d)
{
	unsigned int i, j, k, iMax, jMax, kMax, dim, Indbf, Indn, u1,
	             N, Nbf;
	int          sd,sN;
	double       **GradChiRef_rst, *r, *s, *t;

	N = P+1;

	u1 = 1;
	sd = d;
	sN = N;

	r = malloc(Nn * sizeof *r); // free
	s = NULL; // silence
	t = NULL; // silence
	if (d > 1) s = malloc(Nn * sizeof *s); // free
	if (d > 2) t = malloc(Nn * sizeof *t); // free

	for (i = 0; i < Nn; i++) {
		r[i] = rst[0*Nn+i];
		if (d > 1) s[i] = rst[1*Nn+i];
		if (d > 2) t[i] = rst[2*Nn+i];
	}

	Nbf = pow(N,d);

	GradChiRef_rst = malloc(d * sizeof *GradChiRef_rst); // keep (requires external free)
	for (dim = 0; dim < d; dim++)
		GradChiRef_rst[dim] = malloc(Nn*Nbf * sizeof **GradChiRef_rst); // keep (requires external free)

	Indbf = 0;
	if (d == 1) {
		for (i = 0, iMax = min(max((d)  *N,u1),N); i < iMax; i++) {
			for (Indn = 0; Indn < Nn; Indn++) {
				GradChiRef_rst[0][Indbf*Nn+Indn] = grad_jacobiP(r[Indn],0.0,0.0,(double) i);
			}
			Indbf++;
		}
	} else if (d == 2) {
		for (j = 0, jMax = min(max((d-1)*N,u1),N); j < jMax; j++) {
		for (i = 0, iMax = min(max((d)  *N,u1),N); i < iMax; i++) {
			for (Indn = 0; Indn < Nn; Indn++) {
				GradChiRef_rst[0][Indbf*Nn+Indn] = grad_jacobiP(r[Indn],0.0,0.0,(double) i)*
				                                        jacobiP(s[Indn],0.0,0.0,(double) j);
				GradChiRef_rst[1][Indbf*Nn+Indn] =      jacobiP(r[Indn],0.0,0.0,(double) i)*
				                                   grad_jacobiP(s[Indn],0.0,0.0,(double) j);
			}
			Indbf++;
		}}
	} else if (d == 3) {
		for (k = 0, kMax = (unsigned int) min(max((sd-2)*sN,1),sN); k < kMax; k++) {
		for (j = 0, jMax = min(max((d-1)*N,u1),N); j < jMax; j++) {
		for (i = 0, iMax = min(max((d)  *N,u1),N); i < iMax; i++) {
			for (Indn = 0; Indn < Nn; Indn++) {
				GradChiRef_rst[0][Indbf*Nn+Indn] = grad_jacobiP(r[Indn],0.0,0.0,(double) i)*
				                                        jacobiP(s[Indn],0.0,0.0,(double) j)*
				                                        jacobiP(t[Indn],0.0,0.0,(double) k);
				GradChiRef_rst[1][Indbf*Nn+Indn] =      jacobiP(r[Indn],0.0,0.0,(double) i)*
				                                   grad_jacobiP(s[Indn],0.0,0.0,(double) j)*
				                                        jacobiP(t[Indn],0.0,0.0,(double) k);
				GradChiRef_rst[2][Indbf*Nn+Indn] =      jacobiP(r[Indn],0.0,0.0,(double) i)*
				                                        jacobiP(s[Indn],0.0,0.0,(double) j)*
				                                   grad_jacobiP(t[Indn],0.0,0.0,(double) k);
			}
			Indbf++;
		}}}
	}

	// Transpose GradChiRef_rst
	for (dim = 0; dim < d; dim++)
		mkl_dimatcopy('R','T',Nbf,Nn,1.,GradChiRef_rst[dim],Nn,Nbf);

// array_print_d(Nn,Nbf,GradChiRef_rst[0],'R');

	free(r);
	if (d > 1) free(s);
	if (d > 2) free(t);

	*NbfOut = Nbf;
	return GradChiRef_rst;
}

double **grad_basis_SI(const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut,
                       const unsigned int d)
{
	unsigned int i, j, k, iMax, jMax, kMax, dim, pder, Indbf, Indn, IndGrad, Nbf;
	double       **GradChiRef_rst, *a, *b, *c, a_n, b_n, c_n, jPa, jPb, jPc, djPa, djPb, djPc;
	double       GradChiRef_tmp[3][5];

	if (d < 2 || d > 3)
		printf("Error: grad_basis_SI only supports d = [2,3].\n"), EXIT_MSG;

	// Convert from rst to abc coordinates
	a = malloc(Nn * sizeof *a); // free
	b = malloc(Nn * sizeof *b); // free
	c = malloc(Nn * sizeof *c); // free

	rst_to_abc_SI(Nn,d,rst,a,b,c);

	Nbf = (unsigned int) (factorial_ull(d+P)/(factorial_ull(d)*factorial_ull(P)));

	GradChiRef_rst = malloc(d * sizeof *GradChiRef_rst); // keep (requires external free)
	for (dim = 0; dim < d; dim++)
		 GradChiRef_rst[dim] = calloc(Nn*Nbf , sizeof **GradChiRef_rst); // keep (requires external free)

	Indbf = 0;
	for (i = 0, iMax = P;   i <= iMax; i++) {
	for (j = 0, jMax = P-i; j <= jMax; j++) {
	for (k = 0, kMax = (d-2)*(P-i-j); k <= kMax; k++) {
		for (Indn = 0; Indn < Nn; Indn++) {
			IndGrad = Indbf*Nn+Indn;
			a_n = a[Indn];
			b_n = b[Indn];
			c_n = c[Indn];

			jPa = jacobiP(a_n,0.0,          0.0,(double) i);
			jPb = jacobiP(b_n,2.0*i+1.0,    0.0,(double) j);
			jPc = jacobiP(c_n,2.0*(i+j+1.0),0.0,(double) k);

			djPa = grad_jacobiP(a_n,0.0,          0.0,(double) i);
			djPb = grad_jacobiP(b_n,2.0*i+1.0,    0.0,(double) j);
			djPc = grad_jacobiP(c_n,2.0*(i+j+1.0),0.0,(double) k);

			for (dim = 0; dim < d; dim++)
				for (pder = 0; pder < 5; pder++ )
					GradChiRef_tmp[dim][pder] = 0.0;

			// Obtain contributions from each partial derivative
			GradChiRef_tmp[0][0] = djPa* jPb;
			GradChiRef_tmp[1][0] = djPa* jPb;
			GradChiRef_tmp[1][1] =  jPa*djPb;
			GradChiRef_tmp[1][2] =  jPa* jPb;

			if (d == 3) {
				// djPa = 0 when i = 0
				// djPb = 0 when j = 0
				// djPc = 0 when k = 0
				GradChiRef_tmp[0][0] *= jPc;
				GradChiRef_tmp[1][0] *= jPc;
				GradChiRef_tmp[1][1] *= jPc;
				GradChiRef_tmp[1][2] *= jPc;
				GradChiRef_tmp[2][0]  = djPa* jPb* jPc;
				GradChiRef_tmp[2][1]  =  jPa*djPb* jPc;
				GradChiRef_tmp[2][2]  =  jPa* jPb* jPc;
				GradChiRef_tmp[2][3]  =  jPa* jPb*djPc;
				GradChiRef_tmp[2][4]  =  jPa* jPb* jPc;
			}

			GradChiRef_tmp[0][0] *=  2.0                  ;
			GradChiRef_tmp[1][0] *=  2.0/3.0*sqrt(3.0)*a_n;
			GradChiRef_tmp[1][1] *=  2.0/3.0*sqrt(3.0)    ;
			GradChiRef_tmp[1][2] *= -2.0/3.0*sqrt(3.0)*i  ;
			if (i > 0) {
				GradChiRef_tmp[0][0] *= pow(1.0-b_n,i-1.0);
				GradChiRef_tmp[1][0] *= pow(1.0-b_n,i-1.0);
				GradChiRef_tmp[1][2] *= pow(1.0-b_n,i-1.0);
			} else {
				GradChiRef_tmp[0][0] *= 0.0; // redundant (i = 0 -> djPa = 0)
				GradChiRef_tmp[1][0] *= 0.0; // redundant (i = 0 -> djPa = 0)
				GradChiRef_tmp[1][2] *= 0.0;
			}
			GradChiRef_tmp[1][1] *= pow(1.0-b_n,(double) i);

			if (d == 3) {
				GradChiRef_tmp[0][0] *=  2.0;
				GradChiRef_tmp[1][0] *=  2.0;
				GradChiRef_tmp[1][1] *=  2.0;
				GradChiRef_tmp[1][2] *=  2.0;
				GradChiRef_tmp[2][0] *=  2.0/3.0*sqrt(6.0)*a_n            ;
				GradChiRef_tmp[2][1] *=  1.0/6.0*sqrt(6.0)*(3.0*b_n+1.0)  ;
				GradChiRef_tmp[2][2] *= -1.0/6.0*sqrt(6.0)*(3.0*b_n+1.0)*i;
				GradChiRef_tmp[2][3] *=  1.0/2.0*sqrt(6.0)                ;
				GradChiRef_tmp[2][4] *= -1.0/2.0*sqrt(6.0)*(i+j)          ;
				if (i > 0) {
					GradChiRef_tmp[2][0] *= pow(1.0-b_n,i-1.0);
					GradChiRef_tmp[2][2] *= pow(1.0-b_n,i-1.0);
				} else {
					GradChiRef_tmp[2][0] *= 0.0;
					GradChiRef_tmp[2][2] *= 0.0;
				}
				GradChiRef_tmp[2][1] *= pow(1.0-b_n,(double) i);
				GradChiRef_tmp[2][3] *= pow(1.0-b_n,(double) i);
				GradChiRef_tmp[2][4] *= pow(1.0-b_n,(double) i);

				if (i+j > 0) {
					GradChiRef_tmp[0][0] *= pow(1.0-c_n,i+j-1.0);
					GradChiRef_tmp[1][0] *= pow(1.0-c_n,i+j-1.0);
					GradChiRef_tmp[1][1] *= pow(1.0-c_n,i+j-1.0);
					GradChiRef_tmp[1][2] *= pow(1.0-c_n,i+j-1.0);
					GradChiRef_tmp[2][0] *= pow(1.0-c_n,i+j-1.0);
					GradChiRef_tmp[2][1] *= pow(1.0-c_n,i+j-1.0);
					GradChiRef_tmp[2][2] *= pow(1.0-c_n,i+j-1.0);
					GradChiRef_tmp[2][4] *= pow(1.0-c_n,i+j-1.0);
				} else {
					GradChiRef_tmp[0][0] *= 0.0;
					GradChiRef_tmp[1][0] *= 0.0;
					GradChiRef_tmp[1][1] *= 0.0;
					GradChiRef_tmp[1][2] *= 0.0;
					GradChiRef_tmp[2][0] *= 0.0;
					GradChiRef_tmp[2][1] *= 0.0;
					GradChiRef_tmp[2][2] *= 0.0;
					GradChiRef_tmp[2][4] *= 0.0;
				}
				GradChiRef_tmp[2][3] *= pow(1.0-c_n,(double) (i+j));
			}

			// Sum contributions from all partial derivatives
			for (dim = 0; dim < d; dim++)
				for (pder = 0; pder < 5; pder++)
					GradChiRef_rst[dim][IndGrad] += GradChiRef_tmp[dim][pder];

			// Add scaling constant
			if (d == 2) {
				for (dim = 0; dim < d; dim++)
					GradChiRef_rst[dim][IndGrad] *= 2.0/pow(3.0,0.25);
			} else {
				for (dim = 0; dim < d; dim++)
					GradChiRef_rst[dim][IndGrad] *= 4.0/pow(2.0,0.25);
			}
		}
		Indbf++;
	}}}

	// Transpose GradChiRef_rst
	for (dim = 0; dim < d; dim++)
		mkl_dimatcopy('R','T',Nbf,Nn,1.,GradChiRef_rst[dim],Nn,Nbf);

//	array_print_d(Nn,Nbf,GradChiRef_rst[0],'R');

	free(a);
	free(b);
	free(c);

	*NbfOut = Nbf;
	return GradChiRef_rst;
}

double **grad_basis_PYR(const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut,
                        const unsigned int d)
{
	unsigned int i, j, k, iMax, jMax, kMax, dim, pder, Indbf, Indn, IndGrad, Nbf, mu_ij;
	double       **GradChiRef_rst, *a, *b, *c, a_n, b_n, c_n, jPa, jPb, jPc, djPa, djPb, djPc;
	double       GradChiRef_tmp[3][5];

	if (d != 3)
		printf("Error: grad_basis_PYR only supports d = 3.\n"), EXIT_MSG;

	// Convert from rst to abc coordinates
	a = malloc(Nn * sizeof *a); // free
	b = malloc(Nn * sizeof *b); // free
	c = malloc(Nn * sizeof *c); // free

	rst_to_abc_PYR(Nn,d,rst,a,b,c);

	Nbf = (unsigned int) (1.0/6.0*((P+1)*(P+2)*(2*P+3)));

	GradChiRef_rst = malloc(d * sizeof *GradChiRef_rst); // keep (requires external free)
	for (dim = 0; dim < d; dim++)
		 GradChiRef_rst[dim] = calloc(Nn*Nbf , sizeof **GradChiRef_rst); // keep (requires external free)

	Indbf = 0;
	for (i = 0, iMax = P; i <= iMax; i++) {
	for (j = 0, jMax = P; j <= jMax; j++) {
		mu_ij = max(i,j);
		for (k = 0, kMax = P-mu_ij; k <= kMax; k++) {
			for (Indn = 0; Indn < Nn; Indn++) {
				IndGrad = Indbf*Nn+Indn;
				a_n = a[Indn];
				b_n = b[Indn];
				c_n = c[Indn];

				jPa = jacobiP(a_n,0.0,            0.0,(double) i);
				jPb = jacobiP(b_n,0.0,            0.0,(double) j);
				jPc = jacobiP(c_n,2.0*(mu_ij+1.0),0.0,(double) k);

				djPa = grad_jacobiP(a_n,0.0,            0.0,(double) i);
				djPb = grad_jacobiP(b_n,0.0,            0.0,(double) j);
				djPc = grad_jacobiP(c_n,2.0*(mu_ij+1.0),0.0,(double) k);

				for (dim = 0; dim < d; dim++)
					for (pder = 0; pder < 4; pder++ )
						GradChiRef_tmp[dim][pder] = 0.0;

				// Obtain contributions from each partial derivative
				GradChiRef_tmp[0][0] = djPa* jPb* jPc;
				GradChiRef_tmp[1][1] =  jPa*djPb* jPc;
				GradChiRef_tmp[2][0] = djPa* jPb* jPc;
				GradChiRef_tmp[2][1] =  jPa*djPb* jPc;
				GradChiRef_tmp[2][2] =  jPa* jPb*djPc;
				GradChiRef_tmp[2][3] =  jPa* jPb* jPc;

				GradChiRef_tmp[0][0] *=  2.0          ;
				GradChiRef_tmp[1][1] *=  2.0          ;
				GradChiRef_tmp[2][0] *=  sqrt(2.0)*a_n;
				GradChiRef_tmp[2][1] *=  sqrt(2.0)*b_n;
				GradChiRef_tmp[2][2] *=  sqrt(2.0)    ;
				GradChiRef_tmp[2][3] *= -sqrt(2.0)*mu_ij    ;

				if (mu_ij > 0) {
					GradChiRef_tmp[0][0] *= pow(1.0-c_n,mu_ij-1.0);
					GradChiRef_tmp[1][1] *= pow(1.0-c_n,mu_ij-1.0);
					GradChiRef_tmp[2][0] *= pow(1.0-c_n,mu_ij-1.0);
					GradChiRef_tmp[2][1] *= pow(1.0-c_n,mu_ij-1.0);
					GradChiRef_tmp[2][3] *= pow(1.0-c_n,mu_ij-1.0);
				} else {
					GradChiRef_tmp[0][0] *= 0.0;
					GradChiRef_tmp[1][1] *= 0.0;
					GradChiRef_tmp[2][0] *= 0.0;
					GradChiRef_tmp[2][1] *= 0.0;
					GradChiRef_tmp[2][3] *= 0.0;
				}
				GradChiRef_tmp[2][2] *= pow(1.0-c_n,(double) mu_ij);

				// Sum contributions from all partial derivatives
				for (dim = 0; dim < d; dim++)
					for (pder = 0; pder < 4; pder++)
						GradChiRef_rst[dim][IndGrad] += GradChiRef_tmp[dim][pder];

				// Add scaling constant
				for (dim = 0; dim < d; dim++)
					GradChiRef_rst[dim][IndGrad] *= pow(2.0,1.25);
			}
			Indbf++;
		}
	}}

	// Transpose GradChiRef_rst
	for (dim = 0; dim < d; dim++)
		mkl_dimatcopy('R','T',Nbf,Nn,1.,GradChiRef_rst[dim],Nn,Nbf);

//	array_print_d(Nn,Nbf,GradChiRef_rst[0],'R');

	free(a);
	free(b);
	free(c);

	*NbfOut = Nbf;
	return GradChiRef_rst;
}


/*
 *	Purpose:
 *		Evaluate Jacobi Polynomial of type (alpha,beta) > -1, (alpha+beta <> -1) at the point x for order N.
 *
 *	Comments:
 *		Polynomials are normalized to be orthonormal.
 *
 *	Notation:
 *
 *	References:
 *		Hesthaven (Nodal DG Code): https://github.com/tcew/nodal-dg
 */

double jacobiP(const double x, const double alpha, const double beta, const int N)
{
	int    i, iMax;
	double aold = 0.0, anew = 0.0, bnew = 0.0, h1 = 0.0,
	       gamma0 = 0.0, gamma1 = 0.0,
	       ab = alpha + beta, ab1 = ab + 1.0, a1 = alpha + 1.0, b1 = beta + 1.0,
	       x_bnew, P, *PL;

	PL = malloc((N+1) * sizeof *PL); // free

	// Initial values P_0(x) and P_1(x)
	gamma0 = pow(2,ab1)/(ab1)*gamma_d(a1)*gamma_d(b1)/gamma_d(ab1);

	PL[0] = 1.0/sqrt(gamma0);

	if (N == 0) {
		P = PL[0];
		free(PL);
		return P;
	}

	gamma1 = a1*b1/(ab+3.0)*gamma0;
	PL[1] = (((ab+2.0)*x + (alpha-beta))/2.0) / sqrt(gamma1);

	if (N == 1) {
		P = PL[1];
		free(PL);
		return P;
	}

	// Repeat value in recurrence.
	aold = 2.0/(2.0+ab)*sqrt(a1*b1/(ab+3.0));

	// Forward recurrence using the symmetry of the recurrence.
	for (i = 1, iMax = N-1; i <= iMax; i++) {
		h1 = 2.0*i+ab;
		anew = 2.0/(h1+2.0)*sqrt((i+1)*(i+ab1)*(i+a1)*(i+b1)/(h1+1.0)/(h1+3.0));
		bnew = -(pow(alpha,2)-pow(beta,2))/(h1*(h1+2.0));
		x_bnew = x-bnew;
		PL[i+1] = 1.0/anew * (-aold*PL[i-1] + x_bnew*PL[i]);
		aold = anew;
	}

	P = PL[N];
	free(PL);

	return P;
}

double grad_jacobiP(const double x, const double alpha, const double beta, const int N)
{
	double dP;

	if (N == 0)
		dP = 0.0;
	else
		dP = sqrt(N*(N+alpha+beta+1))*jacobiP(x,alpha+1.0,beta+1.0,N-1);

	return dP;
}


/*
 *	Purpose:
 *		Convert from rst to abc coordinates using Duffy-type transforms.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void rst_to_abc_SI(const unsigned int Nn, const unsigned int d, const double *rst, double *a, double *b, double *c)
{
	unsigned int n;
	double r_n, s_n, t_n;
	double *r, *s, *t;

	if (d < 2)
		printf("Error: Conversion from rst to abc SI coordinates should only be used for d >= 2.\n"), exit(1);

	r = malloc(Nn * sizeof *r); // free
	s = malloc(Nn * sizeof *s); // free
	t = malloc(Nn * sizeof *t); // free

	for (n = 0; n < Nn; n++) {
		r[n] = rst[0*Nn+n];
		s[n] = rst[1*Nn+n];
		if (d > 2)
			t[n] = rst[2*Nn+n];
		else
			t[n] = -1.0/sqrt(6.0);
	}

	for (n = 0; n < Nn; n++) {
		r_n = r[n];
		s_n = s[n];
		t_n = t[n];

		if (fabs(2*sqrt(3.0)*s_n+sqrt(6.0)*t_n-3.0) > 1e2*EPS)
			a[n] = 6.0*r_n/(3.0-2.0*sqrt(3.0)*s_n-sqrt(6.0)*t_n);
		else // On top line of the regular TET / At the top of the regular TRI
			a[n] = 0.0;

		if (fabs(sqrt(6.0)*t_n-3.0) > 1e2*EPS)
			b[n] = 1.0/3.0*(8.0*sqrt(3.0)*s_n/(3.0-sqrt(6.0)*t_n)-1.0);
		else // At the top of the regular TET
			b[n] = 0.0;

		c[n] = 0.5*(sqrt(6.0)*t_n-1.0);
	}

	free(r);
	free(s);
	free(t);
}

void rst_to_abc_PYR(const unsigned int Nn, const unsigned int d, const double *rst, double *a, double *b, double *c)
{
	unsigned int n;
	double r_n, s_n, t_n;
	double *r, *s, *t;

	if (d != 3)
		printf("Error: Conversion from rst to abc PYR coordinates should only be used for d = 3.\n"), exit(1);

	r = malloc(Nn * sizeof *r); // free
	s = malloc(Nn * sizeof *s); // free
	t = malloc(Nn * sizeof *t); // free

	for (n = 0; n < Nn; n++) {
		r[n] = rst[0*Nn+n];
		s[n] = rst[1*Nn+n];
		t[n] = rst[2*Nn+n];
	}

	for (n = 0; n < Nn; n++) {
		r_n = r[n];
		s_n = s[n];
		t_n = t[n];

		if (fabs(0.8*sqrt(2.0)-t_n) > 100*EPS) {
			a[n] = r_n/(0.8-1.0/sqrt(2.0)*t_n);
			b[n] = s_n/(0.8-1.0/sqrt(2.0)*t_n);
		} else { // At the top of the pyramid
			a[n] = 0.0;
			b[n] = 0.0;
		}

		c[n] = sqrt(2.0)*t_n-0.6;
	}

	free(r);
	free(s);
	free(t);
}
