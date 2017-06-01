// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "bases.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <stdbool.h>

#include "mkl.h"

#include "Parameters.h"
#include "Macros.h"
#include "S_ELEMENT.h"

#include "math_functions.h"
#include "matrix_functions.h"
#include "element_functions.h"
#include "setup_operators_support.h"
#include "array_sort.h"
#include "array_norm.h"
#include "array_print.h"

#include "matrix_structs.h"

/*
 *	Purpose:
 *		Return the "matrix" ChiRef_rst representing the orthonormal or Bezier basis functions evaluated at the provided
 *		quadrature nodes for polynomial order P.
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
void rst_to_barycentric_SI (const unsigned int Nn, const unsigned int d, const double *rst, double *BCoords);

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

double *basis_TP_Bezier(const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut,
                        const unsigned int d)
{
	/*
	 *	Comments:
	 *		The recursive definition of the basis functions is used in the computation below.
	 *		The polynomials are defined on the 1D reference element where r \in [-1,1].
	 *
	 *	References:
	 *		Prautzsch(2002)-Bezier_and_B-Spline_Techniques (Ch. 2.1)
	 */

	unsigned int N, Nbf;
	double       *ChiBez_rst;

	N   = P+1;
	Nbf = pow(N,d);

	ChiBez_rst = calloc(Nn*Nbf , sizeof *ChiBez_rst); // keep

	// r-direction
	const double *r_ptr = &rst[0*Nn];
	for (size_t i = 0, iMax = N; i < iMax; i++) {
	for (size_t Indbf = i; Indbf+1 ; Indbf--) {
		for (size_t n = 0; n < Nn; n++) {
			double r = r_ptr[n];

			if (i == 0) {
				ChiBez_rst[Indbf*Nn+n] = 1.0;
			} else {
				double b0, b1;
				b0 = (1.0-r)/2.0;
				b1 = (1.0+r)/2.0;

				if (Indbf != i)
					ChiBez_rst[Indbf*Nn+n] *= b0;
				if (Indbf != 0)
					ChiBez_rst[Indbf*Nn+n] += b1*ChiBez_rst[(Indbf-1)*Nn+n];
			}
		}
	}}

	// s-direction
	unsigned int u1 = 1;
	if (d > 1) {
		const double *s_ptr = &rst[1*Nn];
		// j = 0 gives multiplication by 1.0
		for (size_t j = 1, jMax = min(max((d-1)*N,u1),N); j < jMax; j++) {
		for (size_t Indbf = j; Indbf+1 ; Indbf--) {
		for (size_t i = 0, iMax = N; i < iMax; i++) {
			for (size_t n = 0; n < Nn; n++) {
				double s = s_ptr[n];

				double b0, b1;
				b0 = (1.0-s)/2.0;
				b1 = (1.0+s)/2.0;

				if (Indbf != j)
					ChiBez_rst[(Indbf*N+i)*Nn+n] *= b0;
				if (Indbf != 0)
					ChiBez_rst[(Indbf*N+i)*Nn+n] += b1*ChiBez_rst[((Indbf-1)*N+i)*Nn+n];
			}
		}}}
	}

	// t-direction
	if (d > 2) {
		int          sd = d, sN = N;
		const double *t_ptr = &rst[2*Nn];
		// k = 0 gives multiplication by 1.0
		for (size_t k = 1, kMax = (size_t) min(max((sd-2)*sN,1),sN); k < kMax; k++) {
		for (size_t Indbf = k; Indbf+1 ; Indbf--) {
		for (size_t j = 0, jMax = min(max((d-1)*N,u1),N); j < jMax; j++) {
		for (size_t i = 0, iMax = N; i < iMax; i++) {
			for (size_t n = 0; n < Nn; n++) {
				double t = t_ptr[n];

				double b0, b1;
				b0 = (1.0-t)/2.0;
				b1 = (1.0+t)/2.0;
				if (Indbf != k)
					ChiBez_rst[(Indbf*N*N+j*N+i)*Nn+n] *= b0;
				if (Indbf != 0)
					ChiBez_rst[(Indbf*N*N+j*N+i)*Nn+n] += b1*ChiBez_rst[((Indbf-1)*N*N+j*N+i)*Nn+n];
			}
		}}}}
	}

	// Transpose ChiBez_rst
	mkl_dimatcopy('R','T',Nbf,Nn,1.,ChiBez_rst,Nn,Nbf);

	*NbfOut = Nbf;
	return ChiBez_rst;
}

static void set_Exp_entry(unsigned int *Nperms, unsigned int *Exp_ptr, const unsigned int *ExponentFound,
                          const unsigned int d)
{
	for (size_t i = 0; i < d+1; i++)
		Exp_ptr[i] = ExponentFound[i];

	// Set Nperms based on number of unique entries
	unsigned int Nunique = 1, Indunique;

	for (size_t i = 1; i < d+1; i++) {
		if (Exp_ptr[i] != Exp_ptr[i-1]) {
			Indunique = i;
			Nunique++;
		}
	}

	if (d == 2) {
		if      (Nunique == 1) *Nperms = 1;
		else if (Nunique == 2) *Nperms = 3;
		else if (Nunique == 3) *Nperms = 6;
		else
			printf("Error: Unsupported.\n"), EXIT_MSG;
	} else if (d == 3) {
		if        (Nunique == 1) { *Nperms = 1;
		} else if (Nunique == 2) {
			// Note the two configurations in this case
			if      (Indunique == 1 || Indunique == 3) *Nperms = 4;
			else if (Indunique == 2)                   *Nperms = 6;
			else
				printf("Error: Unsupported (%d).\n",Indunique), EXIT_MSG;
		} else if (Nunique == 3) { *Nperms = 12;
		} else if (Nunique == 4) { *Nperms = 24;
		} else { printf("Error: Unsupported.\n"), EXIT_MSG; }
	} else {
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}
}

void get_BCoord_Exponents(const unsigned int P, const unsigned int d, unsigned int *NExp, unsigned int **NpermsOut,
                          unsigned int **ExponentsOut)
{
	/*
	 *	Comments:
	 *		The exponents are ordered such that Exp[0] >= Exp[1] >= ... >= Exp[d+1].
	 *		If it is found to be slow, a recursive implementation of this function may be possible.
	 */

	unsigned int Nbf, *Nperms, *Exponents, *Exp_current;

	Nbf = (unsigned int) (factorial_ull(d+P)/(factorial_ull(d)*factorial_ull(P)));

	Nperms      = malloc(Nbf       * sizeof *Nperms);      // keep
	Exponents   = malloc(Nbf*(d+1) * sizeof *Exponents);   // keep
	Exp_current = malloc((d+1)     * sizeof *Exp_current); // free


	unsigned int ExponentFound[d+1], Ind = 0;
	if (d == 2) {
		for (size_t i = 0, iMax = P; i <= iMax; i++) {
			unsigned int sum_i = i;
			for (size_t j = 0, jMax = (d > 0)*i; j <= jMax; j++) {
				unsigned int sum_j = sum_i + j;
				if (sum_j > P)
					break;

				size_t k = P-sum_j;
				if (k > j)
					continue;
//				printf("%zu %zu %zu\n",i,j,k);

				// Found a valid combination
				ExponentFound[0] = i;
				ExponentFound[1] = j;
				ExponentFound[2] = k;
				set_Exp_entry(&Nperms[Ind],&Exponents[Ind*(d+1)],ExponentFound,d);

				Ind++;
			}
		}
	} else if (d == 3) {
		for (size_t i = 0, iMax = P; i <= iMax; i++) {
			unsigned int sum_i = i;
			for (size_t j = 0, jMax = (d > 0)*i; j <= jMax; j++) {
				unsigned int sum_j = sum_i + j;
				if (sum_j > P)
					break;
				for (size_t k = 0, kMax = (d > 1)*j; k <= kMax; k++) {
					unsigned int sum_k = sum_j + k;
					if (sum_k > P)
						break;

					size_t l = P-sum_k;
					if (l > k)
						continue;
//					printf("%zu %zu %zu %zu\n",i,j,k,l);

					// Found a valid combination
					ExponentFound[0] = i;
					ExponentFound[1] = j;
					ExponentFound[2] = k;
					ExponentFound[3] = l;
					set_Exp_entry(&Nperms[Ind],&Exponents[Ind*(d+1)],ExponentFound,d);

					Ind++;
				}
			}
		}
	} else {
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}

	// Check that sum of Nperms == Nbf (Simplex)
	unsigned int sum = 0;
	for (size_t n = 0; n < Ind; n++)
		sum += Nperms[n];

	if (sum != Nbf)
		printf("Error: Did not find the correct number of permutations (%d %d).\n",Nbf,sum), EXIT_MSG;

	*NExp         = Ind;
	*NpermsOut    = Nperms;
	*ExponentsOut = Exponents;
	free(Exp_current);
}

static unsigned int get_Index_Exp_pm1(const unsigned int *Exp_Indices, const unsigned int *basis_exp,
                                      const unsigned int reduced_ind, const unsigned int d, const unsigned int Nbf_p)
{
	/*
	 *	Comments:
	 *		It is assumed that Exp_Indices is sorted (using array_sort) to find the column index efficiently.
	 */

	unsigned int reduced_exp[d+1], IndCol;

	for (size_t i = 0; i < d+1; i++) {
		if (i != reduced_ind)
			reduced_exp[i] = basis_exp[i];
		else
			reduced_exp[i] = basis_exp[i]-1;
	}

	IndCol = 0;
	for (size_t n = 0; n < Nbf_p; n++) {
		for ( ; IndCol <= d && Exp_Indices[n*(d+2)+IndCol] == reduced_exp[IndCol]; IndCol++)
			;
		if (IndCol > d)
			return Exp_Indices[n*(d+2)+IndCol];
	}
	printf("Error: Did not find index.\n"), EXIT_MSG;
	return UINT_MAX;
}

double *basis_SI_Bezier(const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut,
                        const unsigned int d)
{
	/*
	 *	Comments:
	 *		The recursive definition of the basis functions is used in the computation below.
	 *		The ordering of the basis is such that the symmetries of the simplex are maintained and no extra memory is
	 *		needed for the recursive definition of the basis.
	 *
	 *		See Chan(2016) for a generalized definition of Bernstein-Bezier basis functions for the simplex (and
	 *		pyramid) which may also exploit sum factorization ideas.
	 *		Likely change the 2D basis definition here to the tensor-product form (section 2, Chan(2016)). ToBeDeleted
	 *
	 *	References:
	 *		Prautzsch(2002)-Bezier_and_B-Spline_Techniques (Ch. 10.1)
	 *		Chan(2016)-A_Short_Note_on_a_Bernstein-Bezier_Basis_for_the_Pyramid
	 */

	if (d < 2 || d > 3)
		printf("Error: basis_SI_Bezier only supports d = [2,3].\n"), EXIT_MSG;

	unsigned int Nbf, *Exp_Indices[2];
	double       *ChiBez_rst, *BCoords;

	Nbf = (unsigned int) (factorial_ull(d+P)/(factorial_ull(d)*factorial_ull(P)));

	// Convert from rst to barycentric coordinates
	BCoords = malloc(Nn*(d+1) * sizeof *BCoords); // free
	rst_to_barycentric_SI(Nn,d,rst,BCoords);

	Exp_Indices[0] = calloc(Nbf*(d+2) , sizeof *Exp_Indices[0]); // free
	Exp_Indices[1] = calloc(Nbf*(d+2) , sizeof *Exp_Indices[1]); // free

	ChiBez_rst = calloc(Nn*Nbf , sizeof *ChiBez_rst); // keep

	// p = 0 basis
	for (size_t n = 0; n < Nn; n++)
		ChiBez_rst[n] = 1.0;

	// Exp_Indices == zeros(Nbf,d+1) for p = 0.

	if (d == 2) {
		unsigned int permutationsTRI[18]  = { 0, 1, 2, 2, 0, 1, 1, 2, 0, 0, 2, 1, 1, 0, 2, 2, 1, 0};

		for (size_t p = 1; p <= P; p++) {
			unsigned int Nbf_p, Indbf;

			Nbf_p = (unsigned int) (factorial_ull(d+p)/(factorial_ull(d)*factorial_ull(p)));
			Indbf = Nbf_p - 1;

			unsigned int NExp, *Nperms, *Exponents;
			get_BCoord_Exponents(p,d,&NExp,&Nperms,&Exponents); // free

			for (size_t e = 0; e < NExp; e++) {
				unsigned int *Exp_ptr = &Exponents[e*(d+1)];

				// Find number of permutations available for current Exp exponents (note descending values)
				for (size_t perm = Nperms[e]; perm ; perm--) {
					unsigned int basis_exp[d+1], *perm_ptr;

					perm_ptr = &permutationsTRI[(perm-1)*(d+1)];
					for (size_t i = 0; i < d+1; i++)
						basis_exp[i] = Exp_ptr[perm_ptr[i]];

					// Find the indices of the lower order bases to be used and sort them such that the highest index is
					// used first. This is required as memory is overwritten in the recursive basis computation.
					unsigned int IndBCoord[d+1], Indpm1[d+1];
					for (size_t i = 0; i <= d; i++) {
						IndBCoord[i] = i;
						if (basis_exp[i] == 0) {
							Indpm1[i] = 0;
							continue;
						}
						Indpm1[i] = get_Index_Exp_pm1(Exp_Indices[1],basis_exp,i,d,Nbf_p);
					}
					array_sort_ui(1,d+1,Indpm1,IndBCoord,'R','N');

					unsigned int Entered = 0;
					for (size_t i = d; i+1; i--) { // Note: Use the highest Indpm1 first
						if (basis_exp[IndBCoord[i]] == 0)
							continue;

						// Find indices of lower order basis functions used to construct the current basis function.

						if (Entered && Indpm1[i] >= Indbf)
							printf("Error: Using modified values.\n"), EXIT_MSG;

						if (!Entered) {
							for (size_t n = 0; n < Nn; n++)
								ChiBez_rst[Indbf*Nn+n] =  BCoords[IndBCoord[i]*Nn+n]*ChiBez_rst[Indpm1[i]*Nn+n];
						} else {
							for (size_t n = 0; n < Nn; n++)
								ChiBez_rst[Indbf*Nn+n] += BCoords[IndBCoord[i]*Nn+n]*ChiBez_rst[Indpm1[i]*Nn+n];
						}
						Entered = 1;

						for (size_t j = 0; j <= d; j++)
							Exp_Indices[0][Indbf*(d+2)+j] = basis_exp[j];
						Exp_Indices[0][Indbf*(d+2)+(d+1)] = Indbf;

					}
					Indbf--;
				}
			}
			free(Nperms);
			free(Exponents);

			// Sort Exp_Indices
			unsigned int Ind_ui[Nbf_p];
			for (size_t i = 0; i < Nbf_p; i++)
				Ind_ui[Nbf_p] = i;
			for (size_t i = 0, iMax = (d+2)*Nbf_p; i < iMax; i++)
				Exp_Indices[1][i] = Exp_Indices[0][i];

			array_sort_ui(Nbf_p,d+2,Exp_Indices[1],Ind_ui,'R','T');
		}
	} else if (d == 3) {
		printf("Add support.\n"), EXIT_MSG;
	}

	free(BCoords);
	free(Exp_Indices[0]);
	free(Exp_Indices[1]);

	// Transpose ChiBez_rst
	mkl_dimatcopy('R','T',Nbf,Nn,1.,ChiBez_rst,Nn,Nbf);

	*NbfOut = Nbf;
	return ChiBez_rst;
}

/*
 *	Purpose:
 *		Convert from rst to barycentric coordinates.
 *
 *	Comments:
 *		See comments in 'cubature_PYR' for the algorithm to determine BCoords.
 *
 *	Notation:
 *
 *	References:
 */

void rst_to_barycentric_SI(const unsigned int Nn, const unsigned int d, const double *rst, double *BCoords)
{
	if (d < 2 || d > 3)
		printf("Error: Conversion from rst to barycentric SI coordinates only supported for d = [2,3].\n"), EXIT_MSG;

	unsigned int Nve;
	double       *rst_vV, *A_Nve, *A_Nn, *A_Nve_Inv, *I_Nve;

	struct S_ELEMENT *ELEMENT = NULL;

	if (d == 2)
		ELEMENT = get_ELEMENT_type(TRI);
	else if (d == 3)
		ELEMENT = get_ELEMENT_type(TET);

	Nve = ELEMENT->Nve;

	rst_vV = get_rst_vV(ELEMENT); // free

	A_Nve = malloc(Nve*(d+1) * sizeof *A_Nve); // free
	A_Nn  = malloc(Nn*(d+1)  * sizeof *A_Nve); // free

	for (size_t n = 0; n < Nve; n++)
		A_Nve[n] = 1.0;
	for (size_t n = 0, nMax = Nve*d; n < nMax; n++)
		A_Nve[Nve+n] = rst_vV[n];
	free(rst_vV);

	I_Nve     = identity_d(Nve);                // free
	A_Nve_Inv = inverse_d(Nve,Nve,A_Nve,I_Nve); // free
	free(I_Nve);
	free(A_Nve);

	for (size_t n = 0; n < Nn; n++)
		A_Nn[n] = 1.0;
	for (size_t n = 0, nMax = Nn*d; n < nMax; n++)
		A_Nn[Nn+n] = rst[n];

	mm_d(CBCM,CBNT,CBNT,Nn,Nve,Nve,1.0,0.0,A_Nn,A_Nve_Inv,BCoords);

	free(A_Nve_Inv);
	free(A_Nn);
}

struct S_MATRIX *basis_mat (unsigned int const P, struct S_CUBATURE const *const CUBDATA, basis_tdef basis)
{
	struct S_MATRIX *ChiRef_rst = malloc(sizeof *ChiRef_rst); // keep

	unsigned int Nbf = 0;
	ChiRef_rst->values = basis(P,CUBDATA->rst,CUBDATA->Nn,&Nbf,CUBDATA->d); // keep

	ChiRef_rst->format = 'D';
	ChiRef_rst->layout = 'R';
	ChiRef_rst->NRows  = CUBDATA->Nn;
	ChiRef_rst->NCols  = Nbf;

	return ChiRef_rst;
}
