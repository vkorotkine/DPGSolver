#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "parameters.h"
#include "functions.h"

#include "mkl.h"

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
 *
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
	s = malloc(0  * sizeof *s); // silence
	t = malloc(0  * sizeof *t); // silence
	if (d > 1) {
		free(s);
		s = malloc(Nn * sizeof *s); // free
	}
	if (d > 2) {
		free(t);
		t = malloc(Nn * sizeof *t); // free
	}

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
