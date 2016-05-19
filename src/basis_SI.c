// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Return the "matrix" ChiRef_rst representing the orthonormal basis functions evaluated at the provided quadrature
 *		nodes for polynomial order P.
 *
 *	Comments:
 *		The "matrix" is returned as a 1D array.
 *
 *	Notation:
 *
 *	References:
 *		Hesthaven(2008)-Nodal_Discontinuous_Galerkin_Methods (Appendix A.1)
 */

double *basis_SI(const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut,
                 const unsigned int d)
{
	unsigned int i, j, k, iMax, jMax, kMax, Indbf, Indn, Nbf;
	double       *ChiRef_rst, *a, *b, *c;

	if (d < 2 || d > 3)
		printf("Error: basis_SI only supports d = [2,3].\n"), exit(1);

	// Convert from rst to abc coordinates
	a = malloc(Nn * sizeof *a); // free
	b = malloc(Nn * sizeof *b); // free
	c = malloc(Nn * sizeof *c); // free

	rst_to_abc_SI(Nn,d,rst,a,b,c);

	Nbf = factorial_ull(d+P)/(factorial_ull(d)*factorial_ull(P));

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
