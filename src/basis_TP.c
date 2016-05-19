// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "parameters.h"
#include "functions.h"

#include "mkl.h"

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
 */

double *basis_TP(const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut,
                 const unsigned int d)
{
	unsigned int i, j, k, iMax, jMax, kMax, Indbf, Indn, u1,
	             N, Nbf;
	int          sd, sN;
	double       *ChiRef_rst, *r, *s, *t;

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
