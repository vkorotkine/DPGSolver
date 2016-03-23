#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//#include "parameters.h"
#include "functions.h"

//#include "petscsys.h"
//#include "mkl.h"

/*
 *	Purpose:
 *		Return the "matrix" Chi_v representing the orthonormal basis functions evaluated at the provided quadrature
 *		nodes for polynomial order P.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 *
 */

double *basis_TP(const int P, const double *xir, const int Nvn, const int d)
{
	int    i, j, k, Indbf, 
	       N, Nbf;
	double *Chi_v;

    N = P+1;
	Nbf = pow(N,d);

	Chi_v = malloc(Nvn*Nbf * sizeof *Chi_v); // keep

// ToBeDeleted:
for (i = 0; i < Nvn*Nbf; i++) Chi_v[i] = i;

	Indbf = 0;
	if (d == 1) {
		for (i = 0; i < N; i++) {

			Indbf++;
		}
	} else if (d == 2) {
		for (j = 0; j < N; j++) {
		for (i = 0; i < N; i++) {

			Indbf++;
		}}
	} else if (d == 3) {
		for (k = 0; k < N; k++) {
		for (j = 0; j < N; j++) {
		for (i = 0; i < N; i++) {

			Indbf++;
		}}}
	}

double jtest;
jtest = jacobiP(xir[0],0.,0.,N);


	return Chi_v;
}
