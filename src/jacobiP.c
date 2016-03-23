#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//#include "parameters.h"
//#include "functions.h"

//#include "petscsys.h"
//#include "mkl.h"

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
 *
 */

double jacobiP(const double x, const double alpha, const double beta, const int N)
{
	double gamma0,
	       *PL;

	PL = malloc((N+1) * sizeof *PL); // tbd

	// Initial values P_0(x) and P_1(x)
	gamma0 = pow(2,alpha+beta+1)/(alpha+beta+1); // * gamma function...

}
