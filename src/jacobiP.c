#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "functions.h"

/*
 *	Purpose:
 *		Evaluate Jacobi Polynomial of type (alpha,beta) > -1, (alpha+beta <> -1) at the point x for order N.
 *
 *	Comments:
 *		Polynomials are normalized to be orthonormal.
 *		If, whille profiling, it is determined that a lot of time is spent in this routine, change it to allow a vector
 *		input for x.
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
