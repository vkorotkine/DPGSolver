// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "functions.h"

/*
 *	Purpose:
 *		Evaluate the derivative of the Jacobi Polynomial of type (alpha,beta) > -1, (alpha+beta <> -1) at the point x
 *		for order N.
 *
 *	Comments:
 *		If, whille profiling, it is determined that a lot of time is spent in this routine, change it to allow a vector
 *		input for x.
 *
 *	Notation:
 *
 *	References:
 *		Hesthaven (Nodal DG Code): https://github.com/tcew/nodal-dg
 */

double grad_jacobiP(const double x, const double alpha, const double beta, const int N)
{
	double dP;

	if (N == 0)
		dP = 0.0;
	else
		dP = sqrt(N*(N+alpha+beta+1))*jacobiP(x,alpha+1.0,beta+1.0,N-1);

	return dP;
}
