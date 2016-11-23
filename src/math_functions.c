// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "math_functions.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Parameters.h"
#include "Macros.h"

/*
 *	Purpose:
 *		Provide several standard math functions.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 *		Press(1992-2nd)_Numerical recipes in C- the art of scientific computing (Ch. 6.1)
 */

long long unsigned int factorial_ull(const unsigned int n)
{
	static unsigned int ntop = 0;
	static double       a[21] = { 1.0 };
	unsigned int i;

	// Note: large values overflow
	if (n > 20)
		printf("Large inputs (n > 20) result in overflow for factorial_ull, n = %d.\n",n), exit(1);

	// As ntop and a are static variables, multiple calls do not result in recomputation.
	while (ntop < n) {
		i = ntop++;
		a[ntop] = a[i]*ntop;
	}
	return (long long unsigned int) a[n];
}

double factorial_d(const unsigned int n)
{
	static unsigned int ntop = 0;
	static double       a[33] = { 1.0 };
	unsigned int i;

	// Note: large values overflow
	if (n > 32)
		printf("Large inputs (n > 32) result in overflow for factorial_d, n = %d.\n",n), exit(1);
	// Check exact value if required (ToBeDeleted)

	// As ntop and a are static variables, multiple calls do not result in recomputation.
	while (ntop < n) {
		i = ntop++;
		a[ntop] = a[i]*ntop;
	}
	return a[n];
}

double gamma_d(const double x)
{
	if (x <= 0.0)
		printf("Error: Input to gamma_d must be greater than 0.0.\n"), exit(1);

	if (fabs(floor(x)-x) < EPS) {
		// Unsigned integer case
		return factorial_d((unsigned int) (x-1.0));
	} else {
		static double cof[6] = { 76.18009172947146,
		                        -86.50532032941677,
		                         24.01409824083091,
		                        -1.231739572450155,
		                         0.1208650973866179e-2,
		                        -0.5395239384953e-5     };
		unsigned int j;
		double z = x, y = x, tmp = x + 5.5,
		       ser = 1.000000000190015;
		tmp -= (z+0.5)*log(tmp);
		for (j = 0; j <= 5; j++) {
			ser += cof[j]/++y;
		}

		return exp(-tmp+log(2.5066282746310005*ser/z));
	}
}
