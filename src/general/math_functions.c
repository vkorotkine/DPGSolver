// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/** \file
 */

#include "math_functions.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "macros.h"
#include "definitions_tol.h"

// Interface functions ********************************************************************************************** //

long long unsigned int factorial_ull (const unsigned int n)
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

double factorial_d (const unsigned int n)
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

double gamma_d (const double x)
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

bool equal_d (const double x0, const double x1, const double tol)
{
	if ((fabs(x0) < tol && fabs(x0-x1) < tol) ||
	    (fabs((x0-x1)/x0) < tol))
		return true;
	return false;
}

double norm_d (const ptrdiff_t n_entries, const double*const data, const char*const norm_type)
{
	double norm = 0.0;
	if (strstr(norm_type,"L2")) {
		for (ptrdiff_t i = 0; i < n_entries; ++i)
			norm += data[i]*data[i];
		return sqrt(norm);
	}
	EXIT_UNSUPPORTED;
}

double norm_diff_d
	(const ptrdiff_t n_entries, const double*const data_0, const double*const data_1, const char*const norm_type)
{
	double norm_num = 0.0,
	       norm_den = 0.0;

	if (strstr(norm_type,"Inf")) {
		for (ptrdiff_t i = 0; i < n_entries; ++i) {
			const double diff = fabs(data_0[i]-data_1[i]);
			if (diff > norm_num)
				norm_num = diff;

			const double max = max_abs_d(data_0[i],data_1[i]);
			if (max > norm_den)
				norm_den = max;
		}
	} else {
		EXIT_UNSUPPORTED;
	}

	return ( fabs(norm_den) > EPS ? norm_num/norm_den : norm_num );
}

double max_abs_d (const double a, const double b)
{
	const double a_abs = fabs(a),
	             b_abs = fabs(b);
	return ( a_abs > b_abs ? a_abs : b_abs );
}

int min_i (const int a, const int b)
{
	return ( a < b ? a : b );
}

int max_i (const int a, const int b)
{
	return ( a > b ? a : b );
}
