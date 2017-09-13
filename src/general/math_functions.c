// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/** \file
 */

#include "math_functions.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "jacobi.h"
#include "gsl/gsl_sf_gamma.h"

#include "macros.h"
#include "definitions_tol.h"

// Interface functions ********************************************************************************************** //

double jac_jacobi_normalized (double x, int n, double a, double b)
{
	const double scale = pow(2.0,a+b+1.0)/(2.0*n+a+b+1.0)
	                    *gsl_sf_gamma(n+a+1.0)*gsl_sf_gamma(n+b+1.0)/(gsl_sf_fact(n)*gsl_sf_gamma(n+a+b+1.0));

	return 1.0/sqrt(scale)*jac_jacobi(x,n,a,b);
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
