/* {{{
This file is part of DPGSolver.

DPGSolver is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or any later version.

DPGSolver is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along with DPGSolver.  If not, see
<http://www.gnu.org/licenses/>.
}}} */
/** \file
 */

#include "math_functions.h"

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "jacobi.h"
#include "gsl/gsl_sf_gamma.h"

#include "macros.h"
#include "definitions_tol.h"

// Static function declarations ************************************************************************************* //

/** \brief Compute the normalization for the Jacobi polynomials such that they are orthonormal.
 *  \return See brief.
 *
 *  The normalization was taken from the [Mathematica page on Jacobi polynomials][mathematica_jacobi].
 *
 *  <!-- References: -->
 *  [mathematica_jacobi]: http://mathworld.wolfram.com/JacobiPolynomial.html
 */
double compute_jacobi_normalization
	(const int n,    ///< Defined for \ref jac_jacobi_normalized.
	 const double a, ///< Defined for \ref jac_jacobi_normalized.
	 const double b  ///< Defined for \ref jac_jacobi_normalized.
	);

// Interface functions ********************************************************************************************** //

double jac_jacobi_normalized (const double x, const int n, const double a, const double b)
{
	return compute_jacobi_normalization(n,a,b)*jac_jacobi(x,n,a,b);
}

double jac_djacobi_normalized (const double x, const int n, const double a, const double b)
{
	return compute_jacobi_normalization(n,a,b)*jac_djacobi(x,n,a,b);
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
	if (strcmp(norm_type,"L2") == 0) {
		for (ptrdiff_t i = 0; i < n_entries; ++i)
			norm += data[i]*data[i];
		return sqrt(norm);
	} else if (strcmp(norm_type,"Inf") == 0) {
		for (ptrdiff_t i = 0; i < n_entries; ++i) {
			if (fabs(data[i]) > norm)
				norm = fabs(data[i]);
		}
		return norm;
	}
	EXIT_UNSUPPORTED;
}

double complex norm_c (const ptrdiff_t n_entries, const double complex*const data, const char*const norm_type)
{
	double complex norm = 0.0;
	if (strcmp(norm_type,"L2") == 0) {
		for (ptrdiff_t i = 0; i < n_entries; ++i)
			norm += data[i]*data[i];
		return csqrt(norm);
	} else if (strcmp(norm_type,"Inf") == 0) {
		for (ptrdiff_t i = 0; i < n_entries; ++i) {
			if (cabs(data[i]) > cabs(norm))
				norm = cabs(data[i]);
		}
		return norm;
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

	return ( fabs(norm_den) > 1e2*EPS ? norm_num/norm_den : norm_num );
}

double max_abs_d (const double a, const double b)
{
	const double a_abs = fabs(a),
	             b_abs = fabs(b);
	return ( a_abs > b_abs ? a_abs : b_abs );
}

double binomial_coef (const int num, const int den)
{
	assert(num >= den);

	return gsl_sf_fact((unsigned)num)/(gsl_sf_fact((unsigned)(num-den))*gsl_sf_fact((unsigned)den));
}

void z_yxpz (const int n, const double* x, const double* y, double* z)
{
	for (int i = 0; i < n; ++i)
		z[i] += y[i]*x[i];
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

double compute_jacobi_normalization (const int n, const double a, const double b)
{
	const double scale = pow(2.0,a+b+1.0)/(2.0*n+a+b+1.0)
	                    *gsl_sf_gamma(n+a+1.0)*gsl_sf_gamma(n+b+1.0)/(gsl_sf_fact((unsigned)n)*gsl_sf_gamma(n+a+b+1.0));
	return 1.0/sqrt(scale);
}
