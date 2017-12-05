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

// Templated functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "def_templates_math_functions_d.h"
#include "def_templates_math_d.h"
#include "math_functions_T.c"
#include "undef_templates_type.h"
#include "undef_templates_math_functions.h"
#include "undef_templates_math.h"

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

double binomial_coef (const int num, const int den)
{
	assert(num >= den);

	return gsl_sf_fact((unsigned)num)/(gsl_sf_fact((unsigned)(num-den))*gsl_sf_fact((unsigned)den));
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

double compute_jacobi_normalization (const int n, const double a, const double b)
{
	const double scale = pow(2.0,a+b+1.0)/(2.0*n+a+b+1.0)
	                    *gsl_sf_gamma(n+a+1.0)*gsl_sf_gamma(n+b+1.0)/(gsl_sf_fact((unsigned)n)*gsl_sf_gamma(n+a+b+1.0));
	return 1.0/sqrt(scale);
}
