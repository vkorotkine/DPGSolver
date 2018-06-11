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

#ifndef DPG__math_functions_h__INCLUDED
#define DPG__math_functions_h__INCLUDED
/** \file
 *  \brief Provides several standard math functions.
 */

#include <stdbool.h>
#include <stddef.h>

#include "def_templates_type_d.h"
#include "math_functions_T.h"
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "math_functions_T.h"
#include "undef_templates_type.h"

/** \brief Evaluates an orthonormalized Jacobi polynomial on the standard unit interval ([-1,1]).
 *  \return See brief. */
double jac_jacobi_normalized
	(const double x, ///< The coordinate at which to evaluate the polynomial.
	 const int n,    ///< The order of the polynomial.
	 const double a, ///< alpha.
	 const double b  ///< beta.
	);

/** \brief Evaluates the derivative of an orthonormalized Jacobi polynomial on the standard unit interval ([-1,1]).
 *  \return See brief. */
double jac_djacobi_normalized
	(const double x, ///< Defined for \ref jac_jacobi_normalized.
	 const int n,    ///< Defined for \ref jac_jacobi_normalized.
	 const double a, ///< Defined for \ref jac_jacobi_normalized.
	 const double b  ///< Defined for \ref jac_jacobi_normalized.
	);

/** \brief Compute the binomial coefficient.
 *  \return See brief. */
double binomial_coef
	(const int num, ///< The numerator.
	 const int den  ///< The denominator.
	);

#endif // DPG__math_functions_h__INCLUDED
