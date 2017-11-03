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

/** \brief Compares input values for approximate equality using the relative infinity norm.
 *  \return `true` if the difference is less than the tolerance. */
bool equal_d
	(const double x0, ///< Input 0.
	 const double x1, ///< Input 1.
	 const double tol ///< The tolerance.
	);

/** \brief Computes the norm of the input `double*` data with the specified norm type.
 *  \return See brief. */
double norm_d
	(const ptrdiff_t n_entries, ///< The number of entries.
	 const double*const data,   ///< The data.
	 const char*const norm_type ///< The norm type. Options: "L2", "Inf".
	);

/** \brief Computes the relative norm of the difference between the input `double*` data with the specified norm type.
 *  \return See brief. */
double norm_diff_d
	(const ptrdiff_t n_entries, ///< The number of entries.
	 const double*const data_0, ///< The data for input 0.
	 const double*const data_1, ///< The data for input 1.
	 const char*const norm_type ///< The norm type. Options: "Inf".
	);

/** \brief Compute the maximum absolute value of the two inputs and return it.
 *  \return See brief. */
double max_abs_d
	(const double a, ///< Input 0.
	 const double b  ///< Input 1.
	);

/** \brief Compute the binomial coefficient.
 *  \return See brief. */
double binomial_coef
	(const int num, ///< The numerator.
	 const int den  ///< The denominator.
	);

/// \brief Variation on the axpy (BLAS 1) function: z = y*x + z.
void z_yxpz
	(const int n,     ///< The number of entries.
	 const double* x, ///< Input x.
	 const double* y, ///< Input y.
	 double* z        ///< Location to store the sum.
	);

#endif // DPG__math_functions_h__INCLUDED
