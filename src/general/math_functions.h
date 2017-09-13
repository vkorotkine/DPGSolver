// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__math_functions_h__INCLUDED
#define DPG__math_functions_h__INCLUDED
/** \file
 *  \brief Provides several standard math functions.
 *
 *  <!-- References: -->
 *  Press(1992,2nd edition) - Numerical Recipes in C - The Art of Scientific Computing (Ch. 6.1)
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
	 const char*const norm_type ///< The norm type. Options: "L2".
	);

/** \brief Computes the relative norm of the difference between the input `double*` data with the specified norm type.
 *  \return See brief. */
double norm_diff_d
	(const ptrdiff_t n_entries, ///< The number of entries.
	 const double*const data_0, ///< The data for input 0.
	 const double*const data_1, ///< The data for input 1.
	 const char*const norm_type ///< The norm type. Options: "L2".
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

#endif // DPG__math_functions_h__INCLUDED
