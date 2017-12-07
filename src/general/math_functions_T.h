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
 *  \brief Provides several standard math functions (templated).
 */

#include <stdbool.h>
#include <stddef.h>
#include <complex.h>

/** \brief Compares input values for approximate equality using the relative infinity norm.
 *  \return `true` if the difference is less than the tolerance. */
bool equal_T
	(const Type x0, ///< Input 0.
	 const Type x1, ///< Input 1.
	 const Real tol ///< The tolerance.
	);

/** \brief Computes the norm of the input `Type*` data with the specified norm type.
 *  \return See brief. */
Type norm_T
	(const ptrdiff_t n_entries, ///< The number of entries.
	 const Type*const data,     ///< The data.
	 const char*const norm_type ///< The norm type. Options: "L2", "Inf".
	);

/** \brief Computes the relative norm of the difference between the input `Type*` data with the specified norm type.
 *  \return See brief. */
Type norm_diff_T
	(const ptrdiff_t n_entries, ///< The number of entries.
	 const Type*const data_0,   ///< The data for input 0.
	 const Type*const data_1,   ///< The data for input 1.
	 const char*const norm_type ///< The norm type. Options: "Inf".
	);

/** \brief Compute the maximum absolute value of the two inputs and return it.
 *  \return See brief. */
Type max_abs_T
	(const Type a, ///< Input 0.
	 const Type b  ///< Input 1.
	);

/// \brief Variation on the axpy (BLAS 1) function: z = y*x + z.
void z_yxpz_T
	(const int n,   ///< The number of entries.
	 const Type* x, ///< Input x.
	 const Type* y, ///< Input y.
	 Type* z        ///< Location to store the sum.
	);

/// \brief Version of \ref z_yxpz_T with input types: `Real`, `Type`, `Type`.
void z_yxpz_RTT
	(const int n,   ///< The number of entries.
	 const Real* x, ///< Input x.
	 const Type* y, ///< Input y.
	 Type* z        ///< Location to store the sum.
	);

