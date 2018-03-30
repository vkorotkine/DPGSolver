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

#ifndef DPG__test_support_vector_h__INCLUDED
#define DPG__test_support_vector_h__INCLUDED
/** \file
 *  \brief Provides support functions for testing relating to the containers defined in \ref vector.h.
 */

#include <stdio.h>
#include <stdbool.h>

struct Vector_i;
struct Vector_d;
struct const_Vector_i;
struct const_Vector_d;

// Difference functions ********************************************************************************************* //

/** \brief Check the relative difference between entries in the input \ref Vector_T\*s up to the input tolerance.
 *  \return The `true` if inputs differ; `false` otherwise. */
bool diff_Vector_d
	(const struct Vector_d*const a, ///< Input 0.
	 const struct Vector_d*const b, ///< Input 1.
	 const double tol               ///< The tolerance.
	);

/** \brief `const` version of \ref diff_Vector_d.
 *  \return See brief. */
bool diff_const_Vector_d
	(const struct const_Vector_d*const a, ///< Defined for \ref diff_Vector_d.
	 const struct const_Vector_d*const b, ///< Defined for \ref diff_Vector_d.
	 const double tol                     ///< Defined for \ref diff_Vector_d.
	);

/** \brief Check the difference between entries in the input \ref Vector_T\*s.
 *  \return The `true` if inputs differ; `false` otherwise. */
bool diff_Vector_i
	(const struct Vector_i*const a, ///< Input 0.
	 const struct Vector_i*const b  ///< Input 1.
	);

/** \brief `const` version of \ref diff_Vector_i.
 *  \return See brief. */
bool diff_const_Vector_i
	(const struct const_Vector_i*const a, ///< See brief.
	 const struct const_Vector_i*const b  ///< See brief.
	);

// Printing functions *********************************************************************************************** //

/// \brief Print the relative difference of the input \ref Vector_T\*s, outputting 0 if less than the tolerance.
void print_diff_Vector_d
	(const struct Vector_d*const a, ///< Input 0.
	 const struct Vector_d*const b, ///< Input 1.
	 const double tol               ///< The tolerance.
	);

/// \brief `const` version of \ref print_diff_Vector_d.
void print_diff_const_Vector_d
	(const struct const_Vector_d*const a, ///< Defined for \ref print_diff_Vector_d.
	 const struct const_Vector_d*const b, ///< Defined for \ref print_diff_Vector_d.
	 const double tol                     ///< Defined for \ref print_diff_Vector_d.
	);

/// \brief Print the difference of the input \ref Vector_T\*s.
void print_diff_Vector_i
	(const struct Vector_i*const a, ///< Input 0.
	 const struct Vector_i*const b  ///< Input 1.
	);

/// \brief `const` version of \ref print_diff_Vector_i.
void print_diff_const_Vector_i
	(const struct const_Vector_i*const a, ///< See brief.
	 const struct const_Vector_i*const b  ///< See brief.
	);

#endif // DPG__test_support_vector_h__INCLUDED
