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

#ifndef DPG__complex_multiarray_h__INCLUDED
#define DPG__complex_multiarray_h__INCLUDED
/** \file
 *  \brief Provides \ref Multiarray_c related functions for complex variables.
 *
 *  See \ref multiarray.h and \ref complex_multiarray.h for potentially relevant comments.
 */

#include "complex_multiarray_minimal.h"

#include "complex_multiarray_constructors.h"
#include "complex_multiarray_math.h"
#include "complex_multiarray_print.h"

struct const_Multiarray_d;

// Interface functions ********************************************************************************************** //

/** \brief `complex` version of \ref get_col_Multiarray_d.
 *  \return See brief. */
double complex* get_col_Multiarray_c
	(const ptrdiff_t col,   ///< See brief.
	 struct Multiarray_c* a ///< See brief.
	);

/** \brief `const` version of \ref get_col_Multiarray_c.
 *  \return See brief. */
const double complex* get_col_const_Multiarray_c
	(const ptrdiff_t col,               ///< See brief.
	 const struct const_Multiarray_c* a ///< See brief.
	);

/// \brief `complex` version of \ref set_to_value_Multiarray_d.
void set_to_value_Multiarray_c
	(struct Multiarray_c*const a, ///< See brief.
	 const double complex val     ///< See brief.
	);

/// \brief Set the data of the \ref Multiarray_c\* container from that of the \ref Multiarray_d\* container.
void set_Multiarray_c_Multiarray_d
	(struct Multiarray_c* a,            ///< Multiarray with data to be set.
	 const struct const_Multiarray_d* b ///< Multiarray from which to take data.
	);

#endif // DPG__complex_multiarray_h__INCLUDED
