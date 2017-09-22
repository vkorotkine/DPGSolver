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

#ifndef DPG__multiarray_math_h__INCLUDED
#define DPG__multiarray_math_h__INCLUDED
/** \file
 *  \brief Provides Multiarray_\* math functions.
 */

#include <stdbool.h>

struct Multiarray_d;

/** \brief Transpose the \ref Multiarray_d\* optionally leaving the values of the extents unchanged if
 *         `mem_only = true`.
 *
 *  \note The multiarray must have `order = 2`.
 */
void transpose_Multiarray_d
	(struct Multiarray_d* a, ///< Multiarray to be transposed.
	 const bool mem_only     ///< Flag for whether only the memory should be transposed (with extents unchanged).
	);

#endif // DPG__multiarray_math_h__INCLUDED
