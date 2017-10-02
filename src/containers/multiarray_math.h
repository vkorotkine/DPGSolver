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
#include <stddef.h>

struct Multiarray_d;
struct const_Matrix_d;
struct const_Multiarray_d;

/// \brief Transpose the \ref Multiarray_d\*'s memory by interpreting as a matrix with `ext_0 = extents[0]`.
void transpose_Multiarray_d
	(struct Multiarray_d* a, ///< Multiarray to be transposed.
	 const bool mem_only     ///< Flag for whether only the memory should be transposed (with extents unchanged).
	);

/// \brief Reinterpret the \ref const_Multiarray_d\* as a \ref const_Matrix_d\* having the given input extents.
void reinterpret_const_Multiarray_as_Matrix_d
	(const struct const_Multiarray_d* a, ///< The multiarray.
	 const struct const_Matrix_d* a_M,   ///< The matrix.
	 const ptrdiff_t ext_0,              ///< The value of `ext_0` for the matrix.
	 const ptrdiff_t ext_1               ///< The value of `ext_1` for the matrix.
	);

/// \brief Reinterpret the \ref const_Matrix_d\* as a \ref const_Multiarray_d\* having the given input order/extents.
void reinterpret_const_Matrix_as_Multiarray_d
	(const struct const_Matrix_d* a_M,   ///< The matrix.
	 const struct const_Multiarray_d* a, ///< The multiarray.
	 const int order,                    ///< The value of `order` for the multiarray.
	 const ptrdiff_t* extents            ///< The value of `extents` for the multiarray.
	);

#endif // DPG__multiarray_math_h__INCLUDED
