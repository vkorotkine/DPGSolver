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

#ifndef DPG__complex_multiarray_math_h__INCLUDED
#define DPG__complex_multiarray_math_h__INCLUDED
/** \file
 *  \brief Provides Multiarray_c math functions.
 */

#include <stdbool.h>

struct Multiarray_c;
struct const_Matrix_d;
struct const_Multiarray_c;

/// \brief `complex` version of \ref transpose_Multiarray_d.
void transpose_Multiarray_c
	(struct Multiarray_c* a, ///< Defined for \ref transpose_Multiarray_d.
	 const bool mem_only     ///< Defined for \ref transpose_Multiarray_d.
	);

/// \brief `complex` version of \ref mm_NNC_Multiarray_d.
void mm_NNC_Multiarray_c
	(const double alpha,                      ///< Defined for \ref mm_NNC_Multiarray_d.
	 const double beta,                       ///< Defined for \ref mm_NNC_Multiarray_d.
	 const struct const_Matrix_d*const a,     ///< Defined for \ref mm_NNC_Multiarray_d.
	 const struct const_Multiarray_c*const b, ///< Defined for \ref mm_NNC_Multiarray_d.
	 struct Multiarray_c*const c              ///< Defined for \ref mm_NNC_Multiarray_d.
	);

#endif // DPG__complex_multiarray_math_h__INCLUDED
