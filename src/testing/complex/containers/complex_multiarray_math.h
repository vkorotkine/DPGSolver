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
#include <complex.h>

struct Multiarray_c;
struct const_Vector_i;
struct const_Vector_d;
struct const_Matrix_d;
struct const_Multiarray_c;

/// \brief `complex` version of \ref transpose_Multiarray_d.
void transpose_Multiarray_c
	(struct Multiarray_c* a, ///< See brief.
	 const bool mem_only     ///< See brief.
	);

/// \brief `complex` version of \ref permute_Multiarray_d_V.
void permute_Multiarray_c_V
	(struct Multiarray_c* a,           ///< See brief.
	 const struct const_Vector_i* p_V, ///< See brief.
	 const char perm_layout            ///< See brief.
	);

/// \brief `complex` version of \ref scale_Multiarray_d.
void scale_Multiarray_c
	(struct Multiarray_c* a,  ///< See brief.
	 const double complex val ///< See brief.
	);

/// \brief `complex` version of \ref scale_Multiarray_by_Vector_d.
void scale_Multiarray_c_by_Vector_d
	(const char side,                     ///< See brief.
	 const double alpha,                  ///< See brief.
	 struct Multiarray_c*const a,         ///< See brief.
	 const struct const_Vector_d*const b, ///< See brief.
	 const bool invert_diag               ///< See brief.
	);

/// \brief `complex` version of \ref mm_NNC_Multiarray_d.
void mm_NNC_Multiarray_c
	(const double alpha,                      ///< See brief.
	 const double beta,                       ///< See brief.
	 const struct const_Matrix_d*const a,     ///< See brief.
	 const struct const_Multiarray_c*const b, ///< See brief.
	 struct Multiarray_c*const c              ///< See brief.
	);

#endif // DPG__complex_multiarray_math_h__INCLUDED
