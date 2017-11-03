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

#ifndef DPG__complex_matrix_math_h__INCLUDED
#define DPG__complex_matrix_math_h__INCLUDED
/** \file
 *  \brief Provides Matrix_\* math functions.
 */

#include <stdbool.h>

struct Matrix_c;
struct const_Matrix_d;
struct const_Matrix_c;

/// \brief `complex` version of \ref transpose_Matrix_d.
void transpose_Matrix_c
	(struct Matrix_c* a, ///< Defined for \ref transpose_Matrix_d.
	 const bool mem_only ///< Defined for \ref transpose_Matrix_d.
	);

/// \brief `complex` version of \ref mm_d (`double`,`double complex`,`double complex`).
void mm_c
	(const char trans_a_i,                ///< Defined for \ref mm_d.
	 const char trans_b_i,                ///< Defined for \ref mm_d.
	 const double alpha,                  ///< Defined for \ref mm_d.
	 const double beta,                   ///< Defined for \ref mm_d.
	 const struct const_Matrix_d*const a, ///< Defined for \ref mm_d.
	 const struct const_Matrix_c*const b, ///< Defined for \ref mm_d.
	 struct Matrix_c*const c              ///< Defined for \ref mm_d.
	);

#endif // DPG__complex_matrix_math_h__INCLUDED
