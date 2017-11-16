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
#include <stddef.h>

struct Matrix_c;
struct const_Vector_d;
struct const_Matrix_d;
struct const_Matrix_c;
struct const_Vector_c;

/// \brief `complex` version of \ref transpose_Matrix_d.
void transpose_Matrix_c
	(struct Matrix_c* a, ///< See brief.
	 const bool mem_only ///< See brief.
	);

/// \brief `complex` version of \ref permute_Matrix_d.
void permute_Matrix_c
	(struct Matrix_c* a, ///< See brief.
	 const ptrdiff_t* p  ///< See brief.
	);

/// \brief `complex` version of \ref scale_Matrix_by_Vector_d.
void scale_Matrix_c_by_Vector_d
	(const char side,                     ///< See brief.
	 const double alpha,                  ///< See brief.
	 struct Matrix_c*const a,             ///< See brief.
	 const struct const_Vector_d*const b, ///< See brief.
	 const bool invert_diag               ///< See brief.
	);

/// \brief `complex` version of \ref mm_d (`double`, `double complex`, `double complex`).
void mm_dcc
	(const char trans_a_i,                ///< See brief.
	 const char trans_b_i,                ///< See brief.
	 const double alpha,                  ///< See brief.
	 const double beta,                   ///< See brief.
	 const struct const_Matrix_d*const a, ///< See brief.
	 const struct const_Matrix_c*const b, ///< See brief.
	 struct Matrix_c*const c              ///< See brief.
	);

/// \brief `complex` version of \ref mm_d (`double complex`, `double`, `double complex`).
void mm_cdc
	(const char trans_a_i,                ///< See brief.
	 const char trans_b_i,                ///< See brief.
	 const double alpha,                  ///< See brief.
	 const double beta,                   ///< See brief.
	 const struct const_Matrix_c*const a, ///< See brief.
	 const struct const_Matrix_d*const b, ///< See brief.
	 struct Matrix_c*const c              ///< See brief.
	);

/// \brief `complex` version of \ref mm_d (`double complex`, `double complex`, `double complex`).
void mm_ccc
	(const char trans_a_i,                ///< See brief.
	 const char trans_b_i,                ///< See brief.
	 const double alpha,                  ///< See brief.
	 const double beta,                   ///< See brief.
	 const struct const_Matrix_c*const a, ///< See brief.
	 const struct const_Matrix_c*const b, ///< See brief.
	 struct Matrix_c*const c              ///< See brief.
	);

/// \brief `complex` version of \ref mm_diag_d (`double`, `double complex`, `double complex`).
void mm_diag_c
	(const char side,                     ///< See brief.
	 const double alpha,                  ///< See brief.
	 const double beta,                   ///< See brief.
	 const struct const_Matrix_d*const a, ///< See brief.
	 const struct const_Vector_c*const b, ///< See brief.
	 struct Matrix_c* c,                  ///< See brief.
	 const bool invert_diag               ///< See brief.
	);

#endif // DPG__complex_matrix_math_h__INCLUDED
