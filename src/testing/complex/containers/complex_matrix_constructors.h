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

#ifndef DPG__complex_matrix_constructors_h__INCLUDED
#define DPG__complex_matrix_constructors_h__INCLUDED
/** \file
 *  \brief Provides Matrix_\* constructors and destructors.
 */

#include <stddef.h>
#include <stdbool.h>
#include <complex.h>

struct const_Matrix_d;
struct const_Vector_d;

// Default constructors ********************************************************************************************* //

// Empty constructors *********************************************************************************************** //

/** \brief Constructs an empty \ref Matrix_c\*.
 *  \return Standard. */
struct Matrix_c* constructor_empty_Matrix_c
	(const char layout,     ///< Standard.
	 const ptrdiff_t ext_0, ///< Standard.
	 const ptrdiff_t ext_1  ///< Standard.
	);

// Copy constructors ************************************************************************************************ //

/** \brief Copy constructor for a \ref Matrix_c\* from a \ref Matrix_c\*.
 *  \return Standard. */
struct Matrix_c* constructor_copy_Matrix_c
	(const struct Matrix_c* src ///< The source matrix.
	);

/** \brief `const` version of \ref constructor_copy_Matrix_c_Matrix_d.
 *  \return See brief. */
const struct const_Matrix_c* constructor_copy_const_Matrix_c_Matrix_d
	(const struct const_Matrix_d* src ///< See brief.
	);

// Move constructors ************************************************************************************************ //

/** \brief `complex` version of \ref constructor_move_Matrix_d_d.
 *  \return Standard. */
struct Matrix_c* constructor_move_Matrix_c_c
	(const char layout,        ///< See brief.
	 const ptrdiff_t ext_0,    ///< See brief.
	 const ptrdiff_t ext_1,    ///< See brief.
	 const bool owns_data,     ///< See brief.
	 double complex*const data ///< See brief.
	);

/** \brief `const` version of constructor_move_Matrix_c_c.
 *  \return Standard. */
const struct const_Matrix_c* constructor_move_const_Matrix_c_c
	(const char layout,              ///< See brief.
	 const ptrdiff_t ext_0,          ///< See brief.
	 const ptrdiff_t ext_1,          ///< See brief.
	 const bool owns_data,           ///< See brief.
	 const double complex*const data ///< See brief.
	);

// Special constructors ********************************************************************************************* //

/** \brief `complex` version of \ref constructor_mm_Matrix_d (`double complex`, `double complex`).
 *  \return See brief. */
struct Matrix_c* constructor_mm_Matrix_cc
	(const char trans_a_i,                ///< See brief.
	 const char trans_b_i,                ///< See brief.
	 const double alpha,                  ///< See brief.
	 const struct const_Matrix_c*const a, ///< See brief.
	 const struct const_Matrix_c*const b, ///< See brief.
	 const char layout                    ///< See brief.
	);

/** \brief `const` version of \ref constructor_mm_Matrix_cc.
 *  \return See brief.. */
const struct const_Matrix_c* constructor_mm_const_Matrix_cc
	(const char trans_a_i,                ///< See brief.
	 const char trans_b_i,                ///< See brief.
	 const double alpha,                  ///< See brief.
	 const struct const_Matrix_c*const a, ///< See brief.
	 const struct const_Matrix_c*const b, ///< See brief.
	 const char layout                    ///< See brief.
	);

/** \brief `complex` version of \ref constructor_mm_diag_Matrix_d taking a \ref Vector_d input.
 *  \return See brief. */
struct Matrix_c* constructor_mm_diag_Matrix_c_d
	(const double alpha,                  ///< See brief.
	 const struct const_Matrix_c*const a, ///< See brief.
	 const struct const_Vector_d*const b, ///< See brief.
	 const char side,                     ///< See brief.
	 const bool invert_diag               ///< See brief.
	);

/** \brief `const` version of \ref constructor_mm_diag_Matrix_c_d.
 *  \return See brief. */
const struct const_Matrix_c* constructor_mm_diag_const_Matrix_c_d
	(const double alpha,                  ///< See brief.
	 const struct const_Matrix_c*const a, ///< See brief.
	 const struct const_Vector_d*const b, ///< See brief.
	 const char side,                     ///< See brief.
	 const bool invert_diag               ///< See brief.
	);

// Destructors ****************************************************************************************************** //

/// \brief `complex` version of \ref destructor_Matrix_d.
void destructor_Matrix_c
	(struct Matrix_c* a ///< Defined for \ref destructor_Matrix_d.
	);

/// \brief `const` version of \ref destructor_Matrix_c.
void destructor_const_Matrix_c
	(const struct const_Matrix_c* a ///< Defined for \ref destructor_Matrix_c.
	);

#endif // DPG__complex_matrix_constructors_h__INCLUDED
