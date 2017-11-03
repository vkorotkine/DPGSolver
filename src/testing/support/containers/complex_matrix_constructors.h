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
 *
 *  Matrices are 2D Multiarrays.
 */

#include <stddef.h>
#include <stdbool.h>
#include <complex.h>

struct const_Matrix_d;

// Default constructors ********************************************************************************************* //

// Empty constructors *********************************************************************************************** //

// Copy constructors ************************************************************************************************ //

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
