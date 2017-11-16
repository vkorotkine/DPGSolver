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

#ifndef DPG__complex_vector_constructors_h__INCLUDED
#define DPG__complex_vector_constructors_h__INCLUDED
/** \file
 *  \brief Provides Vector_\* constructors and destructors.
 */

#include <stddef.h>
#include <stdbool.h>
#include <complex.h>

struct const_Vector_c;

// Default constructors ********************************************************************************************* //

// Empty constructors *********************************************************************************************** //

// Zero constructors ************************************************************************************************ //

/** \brief `complex` version of \ref constructor_zero_Vector_d.
 *  \return Standard. */
struct Vector_c* constructor_zero_Vector_c
	(const ptrdiff_t ext_0 ///< See brief.
	);

// Copy constructors ************************************************************************************************ //

// Move constructors ************************************************************************************************ //

/** \brief `complex` version of \ref constructor_move_Vector_d_d.
 *  \return Standard. */
struct Vector_c* constructor_move_Vector_c_c
	(const ptrdiff_t ext_0,    ///< See brief.
	 const bool owns_data,     ///< See brief.
	 double complex*const data ///< See brief.
	);

// Set constructors ************************************************************************************************* //

// Special constructors ********************************************************************************************* //

// Destructors ****************************************************************************************************** //

/// \brief Destructs a \ref Vector_c\*.
void destructor_Vector_c
	(struct Vector_c* a ///< Standard.
	);

/// \brief Destructs a \ref const_Vector_c\*.
void destructor_const_Vector_c
	(const struct const_Vector_c* a ///< Standard.
	);

#endif // DPG__complex_vector_constructors_h__INCLUDED
