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
struct const_Matrix_c;

// Default constructors ********************************************************************************************* //

// Empty constructors *********************************************************************************************** //

/** \brief `complex` version of \ref constructor_empty_Vector_d.
 *  \return Standard. */
struct Vector_c* constructor_empty_Vector_c
	(const ptrdiff_t ext_0 ///< Defined in \ref Vector_d.
	);

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

/** \brief `complex` version of \ref constructor_mv_Vector_d.
 *  \return Standard. */
struct Vector_c* constructor_mv_Vector_c
	(const char trans_a_i,                ///< See brief.
	 const double alpha,                  ///< See brief.
	 const struct const_Matrix_c*const a, ///< See brief.
	 const struct const_Vector_c*const b  ///< See brief.
	);

/** \brief `const` version of \ref constructor_mv_Vector_c.
 *  \return Standard. */
const struct const_Vector_c* constructor_mv_const_Vector_c
	(const char trans_a_i,                ///< See brief.
	 const double alpha,                  ///< See brief.
	 const struct const_Matrix_c*const a, ///< See brief.
	 const struct const_Vector_c*const b  ///< See brief.
	);

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
