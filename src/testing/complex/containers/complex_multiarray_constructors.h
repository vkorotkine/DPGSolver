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

#ifndef DPG__complex_multiarray_constructors_h__INCLUDED
#define DPG__complex_multiarray_constructors_h__INCLUDED
/** \file
 *  \brief Provides **remaining** Multiarray_\* constructors and destructors for complex variable multiarrays.
 *
 *  See \ref multiarray_constructors.h and \ref complex_multiarray_constructors.h for potentially relevant comments.
 */

#include <stddef.h>
#include <stdbool.h>
#include <complex.h>

#include "complex_multiarray_minimal_constructors.h"

struct const_Multiarray_d;
struct const_Multiarray_c;

// Default constructors ********************************************************************************************* //

// Empty constructors *********************************************************************************************** //

/** \brief Constructor for an empty \ref Multiarray_c\* with extents having been previously dynamically allocated.
 *  \return Standard. */
struct Multiarray_c* constructor_empty_Multiarray_c_dyn_extents
	(const char layout,            ///< Defined in \ref Multiarray_c.
	 const int order,              ///< Defined in \ref Multiarray_c.
	 const ptrdiff_t*const extents ///< Defined in \ref Multiarray_c.
	);

// Zero constructors ************************************************************************************************ //

/** \brief Same as \ref constructor_empty_Multiarray_c but with data calloc'ed.
 *  \return Standard. */
struct Multiarray_c* constructor_zero_Multiarray_c
	(const char layout,              ///< Defined in \ref Multiarray_d.
	 const int order,                ///< Defined in \ref Multiarray_d.
	 const ptrdiff_t*const extents_i ///< The input extents.
	);

// Copy constructors ************************************************************************************************ //

/** \brief `complex` version of \ref constructor_copy_const_Multiarray_T.
 *  \return Standard. */
const struct const_Multiarray_c* constructor_copy_const_Multiarray_c
	(const struct const_Multiarray_c* src ///< See brief.
	);

/** \brief `const` version of \ref constructor_copy_Multiarray_c_Multiarray_d.
 *  \return Standard. */
const struct const_Multiarray_c* constructor_copy_const_Multiarray_c_Multiarray_d
	(const struct const_Multiarray_d* src ///< See brief.
	);

// Move constructors ************************************************************************************************ //

/** \brief `complex` version of \ref constructor_move_Multiarray_d_d.
 *  \return See brief. */
struct Multiarray_c* constructor_move_Multiarray_c_c
	(const char layout,               ///< See brief.
	 const int order,                 ///< See brief.
	 const ptrdiff_t*const extents_i, ///< See brief.
	 const bool owns_data,            ///< See brief.
	 double complex*const data        ///< See brief.
	);

/** \brief `const` version of \ref constructor_move_Multiarray_c_c.
 *  \return Standard. */
const struct const_Multiarray_c* constructor_move_const_Multiarray_c_c
	(const char layout,               ///< See brief.
	 const int order,                 ///< See brief.
	 const ptrdiff_t*const extents_i, ///< See brief.
	 const bool owns_data,            ///< See brief.
	 const double complex*const data  ///< See brief.
	);

// Special constructors ********************************************************************************************* //

// Destructors ****************************************************************************************************** //

/// \brief `const` version of \ref destructor_Multiarray_c.
void destructor_const_Multiarray_c
	(const struct const_Multiarray_c* a ///< Standard.
	);

#endif // DPG__complex_multiarray_constructors_h__INCLUDED
