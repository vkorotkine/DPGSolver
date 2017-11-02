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
 *  \brief Provides Multiarray_\* constructors and destructors for complex variable multiarrays.
 *
 *  See \ref multiarray_constructors.h for potentially relevant comments.
 */

#include <stddef.h>
#include <stdbool.h>

// Default constructors ********************************************************************************************* //

// Empty constructors *********************************************************************************************** //

/** \brief Constructor for an empty \ref Multiarray_d\*.
 *  \return Standard. */
struct Multiarray_c* constructor_empty_Multiarray_c
	(const char layout,              ///< Defined in \ref Multiarray_d.
	 const int order,                ///< Defined in \ref Multiarray_d.
	 const ptrdiff_t*const extents_i ///< The input extents.
	);

// Zero constructors ************************************************************************************************ //

// Copy constructors ************************************************************************************************ //

// Move constructors ************************************************************************************************ //

// Special constructors ********************************************************************************************* //

// Destructors ****************************************************************************************************** //

/// \brief Destructs a \ref Multiarray_c\*.
void destructor_Multiarray_c
	(struct Multiarray_c* a ///< Standard.
	);

#endif // DPG__complex_multiarray_constructors_h__INCLUDED
