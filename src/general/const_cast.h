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

#ifndef DPG__const_cast_h__INCLUDED
#define DPG__const_cast_h__INCLUDED
/** \file
 *  \brief Provides functions for casting to const lvalues for standard datatypes.
 *
 *  Variables with standard datatypes with higher levels of dereferencing than 0 should be placed in Multiarray
 *  containers. Const cast functions are thus not provided for higher levels of dereferencing.
 *
 *  The function naming convention is: const_cast_{0}
 *	- {0} : Datatype.
 *		- Options: i, ptrdiff, bool, const_Element.
 */

#include <stddef.h>
#include <stdbool.h>

// Standard data types ********************************************************************************************** //

/// \brief Cast from `int` to `const int`.
void const_cast_i
	(const int* dest, ///< Destination.
	 const int src    ///< Source.
	);

/// \brief Cast from `int*` to `const int*`.
void const_cast_i1
	(const int* dest, ///< Destination.
	 const int* src,  ///< Source.
	 const int n_src  ///< Number of entries in `src`.
	);

/// \brief Cast from `ptrdiff_t` to `const ptrdiff_t`.
void const_cast_ptrdiff
	(const ptrdiff_t* dest, ///< Destination.
	 const ptrdiff_t src    ///< Source.
	);

/// \brief Cast from `bool` to `const bool`.
void const_cast_bool
	(const bool* dest, ///< Destination.
	 const bool src    ///< Source.
	);

/// \brief Cast from `char` to `const char`.
void const_cast_c
	(const char* dest, ///< Destination.
	 const char src    ///< Source.
	);

/// \brief Cast from `char*` to `const char*const`.
void const_cast_c1
	(const char*const* dest, ///< Destination.
	 const char*const  src   ///< Source.
	);

// Custom data types ************************************************************************************************ //

#endif // DPG__const_cast_h__INCLUDED
