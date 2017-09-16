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

#ifndef DPG__allocators_h__INCLUDED
#define DPG__allocators_h__INCLUDED
/** \file
 *	\brief Provides functions using for dynamic memory allocation.
 *
 *	\deprecated Likely not generally useful as multiply dereferenced objects are to be avoided. Variadic arguments
 *	            should be replaced with compound literal inputs if this is used in the future.
 */

#include <stddef.h>

/// \brief Specifies the type of the variable to be used in allocator functions.
enum Variable_Type {
	CHAR_T,
	INT_T,
	DOUBLE_T,
	PTRDIFF_T,
};

/**	\brief Performs allocation using malloc.
 *
 *	The size arguments are passed in correspondence with the level of dereferencing.
 *
 *	Example usage:
 *
 *	The standard allocation:
 *	```
 *	char** a = malloc(N2 * sizeof *a);
 *	for (ptrdiff_t i = 0; i < N2; i++)
 *		a[i] = malloc(N1 * sizeof *a[i]);
 *	```
 *
 *	is equivalent to
 *	```
 *	char** a = mallocator(CHAR_T,2,N1,N2);
 *	```
 */
void* mallocator
	(const enum Variable_Type type, ///< \ref Variable_Type.
	 const int order,               ///< Maximum level of dereferencing of the allocated variable.
	 ...                            ///< Sizes of the allocated variable for each level of dereferencing.
	);

/**	\brief Performs deallocation.
 *
 *	`NULL` can be passed if `order` is 1.
 */
void deallocator
	(void* a,                       ///< Pointer to the variable to be freed.
	 const enum Variable_Type type, ///< Defined in \ref mallocator.
	 const int order,               ///< Defined in \ref mallocator.
	 ...                            ///< Defined in \ref mallocator.
	);

#endif // DPG__allocators_h__INCLUDED
