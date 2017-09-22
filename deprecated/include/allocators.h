// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__allocators_h__INCLUDED
#define DPG__allocators_h__INCLUDED
/// \file

#include <stddef.h>

/// \brief Specifies the type of the variable to be used in allocator functions.
enum TYPE {
	CHAR_T,
	UINT_T,
	DOUBLE_T,
	SIZE_T_T,
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
 *	for (size_t i = 0; i < N2; i++)
 *		a[i] = malloc(N1 * sizeof *a[i]);
 *	```
 *
 *	is equivalent to
 *	```
 *	char** a = mallocator(CHAR_T,2,N1,N2);
 *	```
 */
void* mallocator
	(const enum TYPE type, ///< Type of the variable.
	 const size_t order,   ///< Maximum level of dereferencing of the allocated variable.
	 ...                   ///< Sizes of the allocated variable for each level of dereferencing.
	);

/**	\brief Performs deallocation.
 *
 *	`NULL` can be passed if `order` is 1.
 */
void deallocator
	(void* a,              ///< Pointer to the variable to be freed.
	 const enum TYPE type, ///< Defined in \ref mallocator.
	 const size_t order,   ///< Defined in \ref mallocator.
	 ...                   ///< Defined in \ref mallocator.
	);

#endif // DPG__allocators_h__INCLUDED
