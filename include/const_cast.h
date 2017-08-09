// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__const_cast_h__INCLUDED
#define DPG__const_cast_h__INCLUDED
/**	\file
 *	\brief Provides functions for casting to const lvalues for standard datatypes.
 *
 *	The function naming convention is: const_cast_{0}_{1}
 *		- {0} : Datatype.
 *			- Options: ui
 *		- {1} : Level of dereferencing.
 */

#include <stddef.h>

/// \brief Cast from `unsigned int*` to `const unsigned int*const`.
void const_cast_ui_1
	(const unsigned int*const* dest, ///< Destination.
	 unsigned int* src               ///< Source.
	);

/// \brief Cast from `unsigned int**` to `const unsigned int*const*const`.
void const_cast_ui_2
	(const unsigned int*const*const* dest, ///< Destination.
	 unsigned int** src                    ///< Source.
	);

/// \brief Cast from `size_t` to `const size_t`.
void const_cast_st_0
	(const size_t* dest, ///< Destination.
	 const size_t src    ///< Source.
	);

#endif // DPG__const_cast_h__INCLUDED
