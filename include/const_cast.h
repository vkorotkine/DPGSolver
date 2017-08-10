// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__const_cast_h__INCLUDED
#define DPG__const_cast_h__INCLUDED
/**	\file
 *	\brief Provides functions for casting to const lvalues for standard datatypes.
 *
 *	Variables with standard datatypes with higher levels of dereferencing than 0 should be placed in Multiarray
 *	containers. Const cast functions are thus not provided for higher levels of dereferencing.
 *
 *	The function naming convention is: const_cast_{0}
 *		- {0} : Datatype.
 *			- Options: ui.
 */

#include <stddef.h>

/// \brief Cast from `unsigned int` to `const unsigned int`.
void const_cast_ui
	(const unsigned int* dest, ///< Destination.
	 unsigned int src          ///< Source.
	);

/// \brief Cast from `size_t` to `const size_t`.
void const_cast_st
	(const size_t* dest, ///< Destination.
	 const size_t src    ///< Source.
	);

#endif // DPG__const_cast_h__INCLUDED
