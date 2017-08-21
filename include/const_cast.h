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
 *			- Options: ui, st, bool, const_Element.
 */

#include <stddef.h>
#include <stdbool.h>

#include "Element.h"

// Standard data types ********************************************************************************************** //

/// \brief Cast from `unsigned int` to `const unsigned int`.
void const_cast_ui
	(const unsigned int* dest, ///< Destination.
	 const unsigned int src    ///< Source.
	);

/// \brief Cast from `size_t` to `const size_t`.
void const_cast_st
	(const size_t* dest, ///< Destination.
	 const size_t src    ///< Source.
	);

/// \brief Cast from `bool` to `const bool`.
void const_cast_bool
	(const bool* dest, ///< Destination.
	 const bool src    ///< Source.
	);

// Custom data types ************************************************************************************************ //

/// \brief Cast from \ref const_Element\* to `const` \ref const_Element `*const`.
void const_cast_const_Element
	(const struct const_Element*const* dest, ///< Destination.
	 const struct const_Element*const src    ///< Source.
	);

#endif // DPG__const_cast_h__INCLUDED
