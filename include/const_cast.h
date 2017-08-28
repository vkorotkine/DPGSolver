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
 *			- Options: i, ptrdiff, bool, const_Element.
 */

// Required includes
#include <stddef.h>
#include <stdbool.h>

// Redundant includes

// Forward declarations

// ****************************************************************************************************************** //

// Standard data types ********************************************************************************************** //

/// \brief Cast from `int` to `const int`.
void const_cast_i
	(const int* dest, ///< Destination.
	 const int src    ///< Source.
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

// Custom data types ************************************************************************************************ //

#endif // DPG__const_cast_h__INCLUDED
