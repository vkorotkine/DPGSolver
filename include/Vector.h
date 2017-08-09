// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__Vector_h__INCLUDED
#define DPG__Vector_h__INCLUDED
/**	\file
 *	\brief Provides Vector_\* containers and related functions.
 *
 *	Potentially relevant comments may be found in containers.h.
 *
 *	Vectors are 1D Multiarrays.
 */

#include <stddef.h>
#include <stdbool.h>

/// \brief Vector (`unsigned int`).
struct Vector_ui {
	size_t extents[1];

	bool    owns_data;
	unsigned int* data;
};

// Constructor/Destructor functions ********************************************************************************* //

/// \brief Constructs an empty \ref Vector_ui.
struct Vector_ui* constructor_empty_Vector_ui
	(const size_t n_rows ///< The value of extents[0].
	);

// Helper functions ************************************************************************************************* //

/** \brief Reorder a \ref Vector_ui based on the provided ordering.
 *	\warning This is not currently done in place and should consequently not be used for large Vectors if limited memory
 *	         usage is a concern.
 */
void reorder_Vector_ui
	(struct Vector_ui*const a,
	 unsigned int*const ordering
	);

// Printing functions *********************************************************************************************** //

/// \brief Print a \ref Vector_ui to the terminal.
void print_Vector_ui
	(const struct Vector_ui*const a ///< Standard.
	);

#endif // DPG__Vector_h__INCLUDED
