// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_support_volume_h__INCLUDED
#define DPG__test_support_volume_h__INCLUDED
/**	\file
 *	\brief Provides support functions for testing relating to the \ref Volume containers.
 */

#include <stdio.h>

#include "intrusive.h"

/** \brief Constructor for a \ref Volume from data in the input file.
 *	\return Standard. */
struct Volume* constructor_Volume
	(FILE* file,                                      ///< The file from which to read the data.
	 char* line,                                      ///< The first line of the input for the current volume.
	 const struct const_Intrusive_List*const elements ///< \ref Simulation::elements.
	);

/**	\brief See return.
 *	\return Pointer to the \ref Volume with the input `index`.
 *
 *	\warning Uses an inefficient linear search and should not be used outside of test functions.
 */
struct Volume* get_volume_by_index
	(const struct Intrusive_List*const volumes, ///< The list of volumes.
	 const int index                            ///< The index.
	);

#endif // DPG__test_support_volume_h__INCLUDED
