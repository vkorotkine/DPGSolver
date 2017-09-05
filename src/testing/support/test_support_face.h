// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_support_face_h__INCLUDED
#define DPG__test_support_face_h__INCLUDED
/**	\file
 *	\brief Provides support functions for testing relating to the \ref Face containers.
 */

#include <stdio.h>

#include "intrusive.h"

/** \brief Constructor for a \ref Face from data in the input file.
 *	\return Standard. */
struct Face* constructor_Face
	(FILE* file,                                       ///< The file from which to read the data.
	 char* line,                                       ///< The first line of the input for the current face.
	 const struct const_Intrusive_List*const elements, ///< \ref Simulation::elements.
	 const struct Intrusive_List*const volumes         ///< \ref Simulation::volumes.
	);

#endif // DPG__test_support_face_h__INCLUDED
