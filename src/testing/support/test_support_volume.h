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

#ifndef DPG__test_support_volume_h__INCLUDED
#define DPG__test_support_volume_h__INCLUDED
/**	\file
 *	\brief Provides support functions for testing relating to the \ref Volume containers.
 */

#include <stdio.h>

struct Intrusive_List;
struct const_Intrusive_List;

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
