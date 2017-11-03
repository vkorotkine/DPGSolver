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

#ifndef DPG__test_support_face_h__INCLUDED
#define DPG__test_support_face_h__INCLUDED
/** \file
 *  \brief Provides support functions for testing relating to the \ref Face containers.
 */

#include <stdio.h>

struct Intrusive_List;
struct const_Intrusive_List;

/** \brief Constructor for a \ref Face from data in the input file.
 *  \return Standard.
 *
 *  If the `volumes` input is `NULL`, members of the face are not set.
 */
struct Face* constructor_Face_file
	(FILE* file,                                       ///< The file from which to read the data.
	 char* line,                                       ///< The first line of the input for the current face.
	 const struct const_Intrusive_List*const elements, ///< \ref Simulation::elements.
	 const struct Intrusive_List*const volumes         ///< \ref Simulation::volumes.
	);

/// \brief Constructor for a \ref Solver_Face from a file.
void constructor_file_Solver_Face
	(FILE* file,           ///< The file from which to read the data.
	 char* line,           ///< The first line of the input for the current face.
	 struct Face* face_ptr ///< Pointer to the face.
	);

#endif // DPG__test_support_face_h__INCLUDED
