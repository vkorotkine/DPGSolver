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

#ifndef DPG__test_support_intrusive_h__INCLUDED
#define DPG__test_support_intrusive_h__INCLUDED
/** \file
 *  \brief Provides support functions for testing relating to containers defined in \ref intrusive.h.
 */

struct Intrusive_List;
struct const_Intrusive_List;

/** \brief Constructor for an intrusive list of type specified by the `list_name`.
 *  \return Standard. */
struct Intrusive_List* constructor_file_name_IL
	(const char*const list_name,                       ///< The name of the type of container to form the links.
	 const char*const file_name,                       ///< The name of the file from which to read the data.
	 const struct const_Intrusive_List*const elements, ///< \ref Simulation::elements.
	 const struct Intrusive_List*const volumes         ///< \ref Simulation::volumes.
	);

/// \brief Constructor for supported derived \ref Face members from the input file in the input list.
void constructor_file_name_derived_Faces
	(struct Intrusive_List* faces, ///< The list of derived faces.
	 const char*const file_name    ///< The file name.
	);

#endif // DPG__test_support_intrusive_h__INCLUDED
