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

#ifndef DPG__computational_elements_h__INCLUDED
#define DPG__computational_elements_h__INCLUDED
/** \file
 *  \brief Provides general functions related to the \ref Volume and \ref Face computational elements and their derived
 *         types.
 */

struct Simulation;
struct Mesh;

/** \brief Constructor for computational element lists.
 *
 *  If a list other than the base list is constructed, its corresponding base list is destructed here.
 */
void constructor_computational_element_lists
	(struct Simulation* sim,       ///< \ref Simulation.
	 const struct Mesh*const mesh, ///< The \ref Mesh (only required for base lists).
	 const int list_category       /**< The computational element list category.
	                                *   Options: see \ref definitions_intrusive.h. */
	);

/** \brief Constructor for a list of elements derived from the base \ref Element list.
 *  \return Standard. */
const struct const_Intrusive_List* constructor_derived_Elements
	(struct Simulation* sim, ///< \ref Simulation
	 const int list_name     ///< The derived \ref Intrusive_List::name.
	);

#endif // DPG__computational_elements_h__INCLUDED
