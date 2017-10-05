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

/** \brief Update pointers to computational elements in the derived lists.
 *
 *  Every time a new derived list is created and the accompanying base list is destructed, any pointers to computational
 *  elements in the base list become invalid. This function updates these pointers such that all pointers to base list
 *  links are replaced with pointers to derived list links.
 */
void update_computational_element_list_pointers
	(const struct Simulation* sim /// \ref Simulation.
	);

#endif // DPG__computational_elements_h__INCLUDED
