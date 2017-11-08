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

#include <stddef.h>

struct Simulation;
struct Mesh;

/** \brief Construct derived \ref Volume and \ref Face computational element lists.
 *  \ref Simulation::volumes and \ref Simulation::faces are set to point to the newly created lists.
 */
void constructor_derived_computational_elements
	(struct Simulation* sim,    ///< \ref Simulation.
	 const int derived_category /**< The computational element list category.
	                             *   Options: see \ref definitions_intrusive.h. */
	);

/** \brief Destructor for derived \ref Volume and \ref Face computational element lists.
 *  The appropriate portion of the derived list members are shallow copied to the base list and the derived lists are
 *  then destructed.
 */
void destructor_derived_computational_elements
	(struct Simulation* sim, ///< \ref Simulation.
	 const int base_category /**< The derived computational element list category.
	                          *   Options: see \ref definitions_intrusive.h. */
	);

/** \brief Constructor for a list of derived \ref Element\*s.
 *  \ref Simulation::elements is set to point to the newly created list.
 *  While POINT elements are included as part of the list, derived POINT elements are not constructed (i.e. only the
 *  base POINT element can be used).
 */
void constructor_derived_Elements
	(struct Simulation* sim, ///< \ref Simulation.
	 const int derived_name  ///< The derived \ref Intrusive_List::name.
	);

/** \brief Destructor for a list of derived \ref Element\*s.
 *  The appropriate portion of the derived list elements are shallow copied to the base list and the derived list is
 *  then destructed.
 */
void destructor_derived_Elements
	(struct Simulation* sim, ///< \ref Simulation.
	 const int base_name     ///< The derived \ref Intrusive_List::name.
	);

/** \brief Constructor for a derived \ref Intrusive_Link\* to be inserted in a list.
 *  \return See brief. */
struct Intrusive_Link* constructor_derived_Intrusive_Link
	(struct Intrusive_Link* base, ///< Pointer to the base link.
	 const size_t sizeof_base,    ///< Value of std::sizeof(base).
	 const size_t sizeof_derived  ///< Value of std::sizeof(derived).
	);

#endif // DPG__computational_elements_h__INCLUDED
