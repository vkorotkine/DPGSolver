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
 *  \brief Provides real general functions related to the \ref Volume and \ref Face computational elements
 *         and their derived types.
 */

#include "def_templates_type_d.h"
#include "computational_elements_T.h"
#include "undef_templates_type.h"

#include <stddef.h>
#include <stdbool.h>

struct Simulation;
struct Element;

/** \brief Function pointer to a derived Element constructor function for a element which is part of an
 *         \ref Intrusive_List.
 *  \param element_ptr Pointer to the element link in the list.
 *  \param sim         \ref Simulation.
 */
typedef void (*constructor_derived_Element_fptr)
	(struct Element* element_ptr,
	 const struct Simulation* sim
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

/** \brief Constructor for a derived element as part of a container and whose memory is offset from the start of the
 *         memory of the container.
 *
 *  The usage of this function is motivated by the fact that the cast to the intrusive base element would no longer work
 *  properly when this memory offset is present.
 *
 *  The base element is (generally redundantly) shallow-copied into the derived element such that the derived element
 *  may be used normally as if it were constructed in the standard manner.
 *
 *  This function is used, for example, to construct elements which are part of other elements.
 */
void constructor_offset_derived_Element
	(constructor_derived_Element_fptr cde, /**< \ref constructor_derived_Element_fptr for the element type under
	                                        *   consideration. */
	 const size_t sizeof_base,             ///< The size of the base element.
	 const struct Element*const base_e,    ///< Pointer to the base element.
	 struct Element*const curr_e,          ///< Pointer to the current element.
	 const struct Simulation*const sim     ///< Defined for \ref constructor_derived_Element_fptr.
	);

/** \brief Check if the list of input type is derived from the desired type.
 *  \return `true` if yes; `false` otherwise. */
bool list_is_derived_from
	(const char*const name_desired,    ///< Name of the desired type. Options: "solver".
	 const char list_type,             ///< Type of the list. Options: 'e'lement, 'f'ace, 'v'olume.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

#endif // DPG__computational_elements_h__INCLUDED
