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
/** \file
 */

#include "computational_elements.h"

#include "face_solver_adaptive.h"
#include "face_solver_dg.h"
#include "face_solver_dpg.h"
#include "face_solver_opg.h"
#include "volume_solver_adaptive.h"
#include "volume_solver_dg.h"
#include "volume_solver_dpg.h"
#include "volume_solver_opg.h"

#include "element.h"
#include "element_geometry.h"
#include "element_plotting.h"
#include "element_solution.h"
#include "element_adaptation.h"
#include "element_solver.h"
#include "element_solver_dg.h"
#include "element_solver_dpg.h"
#include "element_solver_opg.h"

#include "intrusive.h"
#include "simulation.h"

// Templated functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "computational_elements_T.c"
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "computational_elements_T.c"
#include "undef_templates_type.h"

// Static function declarations ************************************************************************************* //

/** \brief Function pointer to a derived Element destructor function for a element which is part of an
 *         \ref Intrusive_List.
 *  \param element Pointer to the element link in the list.
 */
typedef void (*destructor_derived_Element_fptr)
	(struct Element* element_ptr
	);

/// Container holding information used for \ref constructor_derived_Elements and \ref destructor_derived_Elements.
struct Derived_Elements_Info {
	size_t sizeof_base,    ///< The size of the base element.
	       sizeof_derived; ///< The size of the derived element.

	constructor_derived_Element_fptr constructor_derived_Element; ///< \ref constructor_derived_Element_fptr.
	destructor_derived_Element_fptr  destructor_derived_Element;  ///< \ref destructor_derived_Element_fptr.
};

/** \brief Return a stack-allocated \ref Derived_Elements_Info container for \ref constructor_derived_Elements.
 *  \return See brief. */
static struct Derived_Elements_Info get_c_Derived_Elements_Info
	(const int derived_name,      ///< Defined for \ref constructor_derived_Elements.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Return a stack-allocated \ref Derived_Elements_Info container for \ref destructor_derived_Elements.
 *  \return See brief. */
static struct Derived_Elements_Info get_d_Derived_Elements_Info
	(const int base_name,   ///< Defined for \ref destructor_derived_Elements.
	 const int derived_name ///< Defined for \ref constructor_derived_Elements.
	);

/// \brief Update the pointers to the \ref Element members in the current computational element lists.
static void update_computational_element_elements
	(struct Simulation* sim ///< \ref Simulation.
	);

/** \brief `const` version of \ref constructor_derived_Intrusive_Link.
 *  \return See brief. */
static const struct const_Intrusive_Link* constructor_derived_const_Intrusive_Link
	(const struct const_Intrusive_Link* base, ///< Defined for \ref constructor_derived_Intrusive_Link.
	 const size_t sizeof_base,                ///< Defined for \ref constructor_derived_Intrusive_Link.
	 const size_t sizeof_derived              ///< Defined for \ref constructor_derived_Intrusive_Link.
	);

/** \brief `const` version of \ref constructor_base_Intrusive_Link.
 *  \return See brief. */
static const struct const_Intrusive_Link* constructor_base_const_Intrusive_Link
	(const struct const_Intrusive_Link* derived, ///< Defined for \ref constructor_base_Intrusive_Link.
	 const size_t sizeof_base                    ///< Defined for \ref constructor_base_Intrusive_Link.
	);

// Interface functions ********************************************************************************************** //

void constructor_derived_Elements (struct Simulation* sim, const int derived_name)
{
	struct Derived_Elements_Info de_i = get_c_Derived_Elements_Info(derived_name,sim);

	const struct const_Intrusive_List* base = sim->elements;

	// Reserve memory for the derived element list.
	const struct const_Intrusive_List* elements = constructor_empty_const_IL(derived_name,base); // moved
	for (const struct const_Intrusive_Link* curr = base->first; curr; curr = curr->next)
		push_back_const_IL(elements,
		                   constructor_derived_const_Intrusive_Link(curr,de_i.sizeof_base,de_i.sizeof_derived));
	set_element_pointers((struct Intrusive_List*)elements);

	sim->elements = elements;

	// Perform construction specific to the derived element list.
	for (const struct const_Intrusive_Link* curr = sim->elements->first; curr; ) {
		const struct const_Intrusive_Link* next = curr->next;
		const struct const_Element* element = (struct const_Element*) curr;
		if (element->present)
			de_i.constructor_derived_Element((struct Element*)curr,sim);
		curr = next;
	}

	// Update pointers to the base elements to point to the derived elements.
	update_computational_element_elements(sim);

	// Destruct the base list.
	destructor_const_IL_base(sim->elements);
}

void destructor_derived_Elements (struct Simulation* sim, const int base_name)
{
	struct Derived_Elements_Info de_i = get_d_Derived_Elements_Info(base_name,sim->elements->name);

	// Perform destruction specific to the derived element list.
	for (const struct const_Intrusive_Link* curr = sim->elements->first; curr; ) {
		const struct const_Intrusive_Link* next = curr->next;
		const struct const_Element* element = (struct const_Element*) curr;
		if (element->present)
			de_i.destructor_derived_Element((struct Element*)curr);
		curr = next;
	}

	// Reserve memory for the base element list.
	const struct const_Intrusive_List* elements_prev = sim->elements;
	sim->elements = constructor_empty_const_IL(base_name,NULL); // keep
	for (const struct const_Intrusive_Link* curr = elements_prev->first; curr; curr = curr->next) {
		push_back_const_IL(sim->elements,constructor_base_const_Intrusive_Link(curr,de_i.sizeof_base));
	}

	// Update pointers to the derived elements to point to the base elements.
	update_computational_element_elements(sim);

	// Destruct the derived list.
	destructor_const_IL(elements_prev,true);
}

struct Intrusive_Link* constructor_derived_Intrusive_Link
	(struct Intrusive_Link* base, const size_t sizeof_base, const size_t sizeof_derived)
{
	struct Intrusive_Link* derived = calloc(1,sizeof_derived); // returned
	memcpy(derived,base,sizeof_base); // shallow copy of the base.

	assert(base->derived == NULL);
	base->derived = derived;

	return derived;
}

void constructor_offset_derived_Element
	(constructor_derived_Element_fptr cde, const size_t sizeof_base, const struct Element*const base_e,
	 struct Element*const curr_e, const struct Simulation*const sim)
{
	memcpy(curr_e,base_e,sizeof_base);

	struct Intrusive_List* elements = constructor_empty_IL(IL_INVALID,NULL); // destructed

	const ptrdiff_t curr_offset = BYTE_DIFF(curr_e,base_e);

	for (const struct const_Intrusive_Link* curr = sim->elements->first; curr; curr = curr->next) {
		struct Element* sim_e = (struct Element*) curr;
		push_back_IL(elements,(struct Intrusive_Link*)BYTE_ADD(curr,curr_offset));
		if (sim_e->type == curr_e->type)
			break;
	}
	set_element_pointers(elements);
	destructor_IL(elements,false);

	cde(curr_e,sim);
}

bool list_is_derived_from (const char*const name_desired, const char list_type, const struct Simulation*const sim)
{
	bool is_derived_from = false;
	switch (list_type) {
	case 'e':
		if (strcmp(name_desired,"solver") == 0) {
			switch (sim->elements->name) {
			case IL_ELEMENT_SOLVER:     // fallthrough
			case IL_ELEMENT_SOLVER_DG:  // fallthrough
			case IL_ELEMENT_SOLVER_DPG: // fallthrough
			case IL_ELEMENT_SOLVER_OPG:
				is_derived_from = true;
				break;
			default:
				; // Do nothing
				break;
			}
		} else {
			EXIT_ERROR("Unsupported: %s\n",name_desired);
		}
		break;
	case 'f':
		if (strcmp(name_desired,"solver") == 0) {
			switch (sim->faces->name) {
			case IL_FACE_SOLVER:     // fallthrough
			case IL_FACE_SOLVER_DG:  // fallthrough
			case IL_FACE_SOLVER_DPG: // fallthrough
			case IL_FACE_SOLVER_OPG:
				is_derived_from = true;
				break;
			default:
				; // Do nothing
				break;
			}
		} else {
			EXIT_ERROR("Unsupported: %s\n",name_desired);
		}
		break;
	case 'v':
		if (strcmp(name_desired,"solver") == 0) {
			switch (sim->volumes->name) {
			case IL_VOLUME_SOLVER:     // fallthrough
			case IL_VOLUME_SOLVER_DG:  // fallthrough
			case IL_VOLUME_SOLVER_DPG: // fallthrough
			case IL_VOLUME_SOLVER_OPG:
				is_derived_from = true;
				break;
			default:
				; // Do nothing
				break;
			}
		} else {
			EXIT_ERROR("Unsupported: %s\n",name_desired);
		}
		break;
	default:
		EXIT_ERROR("Unsupported: %c\n",list_type);
		break;
	}
	return is_derived_from;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct Derived_Elements_Info get_c_Derived_Elements_Info (const int derived_name, const struct Simulation* sim)
{
	struct Derived_Elements_Info de_info;

	switch (derived_name) {
	case IL_ELEMENT_GEOMETRY:
		assert(sim->elements->name == IL_ELEMENT);
		de_info.sizeof_base    = sizeof(struct Element);
		de_info.sizeof_derived = sizeof(struct Geometry_Element);
		de_info.constructor_derived_Element = constructor_derived_Geometry_Element;
		break;
	case IL_ELEMENT_PLOTTING:
		assert(sim->elements->name == IL_ELEMENT);
		de_info.sizeof_base    = sizeof(struct Element);
		de_info.sizeof_derived = sizeof(struct Plotting_Element);
		de_info.constructor_derived_Element = constructor_derived_Plotting_Element;
		break;
	case IL_ELEMENT_SOLUTION:
		assert(sim->elements->name == IL_ELEMENT);
		de_info.sizeof_base    = sizeof(struct Element);
		de_info.sizeof_derived = sizeof(struct Solution_Element);
		de_info.constructor_derived_Element = constructor_derived_Solution_Element;
		break;
	case IL_ELEMENT_ADAPTATION:
		assert(sim->elements->name == IL_ELEMENT);
		de_info.sizeof_base    = sizeof(struct Element);
		de_info.sizeof_derived = sizeof(struct Adaptation_Element);
		de_info.constructor_derived_Element = constructor_derived_Adaptation_Element;
		break;
	case IL_ELEMENT_SOLVER:
		assert(sim->elements->name == IL_ELEMENT);
		de_info.sizeof_base    = sizeof(struct Element);
		de_info.sizeof_derived = sizeof(struct Solver_Element);
		de_info.constructor_derived_Element = constructor_derived_Solver_Element;
		break;
	case IL_ELEMENT_SOLVER_DG:
		assert(sim->elements->name == IL_ELEMENT_SOLVER);
		de_info.sizeof_base    = sizeof(struct Solver_Element);
		de_info.sizeof_derived = sizeof(struct DG_Solver_Element);
		de_info.constructor_derived_Element = constructor_derived_DG_Solver_Element;
		break;
	case IL_ELEMENT_SOLVER_DPG:
		assert(sim->elements->name == IL_ELEMENT_SOLVER);
		de_info.sizeof_base    = sizeof(struct Solver_Element);
		de_info.sizeof_derived = sizeof(struct DPG_Solver_Element);
		de_info.constructor_derived_Element = constructor_derived_DPG_Solver_Element;
		break;
	case IL_ELEMENT_SOLVER_OPG:
		assert(sim->elements->name == IL_ELEMENT_SOLVER);
		de_info.sizeof_base    = sizeof(struct Solver_Element);
		de_info.sizeof_derived = sizeof(struct OPG_Solver_Element);
		de_info.constructor_derived_Element = constructor_derived_OPG_Solver_Element;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",derived_name);
		break;
	}
	return de_info;
}

static struct Derived_Elements_Info get_d_Derived_Elements_Info (const int base_name, const int derived_name)
{
	struct Derived_Elements_Info de_info;

	switch (base_name) {
	case IL_ELEMENT:
		de_info.sizeof_base = sizeof(struct Element);
		break;
	case IL_ELEMENT_SOLUTION:
		de_info.sizeof_base = sizeof(struct Solution_Element);
		break;
	case IL_ELEMENT_SOLVER:
		de_info.sizeof_base = sizeof(struct Solver_Element);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",base_name);
		break;
	}

	switch (derived_name) {
	case IL_ELEMENT_GEOMETRY:
		assert(base_name == IL_ELEMENT);
		de_info.destructor_derived_Element = destructor_derived_Geometry_Element;
		break;
	case IL_ELEMENT_PLOTTING:
		assert(base_name == IL_ELEMENT);
		de_info.destructor_derived_Element = destructor_derived_Plotting_Element;
		break;
	case IL_ELEMENT_SOLUTION:
		assert(base_name == IL_ELEMENT);
		de_info.destructor_derived_Element = destructor_derived_Solution_Element;
		break;
	case IL_ELEMENT_ADAPTATION:
		assert(base_name == IL_ELEMENT);
		de_info.destructor_derived_Element = destructor_derived_Adaptation_Element;
		break;
	case IL_ELEMENT_SOLVER:
		assert(base_name == IL_ELEMENT);
		de_info.destructor_derived_Element = destructor_derived_Solver_Element;
		break;
	case IL_ELEMENT_SOLVER_DG:
		assert(base_name == IL_ELEMENT_SOLVER);
		de_info.destructor_derived_Element = destructor_derived_DG_Solver_Element;
		break;
	case IL_ELEMENT_SOLVER_DPG:
		assert(base_name == IL_ELEMENT_SOLVER);
		de_info.destructor_derived_Element = destructor_derived_DPG_Solver_Element;
		break;
	case IL_ELEMENT_SOLVER_OPG:
		assert(base_name == IL_ELEMENT_SOLVER);
		de_info.destructor_derived_Element = destructor_derived_OPG_Solver_Element;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",derived_name);
		break;
	}
	return de_info;
}

static void update_computational_element_elements (struct Simulation* sim)
{
	const struct const_Intrusive_List* elements = sim->elements;

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Volume* volume = (struct Volume*) curr;
		const_cast_const_Element(&volume->element,get_element_by_type(elements,volume->element->type));
	}

	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Face* face = (struct Face*) curr;
		const_cast_const_Element(&face->element,get_element_by_type(elements,face->element->type));
	}
}

static const struct const_Intrusive_Link* constructor_derived_const_Intrusive_Link
	(const struct const_Intrusive_Link* base, const size_t sizeof_base, const size_t sizeof_derived)
{
	return (const struct const_Intrusive_Link*)
		constructor_derived_Intrusive_Link((struct Intrusive_Link*)base,sizeof_base,sizeof_derived);
}

static const struct const_Intrusive_Link* constructor_base_const_Intrusive_Link
	(const struct const_Intrusive_Link* derived, const size_t sizeof_base)
{
	return (const struct const_Intrusive_Link*)
		constructor_base_Intrusive_Link((struct Intrusive_Link*)derived,sizeof_base);
}
