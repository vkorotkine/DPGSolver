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

#include <assert.h>
#include <string.h>

#include "macros.h"
#include "definitions_elements.h"
#include "definitions_intrusive.h"

#include "volume.h"
#include "solver_volume.h"
#include "face.h"
#include "solver_face.h"

#include "element.h"
#include "geometry_element.h"

#include "intrusive.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

/** \brief Update pointers to computational elements in the derived lists.
 *
 *  Every time a new derived list is created and the accompanying base list is destructed, any pointers to computational
 *  elements in the base list become invalid. This function updates these pointers such that all pointers to base list
 *  links are replaced with pointers to derived list links.
 */
static void update_computational_element_list_pointers
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Update the pointers to the \ref Element members in the current computational element lists.
static void update_computational_element_elements
	(struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Constructor for a derived \ref Intrusive_Link\* to be inserted in a list.
 *  \return See brief. */
struct Intrusive_Link* constructor_derived_Intrusive_Link
	(struct Intrusive_Link* base, ///< Pointer to the base link.
	 const size_t sizeof_base,    ///< Value of std::sizeof(base).
	 const size_t sizeof_derived  ///< Value of std::sizeof(derived).
	);

/** \brief `const` version of \ref constructor_derived_Intrusive_Link.
 *  \return See brief. */
const struct const_Intrusive_Link* constructor_derived_const_Intrusive_Link
	(const struct const_Intrusive_Link* base, ///< Defined for \ref constructor_derived_Intrusive_Link.
	 const size_t sizeof_base,                ///< Defined for \ref constructor_derived_Intrusive_Link.
	 const size_t sizeof_derived              ///< Defined for \ref constructor_derived_Intrusive_Link.
	);

/** \brief Construts a base \ref Element by copying the input amount of memory from the derived element.
 *  \return The \ref Intrusive_Link\* to be inserted in the list. */
static const struct const_Intrusive_Link* constructor_base_Element
	(struct const_Element* derived, ///< Pointer to the current derived \ref Element.
	 const size_t base_size         ///< The size (in bytes) of the memory to copy from the derived element.
	);

// Interface functions ********************************************************************************************** //

void construct_derived_computational_elements (struct Simulation* sim, const int list_category)
{
	// Reserve memory for the derived lists.
	int list_name[2]         = { 0, 0, };
	size_t sizeof_base[2]    = { 0, 0, },
	       sizeof_derived[2] = { 0, 0, };
	switch (list_category) {
	case IL_SOLVER:
		list_name[0] = IL_SOLVER_VOLUME;
		list_name[1] = IL_SOLVER_FACE;
		sizeof_base[0] = sizeof(struct Volume);
		sizeof_base[1] = sizeof(struct Face);
		sizeof_derived[0] = sizeof(struct Solver_Volume);
		sizeof_derived[1] = sizeof(struct Solver_Face);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",list_category);
		break;
	}

	struct Intrusive_List* base[2]    = { sim->volumes, sim->faces, };
	struct Intrusive_List* derived[2] = { NULL, NULL, };

	for (int i = 0; i < 2; ++i) {
		derived[i] = constructor_empty_IL(list_name[i],base[i]);
		for (struct Intrusive_Link* curr = base[i]->first; curr; curr = curr->next)
			push_back_IL(derived[i],constructor_derived_Intrusive_Link(curr,sizeof_base[i],sizeof_derived[i]));
	}

	sim->volumes = derived[0];
	sim->faces   = derived[1];

	// Perform construction specific to the derived lists.
	switch (list_category) {
	case IL_SOLVER:
		construct_Solver_Volumes(sim);
		construct_Solver_Faces(sim);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",list_category);
		break;
	}

	// Update pointers to the base computational elements to point to the derived computational elements.
	update_computational_element_list_pointers(sim);

	// Destruct the base lists.
	destructor_IL_base(sim->volumes);
	destructor_IL_base(sim->faces);
}

void constructor_derived_Elements (struct Simulation* sim, const int list_name)
{
	// Reserve memory for the derived element list.
	size_t sizeof_base    = 0,
	       sizeof_derived = 0;
	switch (list_name) {
	case IL_GEOMETRY_ELEMENT:
		sizeof_base    = sizeof(struct Element);
		sizeof_derived = sizeof(struct Geometry_Element);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",list_name);
		break;
	}

	const struct const_Intrusive_List* base = sim->elements;

	const struct const_Intrusive_List* elements = constructor_empty_const_IL(list_name,base); // moved
	for (const struct const_Intrusive_Link* curr = base->first; curr; curr = curr->next)
		push_back_const_IL(elements,constructor_derived_const_Intrusive_Link(curr,sizeof_base,sizeof_derived));
	set_tp_sub_elements((struct Intrusive_List*)elements);

	sim->elements = elements;

	// Perform construction specific to the derived element list.
	switch (list_name) {
	case IL_GEOMETRY_ELEMENT: constructor_Geometry_Elements(sim); break;
	default:
		EXIT_ERROR("Unsupported: %d\n",list_name);
		break;
	}

	// Update pointers to the base elements to point to the derived elements.
	update_computational_element_elements(sim);

	// Destruct the base list.
	destructor_const_IL_base(sim->elements);
}

void destructor_derived_Elements (struct Simulation* sim, const int base_name)
{
	const int curr_name = sim->elements->name;
	switch (curr_name) {
		case IL_GEOMETRY_ELEMENT: destructor_Geometry_Elements(sim->elements); break;
		default:
			EXIT_ERROR("Unsupported: %d\n",curr_name);
			break;
	}

	size_t base_size = 0;
	switch (base_name) {
		case IL_ELEMENT: base_size = sizeof(struct Element); break;
		default:
			EXIT_ERROR("Unsupported: %d\n",base_name);
			break;
	}

	const struct const_Intrusive_List* elements = constructor_empty_const_IL(base_name,NULL); // moved
	for (const struct const_Intrusive_Link* curr = sim->elements->first; curr; curr = curr->next)
		push_back_const_IL(elements,constructor_base_Element((struct const_Element*)curr,base_size));

	destructor_const_IL(sim->elements);
	sim->elements = elements;
	update_computational_element_elements(sim);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Update pointers in \ref Volume to the derived list links.
static void update_volume_pointers
	(struct Intrusive_Link* link ///< The current link.
	);

/// \brief Update pointers in \ref Face to the derived list links.
static void update_face_pointers
	(struct Intrusive_Link* link ///< The current link.
	);

static void update_computational_element_list_pointers (const struct Simulation* sim)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next)
		update_volume_pointers(curr);

	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next)
		update_face_pointers(curr);
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

struct Intrusive_Link* constructor_derived_Intrusive_Link
	(struct Intrusive_Link* base, const size_t sizeof_base, const size_t sizeof_derived)
{
	struct Intrusive_Link* derived = calloc(1,sizeof_derived); // returned
	memcpy(derived,base,sizeof_base); // shallow copy of the base.

	assert(base->derived == NULL);
	base->derived = derived;

	return derived;
}

const struct const_Intrusive_Link* constructor_derived_const_Intrusive_Link
	(const struct const_Intrusive_Link* base, const size_t sizeof_base, const size_t sizeof_derived)
{
	return (const struct const_Intrusive_Link*)
		constructor_derived_Intrusive_Link((struct Intrusive_Link*)base,sizeof_base,sizeof_derived);
}

static const struct const_Intrusive_Link* constructor_base_Element
	(struct const_Element* derived, const size_t base_size)
{
	void* element = malloc(base_size);
	memcpy(element,derived,base_size);

	return (const struct const_Intrusive_Link*) element;
}

// Level 1 ********************************************************************************************************** //

static void update_volume_pointers (struct Intrusive_Link* link)
{
	struct Volume* volume = (struct Volume*) link;

	for (int i = 0; i < NFMAX; ++i) {
	for (int j = 0; j < NSUBFMAX; ++j) {
		struct Face* face = (struct Face*) volume->faces[i][j];
		if (!face)
			continue;

		const struct Face*const derived = (const struct Face*const)((struct Intrusive_Link*)face)->derived;
		assert(derived != NULL);
		const_cast_Face(&volume->faces[i][j],derived);
	}}
}

static void update_face_pointers (struct Intrusive_Link* link)
{
	struct Face* face = (struct Face*) link;

	for (int i = 0; i < 2; ++i) {
		struct Volume* volume = face->neigh_info[i].volume;
		if (!volume)
			continue;

		struct Volume* derived = (struct Volume*)((struct Intrusive_Link*)volume)->derived;
		assert(derived != NULL);
		face->neigh_info[i].volume = derived;
	}
}
