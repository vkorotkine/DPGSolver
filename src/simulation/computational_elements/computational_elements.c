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

#include "macros.h"
#include "definitions_elements.h"
#include "definitions_intrusive.h"

#include "volume.h"
#include "solver_volume.h"
#include "face.h"
#include "solver_face.h"

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

// Interface functions ********************************************************************************************** //

void constructor_computational_element_lists
	(struct Simulation* sim, const struct Mesh*const mesh, const int list_category)
{
	switch (list_category) {
	case IL_BASE:
		sim->volumes = constructor_Volumes(sim,mesh);
		sim->faces   = constructor_Faces(sim,mesh);
	case IL_SOLVER:
		sim->volumes = constructor_Solver_Volumes(sim);
		sim->faces   = constructor_Solver_Faces(sim);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",list_category);
		break;
	}

	if (list_category != IL_BASE) {
		update_computational_element_list_pointers(sim);

		destructor_IL_base(sim->volumes);
		destructor_IL_base(sim->faces);
	}
}

const struct const_Intrusive_List* constructor_derived_Elements (struct Simulation* sim, const int list_name)
{
	const struct const_Intrusive_List* base = sim->elements;
	assert(base->name);

	const struct const_Intrusive_List* list = NULL;
	switch (list_name) {
	case IL_GEOMETRY_ELEMENT:
		list = constructor_Geometry_Elements(sim);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",list_name);
		break;
	}
	update_computational_element_elements(sim);

	return list;
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

// Level 1 ********************************************************************************************************** //

static void update_volume_pointers (struct Intrusive_Link* link)
{
	struct Volume* volume = (struct Volume*) link;

	for (int i = 0; i < NFMAX; ++i) {
	for (int j = 0; j < NSUBFMAX; ++j) {
		struct Face* face = (struct Face*) volume->faces[i][j];
		if (face) {
			const struct Face*const derived = (const struct Face*const)((struct Intrusive_Link*)face)->derived;
			assert(derived != NULL);
			const_cast_Face(&volume->faces[i][j],derived);
		}
	}}
}

static void update_face_pointers (struct Intrusive_Link* link)
{
	struct Face* face = (struct Face*) link;

	for (int i = 0; i < 2; ++i) {
		struct Volume* volume = face->neigh_info[i].volume;
		if (volume) {
			struct Volume* derived = (struct Volume*)((struct Intrusive_Link*)volume)->derived;
			assert(derived != NULL);
			face->neigh_info[i].volume = derived;
		}
	}
}
