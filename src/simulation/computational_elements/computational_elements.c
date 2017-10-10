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
#include "element_plotting.h"

#include "intrusive.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

/** \brief Function pointer to a derived Volume constructor function for a volume which is part of an
 *         \ref Intrusive_List.
 *  \param volume_ptr Pointer to the volume link in the list.
 *  \param sim        \ref Simulation.
 */
typedef void (*constructor_derived_Volume_fptr)
	(struct Volume* volume_ptr,
	 const struct Simulation* sim
	);

/** \brief Function pointer to a derived Volume destructor function for a volume which is part of an
 *         \ref Intrusive_List.
 *  \param volume Pointer to the volume link in the list.
 */
typedef void (*destructor_derived_Volume_fptr)
	(struct Volume* volume_ptr
	);

/** \brief Function pointer to a derived Face constructor function for a face which is part of an \ref Intrusive_List.
 *  \param face_ptr Pointer to the face link in the list.
 *  \param sim      \ref Simulation.
 */
typedef void (*constructor_derived_Face_fptr)
	(struct Face* face_ptr,
	 const struct Simulation* sim
	);

/** \brief Function pointer to a derived Face destructor function for a face which is part of an \ref Intrusive_List.
 *  \param face Pointer to the face link in the list.
 */
typedef void (*destructor_derived_Face_fptr)
	(struct Face* face_ptr
	);

/** \brief Function pointer to a derived Element constructor function for a element which is part of an
 *         \ref Intrusive_List.
 *  \param element_ptr Pointer to the element link in the list.
 *  \param sim         \ref Simulation.
 */
typedef void (*constructor_derived_Element_fptr)
	(struct Element* element_ptr,
	 const struct Simulation* sim
	);

/** \brief Function pointer to a derived Element destructor function for a element which is part of an
 *         \ref Intrusive_List.
 *  \param element Pointer to the element link in the list.
 */
typedef void (*destructor_derived_Element_fptr)
	(struct Element* element_ptr
	);

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
static struct Intrusive_Link* constructor_derived_Intrusive_Link
	(struct Intrusive_Link* base, ///< Pointer to the base link.
	 const size_t sizeof_base,    ///< Value of std::sizeof(base).
	 const size_t sizeof_derived  ///< Value of std::sizeof(derived).
	);

/** \brief `const` version of \ref constructor_derived_Intrusive_Link.
 *  \return See brief. */
static const struct const_Intrusive_Link* constructor_derived_const_Intrusive_Link
	(const struct const_Intrusive_Link* base, ///< Defined for \ref constructor_derived_Intrusive_Link.
	 const size_t sizeof_base,                ///< Defined for \ref constructor_derived_Intrusive_Link.
	 const size_t sizeof_derived              ///< Defined for \ref constructor_derived_Intrusive_Link.
	);

/** \brief Constructor for a base \ref Intrusive_Link\* to be inserted in a list.
 *  \return See brief. */
static struct Intrusive_Link* constructor_base_Intrusive_Link
	(struct Intrusive_Link* derived, ///< Pointer to the derived link.
	 const size_t sizeof_base        ///< Value of std::sizeof(base).
	);

/** \brief `const` version of \ref constructor_base_Intrusive_Link.
 *  \return See brief. */
static const struct const_Intrusive_Link* constructor_base_const_Intrusive_Link
	(const struct const_Intrusive_Link* derived, ///< Defined for \ref constructor_base_Intrusive_Link.
	 const size_t sizeof_base                    ///< Defined for \ref constructor_base_Intrusive_Link.
	);

// Interface functions ********************************************************************************************** //

void constructor_derived_computational_elements (struct Simulation* sim, const int derived_category)
{
	// Set parameters
// Make external function here after usage is set.
	int list_name[2]         = { 0, 0, };
	size_t sizeof_base[2]    = { 0, 0, },
	       sizeof_derived[2] = { 0, 0, };
	constructor_derived_Volume_fptr constructor_derived_Volume = NULL;
	constructor_derived_Face_fptr   constructor_derived_Face   = NULL;
	switch (derived_category) {
	case IL_SOLVER:
		list_name[0] = IL_SOLVER_VOLUME;
		list_name[1] = IL_SOLVER_FACE;
		sizeof_base[0] = sizeof(struct Volume);
		sizeof_base[1] = sizeof(struct Face);
		sizeof_derived[0] = sizeof(struct Solver_Volume);
		sizeof_derived[1] = sizeof(struct Solver_Face);
		constructor_derived_Volume = constructor_derived_Solver_Volume;
		constructor_derived_Face   = constructor_derived_Solver_Face;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",derived_category);
		break;
	}

	struct Intrusive_List* base[2]    = { sim->volumes, sim->faces, };
	struct Intrusive_List* derived[2] = { NULL, NULL, };

	// Reserve memory for the derived lists.
	for (int i = 0; i < 2; ++i) {
		derived[i] = constructor_empty_IL(list_name[i],base[i]);
		for (struct Intrusive_Link* curr = base[i]->first; curr; curr = curr->next)
			push_back_IL(derived[i],constructor_derived_Intrusive_Link(curr,sizeof_base[i],sizeof_derived[i]));
	}

	sim->volumes = derived[0];
	sim->faces   = derived[1];

	// Perform construction specific to the derived lists.
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; ) {
		struct Intrusive_Link* next = curr->next;
		constructor_derived_Volume((struct Volume*)curr,sim);
		curr = next;
	}

	for (struct Intrusive_Link* curr = sim->faces->first; curr; ) {
		struct Intrusive_Link* next = curr->next;
		constructor_derived_Face((struct Face*)curr,sim);
		curr = next;
	}

	// Update pointers to the base computational elements to point to the derived computational elements.
	update_computational_element_list_pointers(sim);

	// Destruct the base lists.
	destructor_IL_base(sim->volumes);
	destructor_IL_base(sim->faces);
}

void destructor_derived_computational_elements (struct Simulation* sim, const int derived_category)
{
	// Set parameters
// Make external function here after usage is set.
	int base_name[2]         = { 0, 0, };
	size_t sizeof_base[2]    = { 0, 0, };
	destructor_derived_Volume_fptr destructor_derived_Volume = NULL;
	destructor_derived_Face_fptr   destructor_derived_Face   = NULL;
	switch (derived_category) {
	case IL_SOLVER:
		base_name[0] = IL_VOLUME;
		base_name[1] = IL_FACE;
		sizeof_base[0] = sizeof(struct Volume);
		sizeof_base[1] = sizeof(struct Face);
		destructor_derived_Volume = destructor_derived_Solver_Volume;
		destructor_derived_Face   = destructor_derived_Solver_Face;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",derived_category);
		break;
	}

	// Perform destruction specific to the derived lists.
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; ) {
		struct Intrusive_Link* next = curr->next;
		destructor_derived_Volume((struct Volume*)curr);
		curr = next;
	}

	for (struct Intrusive_Link* curr = sim->faces->first; curr; ) {
		struct Intrusive_Link* next = curr->next;
		destructor_derived_Face((struct Face*)curr);
		curr = next;
	}

	// Reserve memory for the base lists.
	int ind_list = -1;

	++ind_list;
	struct Intrusive_List* volumes_prev = sim->volumes;
	sim->volumes = constructor_empty_IL(base_name[ind_list],NULL); // keep
	for (struct Intrusive_Link* curr = volumes_prev->first; curr; curr = curr->next)
		push_back_IL(sim->volumes,constructor_base_Intrusive_Link(curr,sizeof_base[ind_list]));

	++ind_list;
	struct Intrusive_List* faces_prev = sim->faces;
	sim->faces = constructor_empty_IL(base_name[ind_list],NULL); // keep
	for (struct Intrusive_Link* curr = faces_prev->first; curr; curr = curr->next)
		push_back_IL(sim->faces,constructor_base_Intrusive_Link(curr,sizeof_base[ind_list]));

	// Update pointers to the base computational elements to point to the derived computational elements.
	update_computational_element_list_pointers(sim);

	// Destruct the derived lists.
	destructor_IL(volumes_prev);
	destructor_IL(faces_prev);
}

void constructor_derived_Elements (struct Simulation* sim, const int derived_name)
{
	// Set parameters
// Make external function here after usage is set.
	size_t sizeof_base    = 0,
	       sizeof_derived = 0;
	constructor_derived_Element_fptr constructor_derived_Element = NULL;
	switch (derived_name) {
	case IL_GEOMETRY_ELEMENT:
		assert(sizeof(struct Geometry_Element) == sizeof(struct const_Geometry_Element));
		sizeof_base    = sizeof(struct Element);
		sizeof_derived = sizeof(struct Geometry_Element);
		constructor_derived_Element = constructor_derived_Geometry_Element;
		break;
	case IL_PLOTTING_ELEMENT:
		assert(sizeof(struct Plotting_Element) == sizeof(struct const_Plotting_Element));
		sizeof_base    = sizeof(struct Element);
		sizeof_derived = sizeof(struct Plotting_Element);
		constructor_derived_Element = constructor_derived_Plotting_Element;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",derived_name);
		break;
	}

	const struct const_Intrusive_List* base = sim->elements;

	// Reserve memory for the derived element list.
	const struct const_Intrusive_List* elements = constructor_empty_const_IL(derived_name,base); // moved
	for (const struct const_Intrusive_Link* curr = base->first; curr; curr = curr->next)
		push_back_const_IL(elements,constructor_derived_const_Intrusive_Link(curr,sizeof_base,sizeof_derived));
	set_tp_sub_elements((struct Intrusive_List*)elements);

	sim->elements = elements;

	// Perform construction specific to the derived element list.
	for (const struct const_Intrusive_Link* curr = sim->elements->first; curr; ) {
		const struct const_Intrusive_Link* next = curr->next;
		if (((struct Element*)curr)->present)
			constructor_derived_Element((struct Element*)curr,sim);
		curr = next;
	}

	// Update pointers to the base elements to point to the derived elements.
	update_computational_element_elements(sim);

	// Destruct the base list.
	destructor_const_IL_base(sim->elements);
}

void destructor_derived_Elements (struct Simulation* sim, const int derived_name)
{
	// Set parameters
// Make external function here after usage is set.
	int base_name      = -1;
	size_t sizeof_base = 0;
	destructor_derived_Element_fptr destructor_derived_Element = NULL;

	const int curr_name = sim->elements->name;
	assert(curr_name == derived_name);
	switch (curr_name) {
		case IL_GEOMETRY_ELEMENT:
			base_name = IL_ELEMENT;
			sizeof_base = sizeof(struct Element);
			destructor_derived_Element = destructor_derived_Geometry_Element;
			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",curr_name);
			break;
	}

	// Perform destruction specific to the derived element list.
	for (const struct const_Intrusive_Link* curr = sim->elements->first; curr; ) {
		const struct const_Intrusive_Link* next = curr->next;
		if (((struct Element*)curr)->present)
			destructor_derived_Element((struct Element*)curr);
		curr = next;
	}

	// Reserve memory for the base element list.
	const struct const_Intrusive_List* elements_prev = sim->elements;
	sim->elements = constructor_empty_const_IL(base_name,NULL); // keep
	for (const struct const_Intrusive_Link* curr = elements_prev->first; curr; curr = curr->next) {
		push_back_const_IL(sim->elements,constructor_base_const_Intrusive_Link(curr,sizeof_base));
	}

	// Update pointers to the derived elements to point to the base elements.
	update_computational_element_elements(sim);

	// Destruct the derived list.
	destructor_const_IL(elements_prev);
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

static struct Intrusive_Link* constructor_derived_Intrusive_Link
	(struct Intrusive_Link* base, const size_t sizeof_base, const size_t sizeof_derived)
{
	struct Intrusive_Link* derived = calloc(1,sizeof_derived); // returned
	memcpy(derived,base,sizeof_base); // shallow copy of the base.

	assert(base->derived == NULL);
	base->derived = derived;

	return derived;
}

static const struct const_Intrusive_Link* constructor_derived_const_Intrusive_Link
	(const struct const_Intrusive_Link* base, const size_t sizeof_base, const size_t sizeof_derived)
{
	return (const struct const_Intrusive_Link*)
		constructor_derived_Intrusive_Link((struct Intrusive_Link*)base,sizeof_base,sizeof_derived);
}

static struct Intrusive_Link* constructor_base_Intrusive_Link
	(struct Intrusive_Link* derived, const size_t sizeof_base)
{
	struct Intrusive_Link* base = calloc(1,sizeof_base); // returned
	memcpy(base,derived,sizeof_base); // shallow copy of the derived.

	base->derived = NULL;

	// Set the derived link's `derived` member to the base to fix computational element pointers. The derived list
	// is immediately destructed after fixing the pointers so this does not impact any other code.
	derived->derived = base;

	return base;
}

static const struct const_Intrusive_Link* constructor_base_const_Intrusive_Link
	(const struct const_Intrusive_Link* derived, const size_t sizeof_base)
{
	return (const struct const_Intrusive_Link*)
		constructor_base_Intrusive_Link((struct Intrusive_Link*)derived,sizeof_base);
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
