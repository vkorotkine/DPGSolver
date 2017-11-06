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

#include "face.h"
#include "face_solver.h"
#include "face_solver_dg.h"
#include "face_solver_dpg.h"
#include "face_solver_dg_complex.h"
#include "volume.h"
#include "volume_solver.h"
#include "volume_solver_dg.h"
#include "volume_solver_dpg.h"
#include "volume_solver_dg_complex.h"

#include "element.h"
#include "geometry_element.h"
#include "element_plotting.h"
#include "element_solution.h"
#include "element_error.h"
#include "element_solver_dg.h"
#include "element_solver_dpg.h"

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

/// Container holding information used for \ref constructor_derived_computational_elements.
struct Derived_Comp_Elements_Info {
	int list_name[2]; ///< Names of the volume and face computational element lists.

	size_t sizeof_base[2];    ///< Sizes of the base volume and face computational elements.
	size_t sizeof_derived[2]; ///< Sizes of the derived volume and face computational elements.

	constructor_derived_Volume_fptr constructor_derived_Volume; ///< \ref constructor_derived_Volume_fptr.
	constructor_derived_Face_fptr   constructor_derived_Face;   ///< \ref constructor_derived_Face_fptr.

	destructor_derived_Volume_fptr destructor_derived_Volume; ///< \ref destructor_derived_Volume_fptr.
	destructor_derived_Face_fptr   destructor_derived_Face;   ///< \ref destructor_derived_Face_fptr.
};

/// Container holding information used for \ref constructor_derived_Elements and \ref destructor_derived_Elements.
struct Derived_Elements_Info {
	size_t sizeof_base,    ///< The size of the base element.
	       sizeof_derived; ///< The size of the derived element.

	constructor_derived_Element_fptr constructor_derived_Element; ///< \ref constructor_derived_Element_fptr.
	destructor_derived_Element_fptr  destructor_derived_Element;  ///< \ref destructor_derived_Element_fptr.
};

/** \brief Return a stack-allocated \ref Derived_Comp_Elements_Info container for
 *         \ref constructor_derived_computational_elements.
 *  \return See brief. */
static struct Derived_Comp_Elements_Info get_c_Derived_Comp_Elements_Info
	(const int derived_category,  ///< Defined for \ref constructor_derived_computational_elements.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Return a stack-allocated \ref Derived_Comp_Elements_Info container
 *         for \ref destructor_derived_computational_elements.
 *  \return See brief. */
static struct Derived_Comp_Elements_Info get_d_Derived_Comp_Elements_Info
	(const int base_category,   ///< Defined for \ref destructor_derived_computational_elements.
	 const int derived_category ///< The current category of the volume/face lists.
	);

/** \brief Return a stack-allocated \ref Derived_Elements_Info container for \ref constructor_derived_Elements.
 *  \return See brief. */
static struct Derived_Elements_Info get_c_Derived_Elements_Info
	(const int derived_name ///< Defined for \ref constructor_derived_Elements.
	);

/** \brief Return a stack-allocated \ref Derived_Elements_Info container for \ref destructor_derived_Elements.
 *  \return See brief. */
static struct Derived_Elements_Info get_d_Derived_Elements_Info
	(const int base_name,   ///< Defined for \ref destructor_derived_Elements.
	 const int derived_name ///< Defined for \ref constructor_derived_Elements.
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

/// \brief Version of \ref update_computational_element_list_pointers acting only on volume pointers.
static void update_volume_list_pointers
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Version of \ref update_computational_element_list_pointers acting only on face pointers.
static void update_face_list_pointers
	(const struct Simulation* sim ///< \ref Simulation.
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

/** \brief See return.
 *  \return The index of the intrusive list category corresponding to the current volume and face list names. */
static int get_list_category
	(const struct Simulation* sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void constructor_derived_computational_elements (struct Simulation* sim, const int derived_category)
{
	struct Derived_Comp_Elements_Info de_i = get_c_Derived_Comp_Elements_Info(derived_category,sim);

	struct Intrusive_List* base[2]    = { sim->volumes, sim->faces, };
	struct Intrusive_List* derived[2] = { NULL, NULL, };

	// Reserve memory for the derived lists.
	for (int i = 0; i < 2; ++i) {
		derived[i] = constructor_empty_IL(de_i.list_name[i],base[i]);
		for (struct Intrusive_Link* curr = base[i]->first; curr; curr = curr->next) {
			push_back_IL(derived[i],
			             constructor_derived_Intrusive_Link(curr,de_i.sizeof_base[i],de_i.sizeof_derived[i]));
		}
	}

	sim->volumes = derived[0];
	sim->faces   = derived[1];

	// Perform construction specific to the derived lists and update pointers.
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; ) {
		struct Intrusive_Link* next = curr->next;
		de_i.constructor_derived_Volume((struct Volume*)curr,sim);
		curr = next;
	}
	update_face_list_pointers(sim);

	for (struct Intrusive_Link* curr = sim->faces->first; curr; ) {
		struct Intrusive_Link* next = curr->next;
		de_i.constructor_derived_Face((struct Face*)curr,sim);
		curr = next;
	}
	update_volume_list_pointers(sim);

	// Destruct the base lists.
	destructor_IL_base(sim->volumes);
	destructor_IL_base(sim->faces);
}

void destructor_derived_computational_elements (struct Simulation* sim, const int base_category)
{
	const int derived_category = get_list_category(sim);
	struct Derived_Comp_Elements_Info de_i = get_d_Derived_Comp_Elements_Info(base_category,derived_category);

	// Perform destruction specific to the derived lists.
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; ) {
		struct Intrusive_Link* next = curr->next;
		de_i.destructor_derived_Volume((struct Volume*)curr);
		curr = next;
	}

	for (struct Intrusive_Link* curr = sim->faces->first; curr; ) {
		struct Intrusive_Link* next = curr->next;
		de_i.destructor_derived_Face((struct Face*)curr);
		curr = next;
	}

	// Reserve memory for the base lists.
	int ind_list = -1;

	++ind_list;
	struct Intrusive_List* volumes_prev = sim->volumes;
	sim->volumes = constructor_empty_IL(de_i.list_name[ind_list],NULL); // keep
	for (struct Intrusive_Link* curr = volumes_prev->first; curr; curr = curr->next)
		push_back_IL(sim->volumes,constructor_base_Intrusive_Link(curr,de_i.sizeof_base[ind_list]));

	++ind_list;
	struct Intrusive_List* faces_prev = sim->faces;
	sim->faces = constructor_empty_IL(de_i.list_name[ind_list],NULL); // keep
	for (struct Intrusive_Link* curr = faces_prev->first; curr; curr = curr->next)
		push_back_IL(sim->faces,constructor_base_Intrusive_Link(curr,de_i.sizeof_base[ind_list]));

	// Update pointers to the base computational elements to point to the derived computational elements.
	update_computational_element_list_pointers(sim);

	// Destruct the derived lists.
	destructor_IL(volumes_prev);
	destructor_IL(faces_prev);
}

void constructor_derived_Elements (struct Simulation* sim, const int derived_name)
{
	struct Derived_Elements_Info de_i = get_c_Derived_Elements_Info(derived_name);

	const struct const_Intrusive_List* base = sim->elements;

	// Reserve memory for the derived element list.
	const struct const_Intrusive_List* elements = constructor_empty_const_IL(derived_name,base); // moved
	for (const struct const_Intrusive_Link* curr = base->first; curr; curr = curr->next)
		push_back_const_IL(elements,
		                   constructor_derived_const_Intrusive_Link(curr,de_i.sizeof_base,de_i.sizeof_derived));
	set_tp_sub_elements((struct Intrusive_List*)elements);
	set_face_elements((struct Intrusive_List*)elements);

	sim->elements = elements;

	// Perform construction specific to the derived element list.
	for (const struct const_Intrusive_Link* curr = sim->elements->first; curr; ) {
		const struct const_Intrusive_Link* next = curr->next;
		if (((struct Element*)curr)->present)
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
		if (((struct Element*)curr)->present)
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
	destructor_const_IL(elements_prev);
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

static struct Derived_Comp_Elements_Info get_c_Derived_Comp_Elements_Info
	(const int derived_category, const struct Simulation* sim)
{
	struct Derived_Comp_Elements_Info de_info;

	switch (derived_category) {
	case IL_SOLVER:
		assert(sim->volumes->name == IL_VOLUME);
		assert(sim->faces->name   == IL_FACE);
		de_info.list_name[0] = IL_SOLVER_VOLUME;
		de_info.list_name[1] = IL_SOLVER_FACE;
		de_info.sizeof_base[0] = sizeof(struct Volume);
		de_info.sizeof_base[1] = sizeof(struct Face);
		de_info.sizeof_derived[0] = sizeof(struct Solver_Volume);
		de_info.sizeof_derived[1] = sizeof(struct Solver_Face);
		de_info.constructor_derived_Volume = constructor_derived_Solver_Volume;
		de_info.constructor_derived_Face   = constructor_derived_Solver_Face;
		break;
	case IL_SOLVER_DG:
		assert(sim->volumes->name == IL_SOLVER_VOLUME);
		assert(sim->faces->name   == IL_SOLVER_FACE);
		de_info.list_name[0] = IL_VOLUME_SOLVER_DG;
		de_info.list_name[1] = IL_FACE_SOLVER_DG;
		de_info.sizeof_base[0] = sizeof(struct Solver_Volume);
		de_info.sizeof_base[1] = sizeof(struct Solver_Face);
		de_info.sizeof_derived[0] = sizeof(struct DG_Solver_Volume);
		de_info.sizeof_derived[1] = sizeof(struct DG_Solver_Face);
		de_info.constructor_derived_Volume = constructor_derived_DG_Solver_Volume;
		de_info.constructor_derived_Face   = constructor_derived_DG_Solver_Face;
		break;
	case IL_SOLVER_DPG:
		assert(sim->volumes->name == IL_SOLVER_VOLUME);
		assert(sim->faces->name   == IL_SOLVER_FACE);
		de_info.list_name[0] = IL_VOLUME_SOLVER_DPG;
		de_info.list_name[1] = IL_FACE_SOLVER_DPG;
		de_info.sizeof_base[0] = sizeof(struct Solver_Volume);
		de_info.sizeof_base[1] = sizeof(struct Solver_Face);
		de_info.sizeof_derived[0] = sizeof(struct DPG_Solver_Volume);
		de_info.sizeof_derived[1] = sizeof(struct DPG_Solver_Face);
		de_info.constructor_derived_Volume = constructor_derived_DPG_Solver_Volume;
		de_info.constructor_derived_Face   = constructor_derived_DPG_Solver_Face;
		break;
	case IL_SOLVER_DG_COMPLEX:
		assert(sim->volumes->name == IL_VOLUME_SOLVER_DG);
		assert(sim->faces->name   == IL_FACE_SOLVER_DG);
		de_info.list_name[0] = IL_VOLUME_SOLVER_DG_COMPLEX;
		de_info.list_name[1] = IL_FACE_SOLVER_DG_COMPLEX;
		de_info.sizeof_base[0] = sizeof(struct DG_Solver_Volume);
		de_info.sizeof_base[1] = sizeof(struct DG_Solver_Face);
		de_info.sizeof_derived[0] = sizeof(struct Complex_DG_Solver_Volume);
		de_info.sizeof_derived[1] = sizeof(struct Complex_DG_Solver_Face);
		de_info.constructor_derived_Volume = constructor_derived_Complex_DG_Solver_Volume;
		de_info.constructor_derived_Face   = constructor_derived_Complex_DG_Solver_Face;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",derived_category);
		break;
	}
	return de_info;
}

static struct Derived_Comp_Elements_Info get_d_Derived_Comp_Elements_Info
	(const int base_category, const int derived_category)
{
	struct Derived_Comp_Elements_Info de_info;

	switch (base_category) {
	case IL_BASE:
		de_info.list_name[0] = IL_VOLUME;
		de_info.list_name[1] = IL_FACE;
		de_info.sizeof_base[0] = sizeof(struct Volume);
		de_info.sizeof_base[1] = sizeof(struct Face);
		break;
	case IL_SOLVER:
		de_info.list_name[0] = IL_SOLVER_VOLUME;
		de_info.list_name[1] = IL_SOLVER_FACE;
		de_info.sizeof_base[0] = sizeof(struct Solver_Volume);
		de_info.sizeof_base[1] = sizeof(struct Solver_Face);
		break;
	case IL_SOLVER_DG:
		de_info.list_name[0] = IL_VOLUME_SOLVER_DG;
		de_info.list_name[1] = IL_FACE_SOLVER_DG;
		de_info.sizeof_base[0] = sizeof(struct DG_Solver_Volume);
		de_info.sizeof_base[1] = sizeof(struct DG_Solver_Face);
		break;
	case IL_SOLVER_DPG:
		de_info.list_name[0] = IL_VOLUME_SOLVER_DPG;
		de_info.list_name[1] = IL_FACE_SOLVER_DPG;
		de_info.sizeof_base[0] = sizeof(struct DPG_Solver_Volume);
		de_info.sizeof_base[1] = sizeof(struct DPG_Solver_Face);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",base_category);
		break;
	}

	switch (derived_category) {
	case IL_SOLVER:
		assert(base_category == IL_BASE);
		de_info.destructor_derived_Volume = destructor_derived_Solver_Volume;
		de_info.destructor_derived_Face   = destructor_derived_Solver_Face;
		break;
	case IL_SOLVER_DG:
		assert(base_category == IL_SOLVER);
		de_info.destructor_derived_Volume = destructor_derived_DG_Solver_Volume;
		de_info.destructor_derived_Face   = destructor_derived_DG_Solver_Face;
		break;
	case IL_SOLVER_DG_COMPLEX:
		assert(base_category == IL_SOLVER_DG);
		de_info.destructor_derived_Volume = destructor_derived_Complex_DG_Solver_Volume;
		de_info.destructor_derived_Face   = destructor_derived_Complex_DG_Solver_Face;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",derived_category);
		break;
	}
	return de_info;
}

static struct Derived_Elements_Info get_c_Derived_Elements_Info (const int derived_name)
{
	struct Derived_Elements_Info de_info;

	switch (derived_name) {
	case IL_GEOMETRY_ELEMENT:
		assert(sizeof(struct Geometry_Element) == sizeof(struct const_Geometry_Element));
		de_info.sizeof_base    = sizeof(struct Element);
		de_info.sizeof_derived = sizeof(struct Geometry_Element);
		de_info.constructor_derived_Element = constructor_derived_Geometry_Element;
		break;
	case IL_PLOTTING_ELEMENT:
		assert(sizeof(struct Plotting_Element) == sizeof(struct const_Plotting_Element));
		de_info.sizeof_base    = sizeof(struct Element);
		de_info.sizeof_derived = sizeof(struct Plotting_Element);
		de_info.constructor_derived_Element = constructor_derived_Plotting_Element;
		break;
	case IL_SOLUTION_ELEMENT:
		assert(sizeof(struct Solution_Element) == sizeof(struct const_Solution_Element));
		de_info.sizeof_base    = sizeof(struct Element);
		de_info.sizeof_derived = sizeof(struct Solution_Element);
		de_info.constructor_derived_Element = constructor_derived_Solution_Element;
		break;
	case IL_ELEMENT_ERROR:
		assert(sizeof(struct Error_Element) == sizeof(struct const_Error_Element));
		de_info.sizeof_base    = sizeof(struct Solution_Element);
		de_info.sizeof_derived = sizeof(struct Error_Element);
		de_info.constructor_derived_Element = constructor_derived_Error_Element;
		break;
	case IL_ELEMENT_SOLVER_DG:
		assert(sizeof(struct DG_Solver_Element) == sizeof(struct const_DG_Solver_Element));
		de_info.sizeof_base    = sizeof(struct Element);
		de_info.sizeof_derived = sizeof(struct DG_Solver_Element);
		de_info.constructor_derived_Element = constructor_derived_DG_Solver_Element;
		break;
	case IL_ELEMENT_SOLVER_DPG:
		assert(sizeof(struct DPG_Solver_Element) == sizeof(struct const_DPG_Solver_Element));
		de_info.sizeof_base    = sizeof(struct Element);
		de_info.sizeof_derived = sizeof(struct DPG_Solver_Element);
		de_info.constructor_derived_Element = constructor_derived_DPG_Solver_Element;
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
	case IL_SOLUTION_ELEMENT:
		de_info.sizeof_base = sizeof(struct Solution_Element);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",base_name);
		break;
	}

	switch (derived_name) {
	case IL_GEOMETRY_ELEMENT:
		assert(base_name == IL_ELEMENT);
		de_info.destructor_derived_Element = destructor_derived_Geometry_Element;
		break;
	case IL_PLOTTING_ELEMENT:
		assert(base_name == IL_ELEMENT);
		de_info.destructor_derived_Element = destructor_derived_Plotting_Element;
		break;
	case IL_SOLUTION_ELEMENT:
		assert(base_name == IL_ELEMENT);
		de_info.destructor_derived_Element = destructor_derived_Solution_Element;
		break;
	case IL_ELEMENT_ERROR:
		assert(base_name == IL_SOLUTION_ELEMENT);
		de_info.destructor_derived_Element = destructor_derived_Error_Element;
		break;
	case IL_ELEMENT_SOLVER_DG:
		assert(base_name == IL_ELEMENT);
		de_info.destructor_derived_Element = destructor_derived_DG_Solver_Element;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",derived_name);
		break;
	}
	return de_info;
}

static void update_computational_element_list_pointers (const struct Simulation* sim)
{
	update_volume_list_pointers(sim);
	update_face_list_pointers(sim);
}

static void update_volume_list_pointers (const struct Simulation* sim)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next)
		update_volume_pointers(curr);
}

static void update_face_list_pointers (const struct Simulation* sim)
{
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

static int get_list_category (const struct Simulation* sim)
{
	const int v_name = sim->volumes->name,
	          f_name = sim->faces->name;

	int ce_name = -1;
	switch (v_name) {
	case IL_SOLVER_VOLUME:
		assert(f_name == IL_SOLVER_FACE);
		ce_name = IL_SOLVER;
		break;
	case IL_VOLUME_SOLVER_DG:
		assert(f_name == IL_FACE_SOLVER_DG);
		ce_name = IL_SOLVER_DG;
		break;
	case IL_VOLUME_SOLVER_DG_COMPLEX:
		assert(f_name == IL_FACE_SOLVER_DG_COMPLEX);
		ce_name = IL_SOLVER_DG_COMPLEX;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",v_name);
		break;
	}
	return ce_name;
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
