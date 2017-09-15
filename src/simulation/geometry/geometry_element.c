// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/// \file

#include "geometry_element.h"

#include <string.h>

#include "macros.h"
#include "definitions_intrusive.h"
#include "definitions_elements.h"

#include "multiarray.h"

#include "simulation.h"

// Static function declarations ************************************************************************************* //

/// \brief Constructor for an individual \ref Element.
static struct Geometry_Element* constructor_Geometry_Element
	(struct const_Element* element ///< The base \ref Element.
	);

/// \brief Destructor for an individual \ref Element.
static void destructor_Geometry_Element
	(struct Geometry_Element* element ///< Standard.
	);

/// \brief Set up the geometry related operators.
void set_up_geometry_ops
	(struct Simulation* sim,                  ///< \ref Simulation.
	 struct Intrusive_List* geometry_elements ///< \ref Geometry_Element\*s.
	);

// Interface functions ********************************************************************************************** //

struct const_Intrusive_List* constructor_Geometry_Elements (struct Simulation*const sim)
{
	if (sim->elements->name != IL_ELEMENT)
		EXIT_UNSUPPORTED;

	const struct const_Intrusive_List* elements = sim->elements;
	struct Intrusive_List* geometry_elements    = constructor_empty_IL(IL_GEOMETRY_ELEMENT); // returned

	for (const struct Intrusive_Link* curr = elements->first; curr; curr = curr->next)
		push_back_IL(geometry_elements,
		             (struct Intrusive_Link*) constructor_Geometry_Element((struct const_Element*) curr));

	set_up_geometry_ops(sim,geometry_elements);

	return (struct const_Intrusive_List*) geometry_elements;
}

void destructor_Geometry_Elements (struct Intrusive_List* geometry_elements)
{
	for (const struct Intrusive_Link* curr = geometry_elements->first; curr; ) {
		struct Intrusive_Link* next = curr->next;
		destructor_Geometry_Element((struct Geometry_Element*) curr);
		curr = next;
	}
	destructor_IL(geometry_elements);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Set up the operators for the \ref Geometry_Element\*s.
static void set_up_operators_element
	(struct Geometry_Element* element, ///< \ref Geometry_Element.
	 const struct Simulation* sim      ///< \ref Simulation.
	);

static struct Geometry_Element* constructor_Geometry_Element (struct const_Element* element)
{
	struct Geometry_Element* geometry_element = calloc(1,sizeof *geometry_element); // returned

	memcpy(&geometry_element->element,element,sizeof *element); // shallow copy of the base.

	return geometry_element;
}

static void destructor_Geometry_Element (struct Geometry_Element* element)
{
	destructor_Element((struct Element*) element);
}

void set_up_geometry_ops (struct Simulation* sim, struct Intrusive_List* geometry_elements)
{
	for (struct Intrusive_Link* curr = geometry_elements->first; curr; curr = curr->next) {
		set_up_operators_element((struct Geometry_Element*)curr,sim);
		printf("geom_op: %d\n",((struct Element*)curr)->type);
	}
EXIT_UNSUPPORTED;
}

// Level 1 ********************************************************************************************************** //

/// \brief Set up the operators for the \ref Geometry_Element\*s using the standard method.
static void set_up_operators_standard
	(struct Geometry_Element* element, ///< \ref Geometry_Element.
	 const struct Simulation* sim      ///< \ref Simulation.
	);

/// \brief Set up the operators for the \ref Geometry_Element\*s using a tensor-product of lower-dimensional operators.
static void set_up_operators_tensor_product
	(struct Geometry_Element* element, ///< \ref Geometry_Element.
	 const struct Simulation* sim      ///< \ref Simulation.
	);

static void set_up_operators_element (struct Geometry_Element* element, const struct Simulation* sim)
{
	switch (((struct Element*)element)->type) {
	case POINT:
		return;
		break;
	case LINE: case TRI: case TET: case PYR:
		set_up_operators_standard(element,sim);
		break;
	case QUAD: case HEX: case WEDGE:
		set_up_operators_tensor_product(element,sim);
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
}

// Level 2 ********************************************************************************************************** //

static void set_up_operators_standard (struct Geometry_Element* element, const struct Simulation* sim)
{
	const int n_hp = 3;

	const int d = ((struct Element*)element)->d;

// Modify with values from simulation.
UNUSED(sim);
const int max_p = 3;
const int max_pp1 = max_p+1;
const int max_h = 3;

	const ptrdiff_t ext_v0[] = { max_pp1, max_pp1, max_h, },
	                ext_v1[] = { d, max_pp1, max_pp1, max_h, };
UNUSED(ext_v0);

	struct Multiarray_Matrix_d* ED_vg_vc = constructor_empty_Multiarray_Matrix_d(false,n_hp+1,ext_v1); // keep

//	struct const_Vector_Cubature cub_vg = constructor_const_Vector_Cubature(); // destructed

	element->ED_vg_vc = ED_vg_vc;

//	destructor_const_Vector_Cubature(cub_vg);
}

static void set_up_operators_tensor_product (struct Geometry_Element* element, const struct Simulation* sim)
{
UNUSED(element);
UNUSED(sim);
	EXIT_ADD_SUPPORT;
}
