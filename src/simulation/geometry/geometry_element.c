// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
///	\file

#include "geometry_element.h"

#include <string.h>

#include "macros.h"
#include "definitions_intrusive.h"

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

///	\brief Set up the geometry related operators.
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
	struct Intrusive_List* geometry_elements    = constructor_empty_IL(IL_GEOMETRY_ELEMENT);

	for (const struct Intrusive_Link* curr = elements->first; curr; curr = curr->next)
		push_back_IL(geometry_elements,
		             (struct Intrusive_Link*) constructor_Geometry_Element((struct const_Element*) curr));

	set_up_geometry_ops(sim,geometry_elements);

	destructor_IL((struct Intrusive_List*)elements);
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
UNUSED(sim);
UNUSED(geometry_elements);
}
