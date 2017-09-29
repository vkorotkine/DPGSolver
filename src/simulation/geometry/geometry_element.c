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
/// \file

#include "geometry_element.h"

#include <string.h>

#include "macros.h"
#include "definitions_intrusive.h"
#include "definitions_cubature.h"
#include "definitions_element_operators.h"
#include "definitions_elements.h"

#include "multiarray.h"
#include "multiarray_operator.h"

#include "simulation.h"
#include "element_operators.h"
#include "element_operators_tp.h"
#include "const_cast.h"

// Static function declarations ************************************************************************************* //

/** \brief Constructor for an individual \ref Element.
 *  \return Standard. */
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
	struct Intrusive_List* geometry_elements = constructor_empty_IL(IL_GEOMETRY_ELEMENT); // returned

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
	const_cast_void1(&element->derived,geometry_element);

	return geometry_element;
}

static void destructor_Geometry_Element (struct Geometry_Element* element)
{
	destructor_Element((struct Element*) element);
}

void set_up_geometry_ops (struct Simulation* sim, struct Intrusive_List* geometry_elements)
{
	for (struct Intrusive_Link* curr = geometry_elements->first; curr; curr = curr->next)
		set_up_operators_element((struct Geometry_Element*)curr,sim);
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

/** \brief Set up the wedge operators for the \ref Geometry_Element\*s using a tensor-product of lower-dimensional
           operators. */
static void set_up_operators_tp_wedge
	(struct Geometry_Element* element, ///< Defined for \ref set_up_operators_standard.
	 const struct Simulation* sim      ///< Defined for \ref set_up_operators_standard.
	);

static void set_up_operators_standard (struct Geometry_Element* element, const struct Simulation* sim)
{
	struct const_Element* base_element = (struct const_Element*)element;

	element->cv1_vgs_vcs = constructor_operators("cv1","vgs","vcs","H_1_P_1",  sim->p_s_v,base_element,sim); // keep
	element->cv1_vgc_vcc = constructor_operators("cv1","vgc","vcc","H_1_P_PM0",sim->p_s_v,base_element,sim); // keep
}

static void set_up_operators_tensor_product (struct Geometry_Element* element, const struct Simulation* sim)
{
	struct const_Element* base_element = (struct const_Element*)element;

	// \note The wedge operator set up should also work for the n-cube but would be less efficient as it will build
	//       redundant operators.
	switch (base_element->s_type) {
	case ST_TP:
//		set_up_operators_tp_n_cube(element,sim); /// \todo Use the n_cube operator set up when all is working.
		set_up_operators_tp_wedge(element,sim);
		break;
	case ST_WEDGE:
		set_up_operators_tp_wedge(element,sim);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",base_element->s_type);
		break;
	}
}

// Level 3 ********************************************************************************************************** //

static void set_up_operators_tp_wedge (struct Geometry_Element* element, const struct Simulation* sim)
{
	struct const_Element* b_e = (struct const_Element*)element;

	struct Geometry_Element* s_e[2] = { (struct Geometry_Element*) b_e->sub_element[0]->derived,
	                                    (struct Geometry_Element*) b_e->sub_element[1]->derived, };

	struct Operators_TP ops_tp;

	const struct Multiarray_Operator* cv0_vgs_vcs[2] = // destructed
		{ constructor_operators("cv0","vgs","vcs","H_1_P_1",sim->p_s_v,b_e->sub_element[0],sim),
		  constructor_operators("cv0","vgs","vcs","H_1_P_1",sim->p_s_v,b_e->sub_element[1],sim), };
	set_operators_tp(&ops_tp,cv0_vgs_vcs[0],s_e[0]->cv1_vgs_vcs,cv0_vgs_vcs[1],s_e[1]->cv1_vgs_vcs);
	element->cv1_vgs_vcs = constructor_operators_tp("cv1","vgs","vcs","H_1_P_1",b_e,sim,&ops_tp); // keep

	destructor_Multiarray2_Operator(cv0_vgs_vcs);



//	element->cv1_vgc_vcc =
	EXIT_ADD_SUPPORT;
}
