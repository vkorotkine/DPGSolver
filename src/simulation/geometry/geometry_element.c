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

#include <assert.h>
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
	(struct Simulation* sim,                              ///< \ref Simulation.
	 const struct const_Intrusive_List* geometry_elements ///< \ref Geometry_Element\*s.
	);

// Interface functions ********************************************************************************************** //

void constructor_Geometry_Elements (struct Simulation*const sim)
{
	const struct const_Intrusive_List* base = sim->elements;
	assert(base->name == IL_ELEMENT);

	sim->elements = constructor_empty_const_IL(IL_GEOMETRY_ELEMENT,base); // moved
	for (const struct const_Intrusive_Link* curr = base->first; curr; curr = curr->next) {
		push_back_const_IL(sim->elements,
			(const struct const_Intrusive_Link*)constructor_Geometry_Element((struct const_Element*)curr));
	}
	set_tp_sub_elements((struct Intrusive_List*)sim->elements);

	set_up_geometry_ops(sim,sim->elements);
}

void destructor_Geometry_Elements (const struct const_Intrusive_List* geometry_elements)
{
	for (const struct const_Intrusive_Link* curr = geometry_elements->first; curr; ) {
		const struct const_Intrusive_Link* next = curr->next;
		destructor_Geometry_Element((struct Geometry_Element*) curr);
		curr = next;
	}
EXIT_ADD_SUPPORT; // Likely make a general element destructor function which adjusts the pointer in sim.
	destructor_const_IL(geometry_elements);
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

	set_derived_link(element,geometry_element);
	set_derived_link(geometry_element,NULL);

	return geometry_element;
}

static void destructor_Geometry_Element (struct Geometry_Element* element)
{
	destructor_Multiarray_Operator(element->vc0_vgc_vgc);

	destructor_Multiarray_Operator(element->cv1_vgs_vcs);
	destructor_Multiarray_Operator(element->cv1_vgc_vcc);
	destructor_Multiarray_Operator(element->cv1_vgs_vms);
	destructor_Multiarray_Operator(element->cv1_vgc_vmc);
	destructor_Multiarray_Operator(element->vv0_vms_vcs);
	destructor_Multiarray_Operator(element->vv0_vmc_vcc);

	destructor_Multiarray_Operator(element->cv0_vgs_fcs);
	destructor_Multiarray_Operator(element->cv0_vgs_fcc);
	destructor_Multiarray_Operator(element->cv0_vgc_fcs);
	destructor_Multiarray_Operator(element->cv0_vgc_fcc);

	destructor_Multiarray_Operator(element->cv0_vgs_vcs);
	destructor_Multiarray_Operator(element->cv0_vgc_vcc);
	destructor_Multiarray_Operator(element->cv0_vgs_vms);
	destructor_Multiarray_Operator(element->cv0_vgc_vmc);

	EXIT_ADD_SUPPORT; // Add destructors for geometry operators.
	destructor_Element((struct Element*) element);
}

void set_up_geometry_ops (struct Simulation* sim, const struct const_Intrusive_List* geometry_elements)
{
	for (const struct const_Intrusive_Link* curr = geometry_elements->first; curr; curr = curr->next) {
		if (((struct const_Element*)curr)->present)
			set_up_operators_element((struct Geometry_Element*)curr,sim);
	}
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
	struct const_Element* b_e = (struct const_Element*)element;

	element->vc0_vgc_vgc = constructor_operators("vc0","vgc","vgc","H_1_P_1P",sim->p_s_v,b_e,sim); // keep

	element->cv1_vgs_vcs = constructor_operators("cv1","vgs","vcs","H_1_P_1P", sim->p_s_v,b_e,sim); // keep
	element->cv1_vgc_vcc = constructor_operators("cv1","vgc","vcc","H_1_P_PM0",sim->p_s_v,b_e,sim); // keep
	element->cv1_vgs_vms = constructor_operators("cv1","vgs","vms","H_1_P_1",  sim->p_s_v,b_e,sim); // keep
	element->cv1_vgc_vmc = constructor_operators("cv1","vgc","vmc","H_1_P_PM0",sim->p_s_v,b_e,sim); // keep
	element->vv0_vms_vcs = constructor_operators("vv0","vms","vcs","H_1_P_1P", sim->p_s_v,b_e,sim); // keep
	element->vv0_vmc_vcc = constructor_operators("vv0","vmc","vcc","H_1_P_PM0",sim->p_s_v,b_e,sim); // keep

	element->cv0_vgs_fcs = constructor_operators("cv0","vgs","fcs","H_1_P_1P", sim->p_s_f,b_e,sim); // keep
print_Multiarray_Operator(element->cv0_vgs_fcs);
EXIT_UNSUPPORTED;
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
	const struct const_Element* b_e     = (const struct const_Element*)element;
	const struct const_Element* bs_e[2] = { b_e->sub_element[0], b_e->sub_element[1], };
	struct Geometry_Element* s_e[2]     = { (struct Geometry_Element*) bs_e[0],
	                                        (struct Geometry_Element*) bs_e[1], };

	struct Operators_TP ops_tp;

	// tensor-product sub-operators
	if (s_e[0]->cv0_vgs_vcs == NULL) {
		s_e[0]->cv0_vgs_vcs = constructor_operators("cv0","vgs","vcs","H_1_P_1P", sim->p_s_v,bs_e[0],sim); // destructed
		s_e[0]->cv0_vgc_vcc = constructor_operators("cv0","vgc","vcc","H_1_P_PM0",sim->p_s_v,bs_e[0],sim); // destructed
		s_e[0]->cv0_vgs_vms = constructor_operators("cv0","vgs","vms","H_1_P_1",  sim->p_s_v,bs_e[0],sim); // destructed
		s_e[0]->cv0_vgc_vmc = constructor_operators("cv0","vgc","vmc","H_1_P_PM0",sim->p_s_v,bs_e[0],sim); // destructed
	}
	if (s_e[1]->cv0_vgs_vcs == NULL) {
		s_e[1]->cv0_vgs_vcs = constructor_operators("cv0","vgs","vcs","H_1_P_1P", sim->p_s_v,bs_e[1],sim); // destructed
		s_e[1]->cv0_vgc_vcc = constructor_operators("cv0","vgc","vcc","H_1_P_PM0",sim->p_s_v,bs_e[1],sim); // destructed
		s_e[1]->cv0_vgs_vms = constructor_operators("cv0","vgs","vms","H_1_P_1",  sim->p_s_v,bs_e[1],sim); // destructed
		s_e[1]->cv0_vgc_vmc = constructor_operators("cv0","vgc","vmc","H_1_P_PM0",sim->p_s_v,bs_e[1],sim); // destructed
	}

	set_operators_tp(&ops_tp,s_e[0]->cv0_vgs_vcs,s_e[0]->cv1_vgs_vcs,s_e[1]->cv0_vgs_vcs,s_e[1]->cv1_vgs_vcs);
	element->cv1_vgs_vcs = constructor_operators_tp("cv1","vgs","vcs","H_1_P_1P",sim->p_s_v,b_e,sim,&ops_tp); // keep

	set_operators_tp(&ops_tp,s_e[0]->cv0_vgc_vcc,s_e[0]->cv1_vgc_vcc,s_e[1]->cv0_vgc_vcc,s_e[1]->cv1_vgc_vcc);
	element->cv1_vgc_vcc = constructor_operators_tp("cv1","vgs","vcc","H_1_P_PM0",sim->p_s_v,b_e,sim,&ops_tp); // keep

	set_operators_tp(&ops_tp,s_e[0]->cv0_vgs_vms,s_e[0]->cv1_vgs_vms,s_e[1]->cv0_vgs_vms,s_e[1]->cv1_vgs_vms);
	element->cv1_vgs_vms = constructor_operators_tp("cv1","vgs","vms","H_1_P_1",sim->p_s_v,b_e,sim,&ops_tp); // keep

	set_operators_tp(&ops_tp,s_e[0]->cv0_vgc_vmc,s_e[0]->cv1_vgc_vmc,s_e[1]->cv0_vgc_vmc,s_e[1]->cv1_vgc_vmc);
	element->cv1_vgc_vmc = constructor_operators_tp("cv1","vgs","vmc","H_1_P_PM0",sim->p_s_v,b_e,sim,&ops_tp); // keep

	set_operators_tp(&ops_tp,s_e[0]->vv0_vms_vcs,NULL,s_e[1]->vv0_vms_vcs,NULL);
	element->vv0_vms_vcs = constructor_operators_tp("vv0","vms","vcs","H_1_P_1P", sim->p_s_v,b_e,sim,&ops_tp); // keep

	set_operators_tp(&ops_tp,s_e[0]->vv0_vmc_vcc,NULL,s_e[1]->vv0_vmc_vcc,NULL);
	element->vv0_vmc_vcc = constructor_operators_tp("vv0","vmc","vcc","H_1_P_PM0",sim->p_s_v,b_e,sim,&ops_tp); // keep
}
