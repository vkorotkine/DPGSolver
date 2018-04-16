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

#include "element_plotting.h"

#include "macros.h"
#include "definitions_elements.h"

#include "element_operators.h"
#include "element_operators_tp.h"
#include "multiarray_operator.h"
#include "nodes_plotting.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

/// \brief Constructor for a derived \ref Plotting_Element using the standard operators.
static void constructor_derived_Plotting_Element_std
	(struct Element* element,     ///< Defined for \ref constructor_derived_Plotting_Element.
	 const struct Simulation* sim ///< Defined for \ref constructor_derived_Plotting_Element.
	);

/// \brief Constructor for a derived \ref Plotting_Element using the tensor-product of sub-element operators.
static void constructor_derived_Plotting_Element_tp
	(struct Element* element,     ///< Defined for \ref constructor_derived_Plotting_Element.
	 const struct Simulation* sim ///< Defined for \ref constructor_derived_Plotting_Element.
	);

// Interface functions ********************************************************************************************** //

void constructor_derived_Plotting_Element (struct Element* element_ptr, const struct Simulation* sim)
{
	if (element_ptr->type == POINT)
		return;

	switch (element_ptr->type) {
	case LINE: case TRI: case TET: case PYR:
		constructor_derived_Plotting_Element_std(element_ptr,sim);
		break;
	case QUAD: case HEX: case WEDGE:
		constructor_derived_Plotting_Element_tp(element_ptr,sim);
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}

	const int* p_ref = sim->p_ref;
	const int n_p = p_ref[1]+1;

	struct Plotting_Element* element = (struct Plotting_Element*) element_ptr;
	element->n_p = n_p;
	element->p_nodes = calloc((size_t)n_p , sizeof *element->p_nodes); // free
	for (int i = p_ref[0]; i <= p_ref[1]; ++i)
		element->p_nodes[i] = constructor_const_Plotting_Nodes(i,element_ptr->type); // destructed
}

void destructor_derived_Plotting_Element (struct Element* element_ptr)
{
	if (element_ptr->type == POINT)
		return;

	struct Plotting_Element* element = (struct Plotting_Element*) element_ptr;
	destructor_Multiarray_Operator(element->cv0_vgs_vp);
	destructor_Multiarray_Operator(element->cv0_vgc_vp);
	destructor_Multiarray_Operator(element->cv0_vs_vp);
	destructor_Multiarray_Operator(element->cv0_vr_vp);
	destructor_Multiarray_Operator(element->cv0_vt_vp);

	const int n_p = element->n_p;
	for (int i = 0; i < n_p; ++i) {
		if (element->p_nodes[i])
			destructor_const_Plotting_Nodes(element->p_nodes[i]);
	}
	free(element->p_nodes);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void constructor_derived_Plotting_Element_std (struct Element* element_ptr, const struct Simulation* sim)
{
	struct Plotting_Element* element = (struct Plotting_Element*) element_ptr;

	struct const_Element* b_e = (struct const_Element*)element;

	element->cv0_vgs_vp = constructor_operators("cv0","vgs","vpA","H_1_P_1P", b_e,sim); // destructed
	element->cv0_vgc_vp = constructor_operators("cv0","vgc","vpA","H_1_P_PM0",b_e,sim); // destructed
	element->cv0_vs_vp  = constructor_operators("cv0","vsA","vpA","H_1_P_PM0",b_e,sim); // destructed
	element->cv0_vr_vp  = constructor_operators("cv0","vrA","vpA","H_1_P_PM0",b_e,sim); // destructed
	element->cv0_vt_vp  = constructor_operators("cv0","vtA","vpA","H_1_P_PM0",b_e,sim); // destructed
}

static void constructor_derived_Plotting_Element_tp (struct Element* element_ptr, const struct Simulation* sim)
{
	struct Plotting_Element* element = (struct Plotting_Element*) element_ptr;

	const struct const_Element* b_e     = (const struct const_Element*)element;
	const struct const_Element* bs_e[2] = { b_e->sub_element[0], b_e->sub_element[1], };
	struct Plotting_Element* s_e[2]     = { (struct Plotting_Element*) bs_e[0],
	                                        (struct Plotting_Element*) bs_e[1], };

	struct Operators_TP ops_tp;

	set_operators_tp(&ops_tp,s_e[0]->cv0_vgs_vp,NULL,s_e[1]->cv0_vgs_vp,NULL);
	element->cv0_vgs_vp = constructor_operators_tp("cv0","vgs","vpA","H_1_P_1P",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->cv0_vgc_vp,NULL,s_e[1]->cv0_vgc_vp,NULL);
	element->cv0_vgc_vp = constructor_operators_tp("cv0","vgc","vpA","H_1_P_PM0",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->cv0_vs_vp,NULL,s_e[1]->cv0_vs_vp,NULL);
	element->cv0_vs_vp = constructor_operators_tp("cv0","vsA","vpA","H_1_P_PM0",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->cv0_vr_vp,NULL,s_e[1]->cv0_vr_vp,NULL);
	element->cv0_vr_vp = constructor_operators_tp("cv0","vrA","vpA","H_1_P_PM0",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->cv0_vt_vp,NULL,s_e[1]->cv0_vt_vp,NULL);
	element->cv0_vt_vp = constructor_operators_tp("cv0","vtA","vpA","H_1_P_PM0",b_e,sim,&ops_tp); // destructed
}
