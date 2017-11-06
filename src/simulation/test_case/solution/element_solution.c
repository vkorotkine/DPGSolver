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

#include "element_solution.h"

#include "macros.h"
#include "definitions_elements.h"

#include "element_operators.h"
#include "element_operators_tp.h"
#include "multiarray_operator.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

/// \brief Constructor for a derived \ref Solution_Element using the standard operators.
static void constructor_derived_Solution_Element_std
	(struct Element* element,     ///< Defined for \ref constructor_derived_Solution_Element.
	 const struct Simulation* sim ///< Defined for \ref constructor_derived_Solution_Element.
	);

/// \brief Constructor for a derived \ref Solution_Element using the tensor-product of sub-element operators.
static void constructor_derived_Solution_Element_tp
	(struct Element* element,     ///< Defined for \ref constructor_derived_Solution_Element.
	 const struct Simulation* sim ///< Defined for \ref constructor_derived_Solution_Element.
	);

// Interface functions ********************************************************************************************** //

void constructor_derived_Solution_Element (struct Element* element_ptr, const struct Simulation* sim)
{
	switch (element_ptr->type) {
	case LINE: case TRI: case TET: case PYR:
		constructor_derived_Solution_Element_std(element_ptr,sim);
		break;
	case QUAD: case HEX: case WEDGE:
		constructor_derived_Solution_Element_tp(element_ptr,sim);
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
}

void destructor_derived_Solution_Element (struct Element* element_ptr)
{
	struct Solution_Element* element = (struct Solution_Element*) element_ptr;

	destructor_Multiarray2_Operator(element->cv0_vg_vs);
	destructor_Multiarray2_Operator(element->cv0_vg_vc);
	destructor_Multiarray2_Operator_conditional(element->cv0_vs_vc);
	destructor_Multiarray_Operator_conditional(element->cv0_vs_vs);
	destructor_Multiarray_Operator(element->vc0_vs_vs);

	destructor_Multiarray2_Operator_conditional(element->cv0_vg_vr);
	destructor_Multiarray2_Operator_conditional(element->cv0_vg_vr);

	destructor_Multiarray2_Operator_conditional(element->cv0_vg_fs);
	destructor_Multiarray2_Operator_conditional(element->cv0_vg_fs);

	destructor_Multiarray2_Operator_conditional(element->cv0_vg_fr);
	destructor_Multiarray2_Operator_conditional(element->cv0_vg_fr);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void constructor_derived_Solution_Element_std (struct Element* element_ptr, const struct Simulation* sim)
{
	struct Solution_Element* e = (struct Solution_Element*) element_ptr;

	struct const_Element* b_e = (struct const_Element*)e;

	e->cv0_vg_vs[0] = constructor_operators("cv0","vgs","vsA","H_1_P_1P", sim->p_s_v,b_e,sim); // destructed
	e->cv0_vg_vs[1] = constructor_operators("cv0","vgc","vsA","H_1_P_PM0",sim->p_s_v,b_e,sim); // destructed
	e->cv0_vg_vc[0] = constructor_operators("cv0","vgs","vcs","H_1_P_1P", sim->p_s_v,b_e,sim); // destructed
	e->cv0_vg_vc[1] = constructor_operators("cv0","vgc","vcc","H_1_P_PM0",sim->p_s_v,b_e,sim); // destructed
	e->vc0_vs_vs    = constructor_operators("vc0","vsA","vsA","H_1_P_PM0",sim->p_s_v,b_e,sim); // destructed
	switch (sim->method) {
	case METHOD_DG:
		e->cv0_vs_vc[0] = constructor_operators("cv0","vsA","vcs","H_1_P_PM0",sim->p_s_v,b_e,sim); // destructed
		e->cv0_vs_vc[1] = constructor_operators("cv0","vsA","vcc","H_1_P_PM0",sim->p_s_v,b_e,sim); // destructed
		e->cv0_vs_vs    = constructor_operators("cv0","vsA","vsA","H_1_P_PM0",sim->p_s_v,b_e,sim); // destructed
		break;
	case METHOD_DPG:
		/// \todo Add necessary operators.
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",sim->method);
		break;
	}
}

static void constructor_derived_Solution_Element_tp (struct Element* element_ptr, const struct Simulation* sim)
{
	struct Solution_Element* e = (struct Solution_Element*) element_ptr;

	const struct const_Element* b_e     = (const struct const_Element*)e;
	const struct const_Element* bs_e[2] = { b_e->sub_element[0], b_e->sub_element[1], };
	struct Solution_Element* s_e[2]     = { (struct Solution_Element*) bs_e[0],
	                                        (struct Solution_Element*) bs_e[1], };

	struct Operators_TP ops_tp;

	set_operators_tp(&ops_tp,s_e[0]->cv0_vg_vs[0],NULL,s_e[1]->cv0_vg_vs[0],NULL);
	e->cv0_vg_vs[0] = constructor_operators_tp("cv0","vgs","vsA","H_1_P_1P", sim->p_s_v,b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->cv0_vg_vs[1],NULL,s_e[1]->cv0_vg_vs[1],NULL);
	e->cv0_vg_vs[1] = constructor_operators_tp("cv0","vgc","vsA","H_1_P_PM0",sim->p_s_v,b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->cv0_vg_vc[0],NULL,s_e[1]->cv0_vg_vc[0],NULL);
	e->cv0_vg_vc[0] = constructor_operators_tp("cv0","vgs","vsA","H_1_P_1P", sim->p_s_v,b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->cv0_vg_vc[1],NULL,s_e[1]->cv0_vg_vc[1],NULL);
	e->cv0_vg_vc[1] = constructor_operators_tp("cv0","vgc","vsA","H_1_P_PM0",sim->p_s_v,b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->vc0_vs_vs,NULL,s_e[1]->vc0_vs_vs,NULL);
	e->vc0_vs_vs  = constructor_operators_tp("vc0","vsA","vsA","H_1_P_PM0",sim->p_s_v,b_e,sim,&ops_tp); // destructed

	switch (sim->method) {
	case METHOD_DG:
		set_operators_tp(&ops_tp,s_e[0]->cv0_vs_vc[0],NULL,s_e[1]->cv0_vs_vc[0],NULL);
		e->cv0_vs_vc[0] = constructor_operators_tp("cv0","vsA","vsA","H_1_P_PM0",sim->p_s_v,b_e,sim,&ops_tp); // destructed

		set_operators_tp(&ops_tp,s_e[0]->cv0_vs_vc[1],NULL,s_e[1]->cv0_vs_vc[1],NULL);
		e->cv0_vs_vc[1] = constructor_operators_tp("cv0","vsA","vsA","H_1_P_PM0",sim->p_s_v,b_e,sim,&ops_tp); // destructed

		set_operators_tp(&ops_tp,s_e[0]->cv0_vs_vs,NULL,s_e[1]->cv0_vs_vs,NULL);
		e->cv0_vs_vs = constructor_operators_tp("cv0","vsA","vsA","H_1_P_PM0",sim->p_s_v,b_e,sim,&ops_tp); // destructed
		break;
	case METHOD_DPG:
		/// \todo Add necessary operators.
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",sim->method);
		break;
	}
}
