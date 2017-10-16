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
	destructor_Multiarray_Operator(element->cv0_vgs_vss);
	destructor_Multiarray_Operator(element->cv0_vgc_vsc);
	if (element->cv0_vgs_vrs != NULL) {
		destructor_Multiarray_Operator(element->cv0_vgs_vrs);
		destructor_Multiarray_Operator(element->cv0_vgc_vrc);
	}

	if (element->cv0_vgs_vrs != NULL) {
		destructor_Multiarray_Operator(element->cv0_vgs_fss);
		destructor_Multiarray_Operator(element->cv0_vgc_fsc);
	}

	if (element->cv0_vgs_vrs != NULL) {
		destructor_Multiarray_Operator(element->cv0_vgs_frs);
		destructor_Multiarray_Operator(element->cv0_vgc_frc);
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void constructor_derived_Solution_Element_std (struct Element* element_ptr, const struct Simulation* sim)
{
	struct Solution_Element* e = (struct Solution_Element*) element_ptr;

	struct const_Element* b_e = (struct const_Element*)e;

	switch (sim->method) {
	case METHOD_DG:
		e->cv0_vgs_vss = constructor_operators("cv0","vgs","vss","H_1_P_1P", sim->p_s_v,b_e,sim); // destructed
		e->cv0_vgc_vsc = constructor_operators("cv0","vgc","vsc","H_1_P_PM0",sim->p_s_v,b_e,sim); // destructed
		e->cv0_vgs_vrs = NULL;
		e->cv0_vgc_vrc = NULL;
		e->cv0_vgs_fss = NULL;
		e->cv0_vgc_fsc = NULL;
		e->cv0_vgs_frs = NULL;
		e->cv0_vgc_frc = NULL;
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

	switch (sim->method) {
	case METHOD_DG:
		set_operators_tp(&ops_tp,s_e[0]->cv0_vgs_vss,NULL,s_e[1]->cv0_vgs_vss,NULL);
		e->cv0_vgs_vss = constructor_operators_tp("cv0","vgs","vss","H_1_P_1P", sim->p_s_v,b_e,sim,&ops_tp); // destructed

		set_operators_tp(&ops_tp,s_e[0]->cv0_vgc_vsc,NULL,s_e[1]->cv0_vgc_vsc,NULL);
		e->cv0_vgc_vsc = constructor_operators_tp("cv0","vgc","vsc","H_1_P_PM0",sim->p_s_v,b_e,sim,&ops_tp); // destructed

		e->cv0_vgs_vrs = NULL;
		e->cv0_vgc_vrc = NULL;
		e->cv0_vgs_fss = NULL;
		e->cv0_vgc_fsc = NULL;
		e->cv0_vgs_frs = NULL;
		e->cv0_vgc_frc = NULL;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",sim->method);
		break;
	}
}
