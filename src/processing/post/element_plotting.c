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
}

void destructor_derived_Plotting_Element (struct Element* element_ptr)
{
	struct Plotting_Element* element = (struct Plotting_Element*) element_ptr;
	destructor_Multiarray_Operator(element->cv0_vgs_vps);
	destructor_Multiarray_Operator(element->cv0_vgc_vpc);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void constructor_derived_Plotting_Element_std (struct Element* element_ptr, const struct Simulation* sim)
{
	struct Plotting_Element* element = (struct Plotting_Element*) element_ptr;

	struct const_Element* b_e = (struct const_Element*)element;

	element->cv0_vgs_vps = constructor_operators("cv0","vgs","vps","H_1_P_1P",sim->p_s_v,b_e,sim); // destructed
	element->cv0_vgc_vpc = constructor_operators("cv0","vgc","vpc","H_1_P_PM0",sim->p_s_v,b_e,sim); // destructed
}

static void constructor_derived_Plotting_Element_tp (struct Element* element_ptr, const struct Simulation* sim)
{
	struct Plotting_Element* element = (struct Plotting_Element*) element_ptr;

	const struct const_Element* b_e     = (const struct const_Element*)element;
	const struct const_Element* bs_e[2] = { b_e->sub_element[0], b_e->sub_element[1], };
	struct Plotting_Element* s_e[2]     = { (struct Plotting_Element*) bs_e[0],
	                                        (struct Plotting_Element*) bs_e[1], };

	struct Operators_TP ops_tp;

	set_operators_tp(&ops_tp,s_e[0]->cv0_vgs_vps,NULL,s_e[1]->cv0_vgs_vps,NULL);
	element->cv0_vgs_vps = constructor_operators_tp("cv0","vgs","vps","H_1_P_1P",sim->p_s_v,b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->cv0_vgc_vpc,NULL,s_e[1]->cv0_vgc_vpc,NULL);
	element->cv0_vgc_vpc = constructor_operators_tp("cv0","vgc","vpc","H_1_P_PM0",sim->p_s_v,b_e,sim,&ops_tp); // destructed
}
