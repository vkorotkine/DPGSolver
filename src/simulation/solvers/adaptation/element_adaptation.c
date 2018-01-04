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

#include "element_adaptation.h"

#include "macros.h"
#include "definitions_elements.h"

#include "element_operators.h"
#include "element_operators_tp.h"
#include "multiarray_operator.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

/// \brief \ref constructor_derived_Adaptation_Element constructing the standard operators.
static void constructor_derived_Adaptation_Element_std
	(struct Element* element_ptr, ///< See brief.
	 const struct Simulation* sim ///< See brief.
	);

/// \brief \ref constructor_derived_Adaptation_Element constructing the tensor-product of sub-element operators.
static void constructor_derived_Adaptation_Element_tp
	(struct Element* element_ptr, ///< See brief.
	 const struct Simulation* sim ///< See brief.
	);

/// \brief \ref constructor_derived_Adaptation_Element constructing the common operators.
static void constructor_derived_Adaptation_Element_common
	(struct Element* element_ptr, ///< See brief.
	 const struct Simulation* sim ///< See brief.
	);

// Interface functions ********************************************************************************************** //

void constructor_derived_Adaptation_Element (struct Element* element_ptr, const struct Simulation* sim)
{
	switch (element_ptr->type) {
	case LINE: case TRI: case TET: case PYR:
		constructor_derived_Adaptation_Element_std(element_ptr,sim);
		break;
	case QUAD: case HEX: case WEDGE:
		constructor_derived_Adaptation_Element_tp(element_ptr,sim);
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
	constructor_derived_Adaptation_Element_common(element_ptr,sim);
}

void destructor_derived_Adaptation_Element (struct Element* element_ptr)
{
	struct Adaptation_Element* a_e = (struct Adaptation_Element*) element_ptr;

	destructor_Multiarray_Operator(a_e->cc0_vs_vs);
	destructor_Multiarray_Operator_conditional(a_e->cc0_ff_ff);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void constructor_derived_Adaptation_Element_std (struct Element* element_ptr, const struct Simulation* sim)
{
	struct const_Element* e        = (struct const_Element*) element_ptr;
	struct Adaptation_Element* a_e = (struct Adaptation_Element*) element_ptr;

	a_e->cc0_vs_vs = constructor_operators("cc0","vsA","vsA","H_ALL_P_PM1",e,sim); // destructed
}

static void constructor_derived_Adaptation_Element_tp (struct Element* element_ptr, const struct Simulation* sim)
{
	struct Adaptation_Element* a_e = (struct Adaptation_Element*) element_ptr;

	const struct const_Element* e      = (const struct const_Element*) element_ptr;
	const struct const_Element* se[2]  = { e->sub_element[0], e->sub_element[1], };
	struct Adaptation_Element* a_se[2] = { (struct Adaptation_Element*) se[0], (struct Adaptation_Element*) se[1], };

	struct Operators_TP ops_tp;

	set_operators_tp(&ops_tp,a_se[0]->cc0_vs_vs,NULL,a_se[1]->cc0_vs_vs,NULL);
	a_e->cc0_vs_vs = constructor_operators_tp("cc0","vsA","vsA","H_ALL_P_PM1",e,sim,&ops_tp); // destructed
}

static void constructor_derived_Adaptation_Element_common (struct Element* element_ptr, const struct Simulation* sim)
{
	struct const_Element* e        = (struct const_Element*) element_ptr;
	struct Adaptation_Element* a_e = (struct Adaptation_Element*) element_ptr;

	switch (sim->method) {
	case METHOD_DG:
		break; // Do nothing.
	case METHOD_DPG:
		a_e->cc0_ff_ff = constructor_operators("cc0","ffA","ffA","H_ALL_P_PM1",e,sim); // destructed
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",sim->method);
		break;
	}
}
