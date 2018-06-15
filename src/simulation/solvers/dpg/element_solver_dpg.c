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

#include "element_solver_dpg.h"

#include "macros.h"
#include "definitions_elements.h"

#include "multiarray.h"

#include "element_operators.h"
#include "element_operators_tp.h"
#include "multiarray_operator.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

/// \brief \ref constructor_derived_DPG_Solver_Element constructing the standard operators.
static void constructor_derived_DPG_Solver_Element_std
	(struct Element* element_ptr, ///< See brief.
	 const struct Simulation* sim ///< See brief.
	);

/// \brief \ref constructor_derived_DPG_Solver_Element constructing the tensor-product of sub-element operators.
static void constructor_derived_DPG_Solver_Element_tp
	(struct Element* element_ptr, ///< See brief.
	 const struct Simulation* sim ///< See brief.
	);

/// \brief \ref constructor_derived_DPG_Solver_Element constructing the common operators.
static void constructor_derived_DPG_Solver_Element_common
	(struct Element* element_ptr, ///< See brief.
	 const struct Simulation* sim ///< See brief.
	);

// Interface functions ********************************************************************************************** //

void constructor_derived_DPG_Solver_Element (struct Element* element_ptr, const struct Simulation* sim)
{
	if (element_ptr->type == POINT)
		return;

	switch (element_ptr->type) {
	case LINE: case TRI: case TET: case PYR:
		constructor_derived_DPG_Solver_Element_std(element_ptr,sim);
		break;
	case QUAD: case HEX: case WEDGE:
		constructor_derived_DPG_Solver_Element_tp(element_ptr,sim);
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
	constructor_derived_DPG_Solver_Element_common(element_ptr,sim);
}

void destructor_derived_DPG_Solver_Element (struct Element* element_ptr)
{
	if (element_ptr->type == POINT)
		return;

	struct DPG_Solver_Element* dpg_s_e = (struct DPG_Solver_Element*) element_ptr;

	destructor_Multiarray2_Operator(dpg_s_e->cvcv0_vs_vc);
	destructor_Multiarray2_Operator(dpg_s_e->cvcv1_vt_vc);

	destructor_const_Multiarray_Vector_d(dpg_s_e->ones_coef_vt);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void constructor_derived_DPG_Solver_Element_std (struct Element* element_ptr, const struct Simulation* sim)
{
	UNUSED(element_ptr); UNUSED(sim);
}

static void constructor_derived_DPG_Solver_Element_tp (struct Element* element_ptr, const struct Simulation* sim)
{
	UNUSED(element_ptr); UNUSED(sim);
}

static void constructor_derived_DPG_Solver_Element_common (struct Element* element_ptr, const struct Simulation* sim)
{
	struct const_Element* e = (struct const_Element*) element_ptr;
	struct DPG_Solver_Element* dpg_s_e = (struct DPG_Solver_Element*) element_ptr;

	struct Solver_Element* s_e = (struct Solver_Element*) element_ptr;

	dpg_s_e->cvcv0_vs_vc[0] = constructor_operators_tens3(s_e->cv0_vs_vc[0],s_e->cv0_vs_vc[0]); // destructed
	dpg_s_e->cvcv0_vs_vc[1] = constructor_operators_tens3(s_e->cv0_vs_vc[1],s_e->cv0_vs_vc[1]); // destructed

	dpg_s_e->cvcv1_vt_vc[0] = constructor_operators_tens3(s_e->cv0_vs_vc[0],s_e->cv1_vt_vc[0]); // destructed
	dpg_s_e->cvcv1_vt_vc[1] = constructor_operators_tens3(s_e->cv0_vs_vc[1],s_e->cv1_vt_vc[1]); // destructed

	dpg_s_e->ones_coef_vt = constructor_operators_ones_coef("vtA",e,sim); // destructed
}
