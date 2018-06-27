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

#include "element_solver_opg.h"

#include "macros.h"
#include "definitions_elements.h"

#include "multiarray.h"

#include "element_operators.h"
#include "element_operators_tp.h"
#include "multiarray_operator.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

/// \brief \ref constructor_derived_OPG_Solver_Element constructing the standard operators.
static void constructor_derived_OPG_Solver_Element_std
	(struct Element* element_ptr, ///< See brief.
	 const struct Simulation* sim ///< See brief.
	);

/// \brief \ref constructor_derived_OPG_Solver_Element constructing the tensor-product of sub-element operators.
static void constructor_derived_OPG_Solver_Element_tp
	(struct Element* element_ptr, ///< See brief.
	 const struct Simulation* sim ///< See brief.
	);

/// \brief \ref constructor_derived_OPG_Solver_Element constructing the common operators.
static void constructor_derived_OPG_Solver_Element_common
	(struct Element* element_ptr, ///< See brief.
	 const struct Simulation* sim ///< See brief.
	);

// Interface functions ********************************************************************************************** //

void constructor_derived_OPG_Solver_Element (struct Element* element_ptr, const struct Simulation* sim)
{
	if (element_ptr->type == POINT)
		return;

	switch (element_ptr->type) {
	case LINE: case TRI: case TET: case PYR:
		constructor_derived_OPG_Solver_Element_std(element_ptr,sim);
		break;
	case QUAD: case HEX: case WEDGE:
		constructor_derived_OPG_Solver_Element_tp(element_ptr,sim);
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
	constructor_derived_OPG_Solver_Element_common(element_ptr,sim);
}

void destructor_derived_OPG_Solver_Element (struct Element* element_ptr)
{
	if (element_ptr->type == POINT)
		return;

	struct OPG_Solver_Element* opg_s_e = (struct OPG_Solver_Element*) element_ptr;

	destructor_Multiarray2_Operator(opg_s_e->cv0_vt_vc);
	destructor_Multiarray2_Operator(opg_s_e->cv1_vt_vc);
	destructor_Multiarray_Operator(opg_s_e->vc0_vs_vs);
	destructor_Multiarray_Operator_conditional(opg_s_e->cv0_vt_vs);
	destructor_Multiarray_Operator(opg_s_e->cv1_vt_vs);
	destructor_Multiarray2_Operator(opg_s_e->cv0_vg_vt);

	destructor_Multiarray2_Operator(opg_s_e->cv0_vt_fc);
	destructor_Multiarray2_Operator(opg_s_e->cv1_vt_fc);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void constructor_derived_OPG_Solver_Element_std (struct Element* element_ptr, const struct Simulation* sim)
{
	struct const_Element* e = (struct const_Element*) element_ptr;
	struct OPG_Solver_Element* opg_s_e = (struct OPG_Solver_Element*) element_ptr;

	opg_s_e->cv0_vt_vc[0] = constructor_operators("cv0","vtA","vcs","H_CF_P_PM1",e,sim); // destructed
	opg_s_e->cv0_vt_vc[1] = constructor_operators("cv0","vtA","vcc","H_CF_P_PM1",e,sim); // destructed
	opg_s_e->cv1_vt_vc[0] = constructor_operators("cv1","vtA","vcs","H_1_P_PM0",e,sim); // destructed
	opg_s_e->cv1_vt_vc[1] = constructor_operators("cv1","vtA","vcc","H_1_P_PM0",e,sim); // destructed
	opg_s_e->vc0_vs_vs    = constructor_operators("vc0","vsA","vsA","H_1_P_PM0",e,sim); // destructed
	opg_s_e->cv0_vt_vs    = constructor_operators("cv0","vtA","vsA","H_1_P_PM0",e,sim); // destructed
	opg_s_e->cv1_vt_vs    = constructor_operators("cv1","vtA","vsA","H_1_P_PM0",e,sim); // destructed
	opg_s_e->cv0_vg_vt[0] = constructor_operators("cv0","vgs","vtA","H_1_P_1P", e,sim); // destructed
	opg_s_e->cv0_vg_vt[1] = constructor_operators("cv0","vgc","vtA","H_1_P_PM0",e,sim); // destructed

	opg_s_e->cv0_vt_fc[0] = constructor_operators("cv0","vtA","fcs","H_CF_P_PM1",e,sim); // destructed
	opg_s_e->cv0_vt_fc[1] = constructor_operators("cv0","vtA","fcc","H_CF_P_PM1",e,sim); // destructed
}

static void constructor_derived_OPG_Solver_Element_tp (struct Element* element_ptr, const struct Simulation* sim)
{
	struct OPG_Solver_Element* opg_s_e = (struct OPG_Solver_Element*) element_ptr;

	const struct const_Element* e          = (const struct const_Element*) element_ptr;
	const struct const_Element* se[2]      = { e->sub_element[0], e->sub_element[1], };
	struct OPG_Solver_Element* opg_s_se[2] = { (struct OPG_Solver_Element*) se[0],
	                                           (struct OPG_Solver_Element*) se[1], };

	struct Operators_TP ops_tp;

	set_operators_tp(&ops_tp,opg_s_se[0]->cv0_vt_vc[0],NULL,opg_s_se[1]->cv0_vt_vc[0],NULL);
	opg_s_e->cv0_vt_vc[0] = constructor_operators_tp("cv0","vtA","vcs","H_1_P_PM0",e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,opg_s_se[0]->cv0_vt_vc[1],NULL,opg_s_se[1]->cv0_vt_vc[1],NULL);
	opg_s_e->cv0_vt_vc[1] = constructor_operators_tp("cv0","vtA","vcc","H_1_P_PM0",e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,opg_s_se[0]->cv0_vt_vc[0],opg_s_se[0]->cv1_vt_vc[0],
	                         opg_s_se[1]->cv0_vt_vc[0],opg_s_se[1]->cv1_vt_vc[0]);
	opg_s_e->cv1_vt_vc[0] = constructor_operators_tp("cv1","vtA","vcs","H_1_P_PM0",e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,opg_s_se[0]->cv0_vt_vc[1],opg_s_se[0]->cv1_vt_vc[1],
	                         opg_s_se[1]->cv0_vt_vc[1],opg_s_se[1]->cv1_vt_vc[1]);
	opg_s_e->cv1_vt_vc[1] = constructor_operators_tp("cv1","vtA","vcc","H_1_P_PM0",e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,opg_s_se[0]->vc0_vs_vs,NULL,opg_s_se[1]->vc0_vs_vs,NULL);
	opg_s_e->vc0_vs_vs = constructor_operators_tp("vc0","vsA","vsA","H_1_P_PM0",e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,opg_s_se[0]->cv0_vt_vs,opg_s_se[0]->cv1_vt_vs,
			         opg_s_se[1]->cv0_vt_vs,opg_s_se[1]->cv1_vt_vs);
	opg_s_e->cv1_vt_vs = constructor_operators_tp("cv1","vtA","vsA","H_1_P_PM0",e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,opg_s_se[0]->cv0_vg_vt[0],NULL,opg_s_se[1]->cv0_vg_vt[0],NULL);
	opg_s_e->cv0_vg_vt[0] = constructor_operators_tp("cv0","vgs","vtA","H_1_P_1P",e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,opg_s_se[0]->cv0_vg_vt[1],NULL,opg_s_se[1]->cv0_vg_vt[1],NULL);
	opg_s_e->cv0_vg_vt[1] = constructor_operators_tp("cv0","vgc","vtA","H_1_P_PM0",e,sim,&ops_tp); // destructed


	set_operators_tp(&ops_tp,opg_s_se[0]->cv0_vt_vc[0],opg_s_se[0]->cv0_vt_fc[0],
	                         opg_s_se[1]->cv0_vt_vc[0],opg_s_se[1]->cv0_vt_fc[0]);
	opg_s_e->cv0_vt_fc[0] = constructor_operators_tp("cv0","vsA","fcs","H_CF_P_PM1",e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,opg_s_se[0]->cv0_vt_vc[1],opg_s_se[0]->cv0_vt_fc[1],
	                         opg_s_se[1]->cv0_vt_vc[1],opg_s_se[1]->cv0_vt_fc[1]);
	opg_s_e->cv0_vt_fc[1] = constructor_operators_tp("cv0","vsA","fcc","H_CF_P_PM1",e,sim,&ops_tp); // destructed
}

static void constructor_derived_OPG_Solver_Element_common (struct Element* element_ptr, const struct Simulation* sim)
{
	struct const_Element*const e = (struct const_Element*) element_ptr;
	struct OPG_Solver_Element*const opg_s_e = (struct OPG_Solver_Element*) element_ptr;

	// Can be moved to the standard functions once the functionality for 8 tensor-product operator inputs becomes
	// supported.
	opg_s_e->cv1_vt_fc[0] = constructor_operators("cv1","vtA","fcs","H_CF_P_PM1",e,sim); // destructed
	opg_s_e->cv1_vt_fc[1] = constructor_operators("cv1","vtA","fcc","H_CF_P_PM1",e,sim); // destructed
}
