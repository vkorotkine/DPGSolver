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

#include "element_solver_dg.h"

#include "macros.h"
#include "definitions_elements.h"

#include "multiarray.h"

#include "element_operators.h"
#include "element_operators_tp.h"
#include "multiarray_operator.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

/// \brief Constructor for a derived \ref DG_Solver_Element using the standard operators.
static void constructor_derived_DG_Solver_Element_std
	(struct Element* element_ptr, ///< Defined for \ref constructor_derived_DG_Solver_Element.
	 const struct Simulation* sim ///< Defined for \ref constructor_derived_DG_Solver_Element.
	);

/// \brief Constructor for a derived \ref DG_Solver_Element using the tensor-product of sub-element operators.
static void constructor_derived_DG_Solver_Element_tp
	(struct Element* element_ptr, ///< Defined for \ref constructor_derived_DG_Solver_Element.
	 const struct Simulation* sim ///< Defined for \ref constructor_derived_DG_Solver_Element.
	);

/// \brief Constructor for the common members of a derived \ref DG_Solver_Element.
static void constructor_derived_DG_Solver_Element_common
	(struct Element* element_ptr, ///< Defined for \ref constructor_derived_DG_Solver_Element.
	 const struct Simulation* sim ///< Defined for \ref constructor_derived_DG_Solver_Element.
	);

// Interface functions ********************************************************************************************** //

void constructor_derived_DG_Solver_Element (struct Element* element_ptr, const struct Simulation* sim)
{
	switch (element_ptr->type) {
	case LINE: case TRI: case TET: case PYR:
		constructor_derived_DG_Solver_Element_std(element_ptr,sim);
		break;
	case QUAD: case HEX: case WEDGE:
		constructor_derived_DG_Solver_Element_tp(element_ptr,sim);
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
	constructor_derived_DG_Solver_Element_common(element_ptr,sim);
}

void destructor_derived_DG_Solver_Element (struct Element* element_ptr)
{
	struct DG_Solver_Element* element = (struct DG_Solver_Element*) element_ptr;

	destructor_Multiarray_Operator(element->cv0_vs_vcs);
	destructor_Multiarray_Operator(element->cv0_vs_vcc);
	destructor_Multiarray_Operator(element->tw1_vs_vcs);
	destructor_Multiarray_Operator(element->tw1_vs_vcc);

	destructor_Multiarray_Operator(element->cv0_vs_fcs);
	destructor_Multiarray_Operator(element->cv0_vs_fcc);
	destructor_Multiarray_Operator(element->tw0_vs_fcs);
	destructor_Multiarray_Operator(element->tw0_vs_fcc);

	struct const_Element* b_e = (struct const_Element*)element_ptr;
	for (int i = 0; i < 2; ++i) {
		struct DG_Solver_Element* f_element = (struct DG_Solver_Element*) b_e->face_element[i];
		if (f_element && f_element->nc_fcs != NULL) {
			destructor_const_Multiarray_Vector_i(f_element->nc_fcs);
			destructor_const_Multiarray_Vector_i(f_element->nc_fcc);
			f_element->nc_fcs = NULL;
		}
	}

	destructor_Multiarray_Operator(element->tw0_vs_vcs);
	destructor_Multiarray_Operator(element->tw0_vs_vcc);

	destructor_const_Multiarray_Vector_d(element->w_vcs);
	destructor_const_Multiarray_Vector_d(element->w_vcc);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void constructor_derived_DG_Solver_Element_std (struct Element* element_ptr, const struct Simulation* sim)
{
	struct DG_Solver_Element* element = (struct DG_Solver_Element*) element_ptr;

	struct const_Element* b_e = (struct const_Element*)element;

	// H_CF, P_PM1 are needed for cv0_vs_vc*, and tw0_vs_vc* operators as they are used to assemble tensor-product
	// operators.
	element->cv0_vs_vcs = constructor_operators("cv0","vsA","vcs","H_CF_P_PM1",sim->p_s_v,b_e,sim); // destructed
	element->cv0_vs_vcc = constructor_operators("cv0","vsA","vcc","H_CF_P_PM1",sim->p_s_v,b_e,sim); // destructed
	element->tw1_vs_vcs = constructor_operators("tw1","vsA","vcs","H_1_P_PM0", sim->p_s_v,b_e,sim); // destructed
	element->tw1_vs_vcc = constructor_operators("tw1","vsA","vcc","H_1_P_PM0", sim->p_s_v,b_e,sim); // destructed

	element->cv0_vs_fcs = constructor_operators("cv0","vsA","fcs","H_CF_P_PM1",sim->p_s_v,b_e,sim); // destructed
	element->cv0_vs_fcc = constructor_operators("cv0","vsA","fcc","H_CF_P_PM1",sim->p_s_v,b_e,sim); // destructed
	element->tw0_vs_fcs = constructor_operators("tw0","vsA","fcs","H_CF_P_PM1",sim->p_s_v,b_e,sim); // destructed
	element->tw0_vs_fcc = constructor_operators("tw0","vsA","fcc","H_CF_P_PM1",sim->p_s_v,b_e,sim); // destructed

	element->tw0_vs_vcs = constructor_operators("tw0","vsA","vcs","H_CF_P_PM1",sim->p_s_v,b_e,sim); // destructed
	element->tw0_vs_vcc = constructor_operators("tw0","vsA","vcc","H_CF_P_PM1",sim->p_s_v,b_e,sim); // destructed
}

static void constructor_derived_DG_Solver_Element_tp (struct Element* element_ptr, const struct Simulation* sim)
{
	struct DG_Solver_Element* element = (struct DG_Solver_Element*) element_ptr;

	const struct const_Element* b_e     = (const struct const_Element*)element;
	const struct const_Element* bs_e[2] = { b_e->sub_element[0], b_e->sub_element[1], };
	struct DG_Solver_Element* s_e[2]    = { (struct DG_Solver_Element*) bs_e[0],
	                                        (struct DG_Solver_Element*) bs_e[1], };

	struct Operators_TP ops_tp;

	set_operators_tp(&ops_tp,s_e[0]->cv0_vs_vcs,NULL,s_e[1]->cv0_vs_vcs,NULL);
	element->cv0_vs_vcs = constructor_operators_tp("cv0","vsA","vcs","H_1_P_PM0",sim->p_s_v,b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->cv0_vs_vcc,NULL,s_e[1]->cv0_vs_vcc,NULL);
	element->cv0_vs_vcc = constructor_operators_tp("cv0","vsA","vcc","H_1_P_PM0",sim->p_s_v,b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->tw0_vs_vcs,s_e[0]->tw1_vs_vcs,s_e[1]->tw0_vs_vcs,s_e[1]->tw1_vs_vcs);
	element->tw1_vs_vcs = constructor_operators_tp("tw1","vsA","vcs","H_1_P_PM0",sim->p_s_v,b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->tw0_vs_vcc,s_e[0]->tw1_vs_vcc,s_e[1]->tw0_vs_vcc,s_e[1]->tw1_vs_vcc);
	element->tw1_vs_vcc = constructor_operators_tp("tw1","vsA","vcc","H_1_P_PM0",sim->p_s_v,b_e,sim,&ops_tp); // destructed


	set_operators_tp(&ops_tp,s_e[0]->cv0_vs_vcs,s_e[0]->cv0_vs_fcs,s_e[1]->cv0_vs_vcs,s_e[1]->cv0_vs_fcs);
	element->cv0_vs_fcs = constructor_operators_tp("cv0","vsA","fcs","H_CF_P_PM1",sim->p_s_v,b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->cv0_vs_vcc,s_e[0]->cv0_vs_fcc,s_e[1]->cv0_vs_vcc,s_e[1]->cv0_vs_fcc);
	element->cv0_vs_fcc = constructor_operators_tp("cv0","vsA","fcc","H_CF_P_PM1",sim->p_s_v,b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->tw0_vs_vcs,s_e[0]->tw0_vs_fcs,s_e[1]->tw0_vs_vcs,s_e[1]->tw0_vs_fcs);
	element->tw0_vs_fcs = constructor_operators_tp("tw0","vsA","fcs","H_CF_P_PM1",sim->p_s_v,b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->tw0_vs_vcc,s_e[0]->tw0_vs_fcc,s_e[1]->tw0_vs_vcc,s_e[1]->tw0_vs_fcc);
	element->tw0_vs_fcc = constructor_operators_tp("tw0","vsA","fcc","H_CF_P_PM1",sim->p_s_v,b_e,sim,&ops_tp); // destructed


	set_operators_tp(&ops_tp,s_e[0]->tw0_vs_vcs,NULL,s_e[1]->tw0_vs_vcs,NULL);
	element->tw0_vs_vcs = constructor_operators_tp("tw0","vsA","vcs","H_1_P_PM0",sim->p_s_v,b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->tw0_vs_vcc,NULL,s_e[1]->tw0_vs_vcc,NULL);
	element->tw0_vs_vcc = constructor_operators_tp("tw0","vsA","vcc","H_1_P_PM0",sim->p_s_v,b_e,sim,&ops_tp); // destructed
}

static void constructor_derived_DG_Solver_Element_common (struct Element* element_ptr, const struct Simulation* sim)
{
	struct const_Element* b_e = (struct const_Element*)element_ptr;
	for (int i = 0; i < 2; ++i) {
		struct DG_Solver_Element* f_element = (struct DG_Solver_Element*) b_e->face_element[i];
		if (f_element && f_element->nc_fcs == NULL) {
			f_element->nc_fcs =
				constructor_operators_nc(i,"fcs","fcs","H_1_P_PM0",sim->p_s_v,b_e,sim); // destructed
			f_element->nc_fcc =
				constructor_operators_nc(i,"fcc","fcc","H_1_P_PM0",sim->p_s_v,b_e,sim); // destructed
		}
	}

	struct DG_Solver_Element* element = (struct DG_Solver_Element*) element_ptr;

	element->w_vcs = constructor_operators_w("vcs","vcs","H_1_P_PM0",sim->p_s_v,b_e,sim); // destructed
	element->w_vcc = constructor_operators_w("vcc","vcc","H_1_P_PM0",sim->p_s_v,b_e,sim); // destructed
}
