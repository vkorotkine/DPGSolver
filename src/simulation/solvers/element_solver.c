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

#include "element_solver.h"

#include <string.h>

#include "macros.h"
#include "definitions_elements.h"

#include "multiarray.h"

#include "computational_elements.h"
#include "element_operators.h"
#include "element_operators_tp.h"
#include "multiarray_operator.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

/// \brief Constructor for a derived \ref Solver_Element using the standard operators.
static void constructor_derived_Solver_Element_std
	(struct Element* element_ptr, ///< Defined for \ref constructor_derived_Solver_Element.
	 const struct Simulation* sim ///< Defined for \ref constructor_derived_Solver_Element.
	);

/// \brief Constructor for a derived \ref Solver_Element using the tensor-product of sub-element operators.
static void constructor_derived_Solver_Element_tp
	(struct Element* element_ptr, ///< Defined for \ref constructor_derived_Solver_Element.
	 const struct Simulation* sim ///< Defined for \ref constructor_derived_Solver_Element.
	);

/// \brief Constructor for the common members of a derived \ref Solver_Element.
static void constructor_derived_Solver_Element_common
	(struct Element* element_ptr, ///< Defined for \ref constructor_derived_Solver_Element.
	 const struct Simulation* sim ///< Defined for \ref constructor_derived_Solver_Element.
	);

// Interface functions ********************************************************************************************** //

void constructor_derived_Solver_Element (struct Element* element_ptr, const struct Simulation* sim)
{
	struct Solver_Element* s_e = (struct Solver_Element*) element_ptr;

	constructor_offset_derived_Element(constructor_derived_Adaptation_Element,
		sizeof(struct Element),element_ptr,(struct Element*)&s_e->a_e,sim); // destructed
	constructor_offset_derived_Element(constructor_derived_Geometry_Element,
		sizeof(struct Element),element_ptr,(struct Element*)&s_e->g_e,sim); // destructed
	constructor_offset_derived_Element(constructor_derived_Plotting_Element,
		sizeof(struct Element),element_ptr,(struct Element*)&s_e->p_e,sim); // destructed
	constructor_offset_derived_Element(constructor_derived_Solution_Element,
		sizeof(struct Element),element_ptr,(struct Element*)&s_e->s_e,sim); // destructed

	if (element_ptr->type == POINT)
		return;

	switch (element_ptr->type) {
	case LINE: case TRI: case TET: case PYR:
		constructor_derived_Solver_Element_std(element_ptr,sim);
		break;
	case QUAD: case HEX: case WEDGE:
		constructor_derived_Solver_Element_tp(element_ptr,sim);
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
	constructor_derived_Solver_Element_common(element_ptr,sim);
}

void destructor_derived_Solver_Element (struct Element* element_ptr)
{
	if (element_ptr->type == POINT)
		return;

	struct Solver_Element* s_e = (struct Solver_Element*) element_ptr;

	destructor_derived_Adaptation_Element((struct Element*)&s_e->a_e);
	destructor_derived_Geometry_Element((struct Element*)&s_e->g_e);
	destructor_derived_Plotting_Element((struct Element*)&s_e->p_e);
	destructor_derived_Solution_Element((struct Element*)&s_e->s_e);

	destructor_Multiarray2_Operator(s_e->cv0_vs_vc);
	destructor_Multiarray2_Operator(s_e->tw1_vt_vc);
	for (int i = 0; i < 2; ++i)
		destructor_const_Multiarray_Vector_d(s_e->w_vc[i]);

	destructor_Multiarray2_Operator(s_e->cv0_vs_fc);
	destructor_Multiarray2_Operator(s_e->tw0_vt_fc);
	const int n_fe = get_number_of_face_elements((struct const_Element*)element_ptr);
	for (int i = 0; i < n_fe; ++i) {
		destructor_const_Multiarray_Vector_i(s_e->nc_fc[0]);
		destructor_const_Multiarray_Vector_i(s_e->nc_fc[1]);
	}

	destructor_Multiarray2_Operator(s_e->cv0_vg_vc);
	destructor_Multiarray2_Operator(s_e->tw0_vt_vc);

	destructor_Multiarray_Operator(s_e->ccSB0_vs_vs);
	destructor_Multiarray_Operator(s_e->ccBS0_vs_vs);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void constructor_derived_Solver_Element_std (struct Element* element_ptr, const struct Simulation* sim)
{
	struct const_Element* e    = (struct const_Element*) element_ptr;
	struct Solver_Element* s_e = (struct Solver_Element*) element_ptr;

	// H_CF, P_PM1 are needed for cv0_vs_vc* and tw0_vt_vc* operators as they are used to assemble tensor-product
	// operators.
	s_e->cv0_vs_vc[0] = constructor_operators("cv0","vsA","vcs","H_CF_P_PM1",e,sim); // destructed
	s_e->cv0_vs_vc[1] = constructor_operators("cv0","vsA","vcc","H_CF_P_PM1",e,sim); // destructed
	s_e->tw1_vt_vc[0] = constructor_operators("tw1","vtA","vcs","H_1_P_PM0", e,sim); // destructed
	s_e->tw1_vt_vc[1] = constructor_operators("tw1","vtA","vcc","H_1_P_PM0", e,sim); // destructed

	s_e->cv0_vs_fc[0] = constructor_operators("cv0","vsA","fcs","H_CF_P_PM1",e,sim); // destructed
	s_e->cv0_vs_fc[1] = constructor_operators("cv0","vsA","fcc","H_CF_P_PM1",e,sim); // destructed
	s_e->tw0_vt_fc[0] = constructor_operators("tw0","vtA","fcs","H_CF_P_PM1",e,sim); // destructed
	s_e->tw0_vt_fc[1] = constructor_operators("tw0","vtA","fcc","H_CF_P_PM1",e,sim); // destructed

	s_e->cv0_vg_vc[0] = constructor_operators("cv0","vgs","vcs","H_1_P_1P",  e,sim); // destructed
	s_e->cv0_vg_vc[1] = constructor_operators("cv0","vgc","vcc","H_1_P_PM0", e,sim); // destructed
	s_e->tw0_vt_vc[0] = constructor_operators("tw0","vtA","vcs","H_CF_P_PM1",e,sim); // destructed
	s_e->tw0_vt_vc[1] = constructor_operators("tw0","vtA","vcc","H_CF_P_PM1",e,sim); // destructed

	s_e->ccSB0_vs_vs = constructor_operators_bt("ccSB0","vsA","vsA","H_1_P_PM0",e,sim); // destructed
	s_e->ccBS0_vs_vs = constructor_operators_bt("ccBS0","vsA","vsA","H_1_P_PM0",e,sim); // destructed
}

static void constructor_derived_Solver_Element_tp (struct Element* element_ptr, const struct Simulation* sim)
{
	struct Solver_Element* s_e = (struct Solver_Element*) element_ptr;

	const struct const_Element* e     = (const struct const_Element*) element_ptr;
	const struct const_Element* se[2] = { e->sub_element[0], e->sub_element[1], };
	struct Solver_Element* s_se[2]    = { (struct Solver_Element*) se[0], (struct Solver_Element*) se[1], };

	struct Operators_TP ops_tp;

	set_operators_tp(&ops_tp,s_se[0]->cv0_vs_vc[0],NULL,s_se[1]->cv0_vs_vc[0],NULL);
	s_e->cv0_vs_vc[0] = constructor_operators_tp("cv0","vsA","vcs","H_1_P_PM0",e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_se[0]->cv0_vs_vc[1],NULL,s_se[1]->cv0_vs_vc[1],NULL);
	s_e->cv0_vs_vc[1] = constructor_operators_tp("cv0","vsA","vcc","H_1_P_PM0",e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_se[0]->tw0_vt_vc[0],s_se[0]->tw1_vt_vc[0],s_se[1]->tw0_vt_vc[0],s_se[1]->tw1_vt_vc[0]);
	s_e->tw1_vt_vc[0] = constructor_operators_tp("tw1","vtA","vcs","H_1_P_PM0",e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_se[0]->tw0_vt_vc[1],s_se[0]->tw1_vt_vc[1],s_se[1]->tw0_vt_vc[1],s_se[1]->tw1_vt_vc[1]);
	s_e->tw1_vt_vc[1] = constructor_operators_tp("tw1","vtA","vcc","H_1_P_PM0",e,sim,&ops_tp); // destructed


	set_operators_tp(&ops_tp,s_se[0]->cv0_vs_vc[0],s_se[0]->cv0_vs_fc[0],s_se[1]->cv0_vs_vc[0],s_se[1]->cv0_vs_fc[0]);
	s_e->cv0_vs_fc[0] = constructor_operators_tp("cv0","vsA","fcs","H_CF_P_PM1",e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_se[0]->cv0_vs_vc[1],s_se[0]->cv0_vs_fc[1],s_se[1]->cv0_vs_vc[1],s_se[1]->cv0_vs_fc[1]);
	s_e->cv0_vs_fc[1] = constructor_operators_tp("cv0","vsA","fcc","H_CF_P_PM1",e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_se[0]->tw0_vt_vc[0],s_se[0]->tw0_vt_fc[0],s_se[1]->tw0_vt_vc[0],s_se[1]->tw0_vt_fc[0]);
	s_e->tw0_vt_fc[0] = constructor_operators_tp("tw0","vsA","fcs","H_CF_P_PM1",e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_se[0]->tw0_vt_vc[1],s_se[0]->tw0_vt_fc[1],s_se[1]->tw0_vt_vc[1],s_se[1]->tw0_vt_fc[1]);
	s_e->tw0_vt_fc[1] = constructor_operators_tp("tw0","vsA","fcc","H_CF_P_PM1",e,sim,&ops_tp); // destructed


	set_operators_tp(&ops_tp,s_se[0]->cv0_vg_vc[0],NULL,s_se[1]->cv0_vg_vc[0],NULL);
	s_e->cv0_vg_vc[0] = constructor_operators_tp("cv0","vgs","vcs","H_1_P_1P",e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_se[0]->cv0_vg_vc[1],NULL,s_se[1]->cv0_vg_vc[1],NULL);
	s_e->cv0_vg_vc[1] = constructor_operators_tp("cv0","vgc","vcc","H_1_P_PM0",e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_se[0]->tw0_vt_vc[0],NULL,s_se[1]->tw0_vt_vc[0],NULL);
	s_e->tw0_vt_vc[0] = constructor_operators_tp("tw0","vsA","vcs","H_1_P_PM0",e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_se[0]->tw0_vt_vc[1],NULL,s_se[1]->tw0_vt_vc[1],NULL);
	s_e->tw0_vt_vc[1] = constructor_operators_tp("tw0","vsA","vcc","H_1_P_PM0",e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_se[0]->ccSB0_vs_vs,NULL,s_se[1]->ccSB0_vs_vs,NULL);
	s_e->ccSB0_vs_vs = constructor_operators_tp("ccSB0","vsA","vsA","H_1_P_PM0",e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_se[0]->ccBS0_vs_vs,NULL,s_se[1]->ccBS0_vs_vs,NULL);
	s_e->ccBS0_vs_vs = constructor_operators_tp("ccBS0","vsA","vsA","H_1_P_PM0",e,sim,&ops_tp); // destructed
}

static void constructor_derived_Solver_Element_common (struct Element* element_ptr, const struct Simulation* sim)
{
	const struct const_Element* e = (struct const_Element*) element_ptr;
	struct Solver_Element* s_e    = (struct Solver_Element*) element_ptr;

	s_e->w_vc[0] = constructor_operators_w("vcs","vcs","H_1_P_PM0",sim->p_s_v,e,sim); // destructed
	s_e->w_vc[1] = constructor_operators_w("vcc","vcc","H_1_P_PM0",sim->p_s_v,e,sim); // destructed

	const int n_fe = get_number_of_face_elements(e);
	for (int i = 0; i < n_fe; ++i) {
		s_e->nc_fc[0] = constructor_operators_nc(i,"fcs","fcs","H_1_P_PM0",sim->p_s_v,e,sim); // destructed
		s_e->nc_fc[1] = constructor_operators_nc(i,"fcc","fcc","H_1_P_PM0",sim->p_s_v,e,sim); // destructed
	}
}
