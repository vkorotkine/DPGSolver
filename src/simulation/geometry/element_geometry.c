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

#include "element_geometry.h"

#include <assert.h>
#include <string.h>

#include "macros.h"
#include "definitions_intrusive.h"
#include "definitions_element_operators.h"
#include "definitions_elements.h"

#include "multiarray.h"

#include "const_cast.h"
#include "element_operators.h"
#include "element_operators_tp.h"
#include "multiarray_operator.h"
#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Constructor for a derived \ref Geometry_Element using the standard operators.
static void constructor_derived_Geometry_Element_std
	(struct Element* element,     ///< Defined for \ref constructor_derived_Geometry_Element.
	 const struct Simulation* sim ///< Defined for \ref constructor_derived_Geometry_Element.
	);

/// \brief Constructor for a derived \ref Geometry_Element using the tensor-product of sub-element operators.
static void constructor_derived_Geometry_Element_tp
	(struct Element* element,     ///< Defined for \ref constructor_derived_Geometry_Element.
	 const struct Simulation* sim ///< Defined for \ref constructor_derived_Geometry_Element.
	);

/// \brief Constructor for the common members of a derived \ref Geometry_Element.
static void constructor_derived_Geometry_Element_common
	(struct Element* element_ptr, ///< Defined for \ref constructor_derived_Geometry_Element.
	 const struct Simulation* sim ///< Defined for \ref constructor_derived_Geometry_Element.
	);

// Interface functions ********************************************************************************************** //

void constructor_derived_Geometry_Element (struct Element* element_ptr, const struct Simulation* sim)
{
	if (element_ptr->type == POINT)
		return;

	switch (element_ptr->type) {
	case LINE: case TRI: case TET: case PYR:
		constructor_derived_Geometry_Element_std(element_ptr,sim);
		break;
	case QUAD: case HEX: case WEDGE:
		constructor_derived_Geometry_Element_tp(element_ptr,sim);
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
	constructor_derived_Geometry_Element_common(element_ptr,sim);
}

void destructor_derived_Geometry_Element (struct Element* element_ptr)
{
	if (element_ptr->type == POINT)
		return;

	struct Geometry_Element* g_e = (struct Geometry_Element*) element_ptr;
	destructor_Multiarray2_Operator(g_e->vv0_vv_vg);
	destructor_Multiarray_Operator(g_e->vc0_vgc_vgc);

	destructor_Multiarray2_Operator(g_e->cv1_vg_vc);
	destructor_Multiarray2_Operator(g_e->cv1_vg_vm);
	destructor_Multiarray2_Operator(g_e->vv0_vm_vc);

	destructor_Multiarray2_Operator(g_e->cv0_vgs_fc);
	destructor_Multiarray2_Operator(g_e->cv0_vgc_fc);
	destructor_Multiarray2_Operator(g_e->vv0_vms_fc);
	destructor_Multiarray2_Operator(g_e->vv0_vmc_fc);

	destructor_Multiarray_Operator_conditional(g_e->vv0_fv_vgc);
	destructor_Multiarray_Operator_conditional(g_e->vv0_vv_fgc);
	destructor_Multiarray_Operator_conditional(g_e->vv0_vv_fv);
	destructor_Multiarray_Operator_conditional(g_e->vv0_vv_fcc);
	destructor_Multiarray_Operator_conditional(g_e->vv0_fv_fgc);
	destructor_Multiarray_Operator_conditional(g_e->vv0_fgc_vgc);
	destructor_Multiarray_Operator(g_e->vv0_vgc_fgc);

	destructor_Multiarray_Operator(g_e->cv0_vgc_fgc);
	destructor_Multiarray_Operator(g_e->cc0_vgc_fgc);
	destructor_Multiarray_Operator_conditional(g_e->vv0_fgc_fgc);
	destructor_Multiarray_Operator_conditional(g_e->vc0_fgc_fgc);
	destructor_Multiarray2_Operator(g_e->cv0_vg_vv);

	destructor_Multiarray_Operator(g_e->cv0_vgc_fis);
	destructor_Multiarray_Operator(g_e->vc0_fis_fgc);

	const int n_fe = get_number_of_face_elements((struct const_Element*)element_ptr);
	for (int i = 0; i < n_fe; ++i) {
		destructor_conditional_const_Multiarray_Vector_i(g_e->nc_fg[0]);
		destructor_conditional_const_Multiarray_Vector_i(g_e->nc_fg[1]);
	}

	destructor_Multiarray2_Operator_conditional(g_e->cv0_vg_vc);
	destructor_Multiarray2_Operator_conditional(g_e->cv0_vg_vm);

	destructor_Multiarray_Operator_conditional(g_e->cv0_vgs_vcc);
	destructor_Multiarray_Operator_conditional(g_e->cv0_vgc_vcs);
	destructor_Multiarray_Operator_conditional(g_e->vv0_vms_vcc);
	destructor_Multiarray_Operator_conditional(g_e->vv0_vmc_vcs);

	destructor_Multiarray_Operator_conditional(g_e->vv0_vgc_vgc);
	destructor_Multiarray_Operator_conditional(g_e->cv0_vgc_vgc);
	destructor_Multiarray_Operator_conditional(g_e->cc0_vgc_vgc);

	destructor_Multiarray_Operator_conditional(g_e->cv0_vgc_vis);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void constructor_derived_Geometry_Element_std (struct Element* element_ptr, const struct Simulation* sim)
{
	struct Geometry_Element* g_e = (struct Geometry_Element*) element_ptr;

	struct const_Element* b_e = (struct const_Element*)g_e;

	g_e->vv0_vv_vg[0] = constructor_operators("vv0","vvA","vgs","H_1_P_1P", b_e,sim); // destructed
	g_e->vv0_vv_vg[1] = constructor_operators("vv0","vvA","vgc","H_1_P_1P", b_e,sim); // destructed
	g_e->vc0_vgc_vgc  = constructor_operators("vc0","vgc","vgc","H_1_P_PM0",b_e,sim); // destructed

	g_e->cv1_vg_vc[0] = constructor_operators("cv1","vgs","vcs","H_1_P_1P", b_e,sim); // destructed
	g_e->cv1_vg_vc[1] = constructor_operators("cv1","vgc","vcc","H_1_P_PM0",b_e,sim); // destructed
	g_e->cv1_vg_vm[0] = constructor_operators("cv1","vgs","vms","H_1_P_1",  b_e,sim); // destructed
	g_e->cv1_vg_vm[1] = constructor_operators("cv1","vgc","vmc","H_1_P_PM0",b_e,sim); // destructed
	g_e->vv0_vm_vc[0] = constructor_operators("vv0","vms","vcs","H_1_P_1P", b_e,sim); // destructed
	g_e->vv0_vm_vc[1] = constructor_operators("vv0","vmc","vcc","H_1_P_PM1",b_e,sim); // destructed

	g_e->cv0_vgs_fc[0] = constructor_operators("cv0","vgs","fcs","H_1_P_1P", b_e,sim); // destructed
	g_e->cv0_vgs_fc[1] = constructor_operators("cv0","vgs","fcc","H_1_P_1P", b_e,sim); // destructed
	g_e->cv0_vgc_fc[0] = constructor_operators("cv0","vgc","fcs","H_1_P_PM1",b_e,sim); // destructed
	g_e->cv0_vgc_fc[1] = constructor_operators("cv0","vgc","fcc","H_1_P_PM1",b_e,sim); // destructed

	g_e->vv0_vms_fc[0] = constructor_operators("vv0","vms","fcs","H_1_P_1P", b_e,sim); // destructed
	g_e->vv0_vms_fc[1] = constructor_operators("vv0","vms","fcc","H_1_P_1P", b_e,sim); // destructed
	g_e->vv0_vmc_fc[0] = constructor_operators("vv0","vmc","fcs","H_1_P_PM1",b_e,sim); // destructed
	g_e->vv0_vmc_fc[1] = constructor_operators("vv0","vmc","fcc","H_1_P_PM1",b_e,sim); // destructed

	g_e->vv0_vgc_fgc = constructor_operators("vv0","vgc","fgc","H_1_P_PM0", b_e,sim); // destructed

	g_e->cv0_vgc_fgc = constructor_operators("cv0","vgc","fgc","H_CF_P_ALL",b_e,sim); // destructed
	g_e->cc0_vgc_fgc = constructor_operators("cc0","vgc","fgc","H_CF_P_PM1",b_e,sim); // destructed

	g_e->cv0_vg_vv[0] = constructor_operators("cv0","vgs","vvA","H_1_P_1",  b_e,sim); // destructed
	g_e->cv0_vg_vv[1] = constructor_operators("cv0","vgc","vvA","H_1_P_P1", b_e,sim); // destructed

	g_e->cv0_vgc_fis = constructor_operators("cv0","vgc","fis","H_1_P_P1",b_e,sim); // destructed
}

static void constructor_derived_Geometry_Element_tp (struct Element* element_ptr, const struct Simulation* sim)
{
	struct Geometry_Element* g_e = (struct Geometry_Element*) element_ptr;

	const struct const_Element* b_e     = (const struct const_Element*)g_e;
	const struct const_Element* bs_e[2] = { b_e->sub_element[0], b_e->sub_element[1], };
	struct Geometry_Element* s_e[2]     = { (struct Geometry_Element*) bs_e[0],
	                                        (struct Geometry_Element*) bs_e[1], };

	struct Operators_TP ops_tp;

	// tensor-product sub-operators
	for (int i = 0; i < 2; ++i) {
		if (s_e[i]->cv0_vg_vc[0] != NULL)
			continue;

		s_e[i]->cv0_vg_vc[0] = constructor_operators("cv0","vgs","vcs","H_1_P_1P", bs_e[i],sim); // destructed
		s_e[i]->cv0_vg_vc[1] = constructor_operators("cv0","vgc","vcc","H_1_P_PM1",bs_e[i],sim); // destructed
		s_e[i]->cv0_vg_vm[0] = constructor_operators("cv0","vgs","vms","H_1_P_1",  bs_e[i],sim); // destructed
		s_e[i]->cv0_vg_vm[1] = constructor_operators("cv0","vgc","vmc","H_1_P_PM0",bs_e[i],sim); // destructed

		s_e[i]->cv0_vgs_vcc = constructor_operators("cv0","vgs","vcc","H_1_P_1P", bs_e[i],sim); // destructed
		s_e[i]->cv0_vgc_vcs = constructor_operators("cv0","vgc","vcs","H_1_P_PM1",bs_e[i],sim); // destructed
		s_e[i]->vv0_vms_vcc = constructor_operators("vv0","vms","vcc","H_1_P_1P", bs_e[i],sim); // destructed
		s_e[i]->vv0_vmc_vcs = constructor_operators("vv0","vmc","vcs","H_1_P_PM1",bs_e[i],sim); // destructed

		s_e[i]->vv0_vgc_vgc = constructor_operators("vv0","vgc","vgc","H_1_P_PM0",bs_e[i],sim); // destructed

		s_e[i]->cv0_vgc_vgc = constructor_operators("cv0","vgc","vgc","H_CF_P_ALL",bs_e[i],sim); // destructed
		s_e[i]->cc0_vgc_vgc = constructor_operators("cc0","vgc","vgc","H_CF_P_PM1",bs_e[i],sim); // destructed

		s_e[i]->cv0_vgc_vis = constructor_operators("cv0","vgc","vis","H_1_P_P1",bs_e[i],sim); // destructed
	}

	set_operators_tp(&ops_tp,s_e[0]->vv0_vv_vg[0],NULL,s_e[1]->vv0_vv_vg[0],NULL);
	g_e->vv0_vv_vg[0] = constructor_operators_tp("vv0","vvA","vgs","H_1_P_1P",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->vv0_vv_vg[1],NULL,s_e[1]->vv0_vv_vg[1],NULL);
	g_e->vv0_vv_vg[1] = constructor_operators_tp("vv0","vvA","vgc","H_1_P_1P",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->vc0_vgc_vgc,NULL,s_e[1]->vc0_vgc_vgc,NULL);
	g_e->vc0_vgc_vgc = constructor_operators_tp("vc0","vgc","vgc","H_1_P_PM0",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->cv0_vg_vc[0],s_e[0]->cv1_vg_vc[0],s_e[1]->cv0_vg_vc[0],s_e[1]->cv1_vg_vc[0]);
	g_e->cv1_vg_vc[0] = constructor_operators_tp("cv1","vgs","vcs","H_1_P_1P",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->cv0_vg_vc[1],s_e[0]->cv1_vg_vc[1],s_e[1]->cv0_vg_vc[1],s_e[1]->cv1_vg_vc[1]);
	g_e->cv1_vg_vc[1] = constructor_operators_tp("cv1","vgs","vcc","H_1_P_PM0",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->cv0_vg_vm[0],s_e[0]->cv1_vg_vm[0],s_e[1]->cv0_vg_vm[0],s_e[1]->cv1_vg_vm[0]);
	g_e->cv1_vg_vm[0] = constructor_operators_tp("cv1","vgs","vms","H_1_P_1",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->cv0_vg_vm[1],s_e[0]->cv1_vg_vm[1],s_e[1]->cv0_vg_vm[1],s_e[1]->cv1_vg_vm[1]);
	g_e->cv1_vg_vm[1] = constructor_operators_tp("cv1","vgs","vmc","H_1_P_PM0",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->vv0_vm_vc[0],NULL,s_e[1]->vv0_vm_vc[0],NULL);
	g_e->vv0_vm_vc[0] = constructor_operators_tp("vv0","vms","vcs","H_1_P_1P", b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->vv0_vm_vc[1],NULL,s_e[1]->vv0_vm_vc[1],NULL);
	g_e->vv0_vm_vc[1] = constructor_operators_tp("vv0","vmc","vcc","H_1_P_PM0",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->cv0_vg_vc[0],s_e[0]->cv0_vgs_fc[0],s_e[1]->cv0_vg_vc[0],s_e[1]->cv0_vgs_fc[0]);
	g_e->cv0_vgs_fc[0] = constructor_operators_tp("cv0","vgs","fcs","H_1_P_1P",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->cv0_vgs_vcc,s_e[0]->cv0_vgs_fc[1],s_e[1]->cv0_vgs_vcc,s_e[1]->cv0_vgs_fc[1]);
	g_e->cv0_vgs_fc[1] = constructor_operators_tp("cv0","vgs","fcc","H_1_P_1P",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->cv0_vgc_vcs,s_e[0]->cv0_vgc_fc[0],s_e[1]->cv0_vgc_vcs,s_e[1]->cv0_vgc_fc[0]);
	g_e->cv0_vgc_fc[0] = constructor_operators_tp("cv0","vgc","fcs","H_1_P_PM1",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->cv0_vg_vc[1],s_e[0]->cv0_vgc_fc[1],s_e[1]->cv0_vg_vc[1],s_e[1]->cv0_vgc_fc[1]);
	g_e->cv0_vgc_fc[1] = constructor_operators_tp("cv0","vgc","fcc","H_1_P_PM1",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->vv0_vm_vc[0],s_e[0]->vv0_vms_fc[0],s_e[1]->vv0_vm_vc[0],s_e[1]->vv0_vms_fc[0]);
	g_e->vv0_vms_fc[0] = constructor_operators_tp("vv0","vms","fcs","H_1_P_1P",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->vv0_vms_vcc,s_e[0]->vv0_vms_fc[1],s_e[1]->vv0_vms_vcc,s_e[1]->vv0_vms_fc[1]);
	g_e->vv0_vms_fc[1] = constructor_operators_tp("vv0","vms","fcc","H_1_P_1P",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->vv0_vmc_vcs,s_e[0]->vv0_vmc_fc[0],s_e[1]->vv0_vmc_vcs,s_e[1]->vv0_vmc_fc[0]);
	g_e->vv0_vmc_fc[0] = constructor_operators_tp("vv0","vmc","fcs","H_1_P_PM1",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->vv0_vm_vc[1],s_e[0]->vv0_vmc_fc[1],s_e[1]->vv0_vm_vc[1],s_e[1]->vv0_vmc_fc[1]);
	g_e->vv0_vmc_fc[1] = constructor_operators_tp("vv0","vmc","fcc","H_1_P_PM1",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->vv0_vgc_vgc,s_e[0]->vv0_vgc_fgc,s_e[1]->vv0_vgc_vgc,s_e[1]->vv0_vgc_fgc);
	g_e->vv0_vgc_fgc = constructor_operators_tp("vv0","vgc","fgc","H_1_P_PM0",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->cv0_vgc_vgc,s_e[0]->cv0_vgc_fgc,s_e[1]->cv0_vgc_vgc,s_e[1]->cv0_vgc_fgc);
	g_e->cv0_vgc_fgc = constructor_operators_tp("cv0","vgc","fgc","H_CF_P_ALL",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->cc0_vgc_vgc,s_e[0]->cc0_vgc_fgc,s_e[1]->cc0_vgc_vgc,s_e[1]->cc0_vgc_fgc);
	g_e->cc0_vgc_fgc = constructor_operators_tp("cc0","vgc","fgc","H_CF_P_PM1",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->cv0_vg_vv[0],NULL,s_e[1]->cv0_vg_vv[0],NULL);
	g_e->cv0_vg_vv[0] = constructor_operators_tp("cv0","vgs","vvA","H_1_P_1",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->cv0_vg_vv[1],NULL,s_e[1]->cv0_vg_vv[1],NULL);
	g_e->cv0_vg_vv[1] = constructor_operators_tp("cv0","vgc","vvA","H_1_P_P1",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->cv0_vgc_vis,s_e[0]->cv0_vgc_fis,s_e[1]->cv0_vgc_vis,s_e[1]->cv0_vgc_fis);
	g_e->cv0_vgc_fis = constructor_operators_tp("cv0","vgc","fis","H_1_P_P1",b_e,sim,&ops_tp); // destructed
}

static void constructor_derived_Geometry_Element_common (struct Element* element_ptr, const struct Simulation* sim)
{
	struct Geometry_Element*const g_e = (struct Geometry_Element*) element_ptr;
	struct const_Element*const e      = (struct const_Element*)g_e;

	if (e->d > 1) {
		g_e->vv0_vv_fgc = constructor_operators("vv0","vvA","fgc","H_1_P_1P",e,sim); // destructed
		g_e->vv0_vv_fv  = constructor_operators("vv0","vvA","fvA","H_1_P_1", e,sim); // destructed
		g_e->vv0_vv_fcc = constructor_operators("vv0","vvA","fcc","H_1_P_1P",e,sim); // destructed

		g_e->vv0_fv_vgc  = constructor_operators("vv0","fvA","vgc","H_1_P_1P", e,sim); // destructed
		g_e->vv0_fv_fgc  = constructor_operators("vv0","fvA","fgc","H_1_P_1P", e,sim); // destructed
		g_e->vv0_fgc_vgc = constructor_operators("vv0","fgc","vgc","H_1_P_PM0",e,sim); // destructed

		g_e->vv0_fgc_fgc = constructor_operators("vv0","fgc","fgc","H_1_P_ALL",e,sim); // destructed
		g_e->vc0_fgc_fgc = constructor_operators("vc0","fgc","fgc","H_1_P_PM0",e,sim); // destructed

		g_e->vc0_fis_fgc = constructor_operators("vc0","fis","fgc","H_1_P_1P",e,sim); // destructed

		g_e->nc_fg[0] = constructor_operators_nc("fgs","fgs","H_1_P_PM0",e,sim); // destructed
		g_e->nc_fg[1] = constructor_operators_nc("fgc","fgc","H_1_P_PM0",e,sim); // destructed
	} else if (e->d > 2) {
		EXIT_ADD_SUPPORT;
	}
}
