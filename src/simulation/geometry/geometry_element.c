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

#include "geometry_element.h"

#include <assert.h>
#include <string.h>

#include "macros.h"
#include "definitions_intrusive.h"
#include "definitions_element_operators.h"
#include "definitions_elements.h"

#include "multiarray.h"
#include "multiarray_operator.h"

#include "simulation.h"
#include "element_operators.h"
#include "element_operators_tp.h"
#include "const_cast.h"

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

// Interface functions ********************************************************************************************** //

void constructor_derived_Geometry_Element (struct Element* element_ptr, const struct Simulation* sim)
{
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
}

void destructor_derived_Geometry_Element (struct Element* element_ptr)
{
	struct Geometry_Element* element = (struct Geometry_Element*) element_ptr;
	destructor_Multiarray_Operator(element->vc0_vgc_vgc);

	destructor_Multiarray2_Operator(element->cv1_vg_vc);
	destructor_Multiarray2_Operator(element->cv1_vg_vm);
	destructor_Multiarray2_Operator(element->vv0_vm_vc);

	destructor_Multiarray2_Operator(element->cv0_vgs_fc);
	destructor_Multiarray2_Operator(element->cv0_vgc_fc);
	destructor_Multiarray2_Operator(element->vv0_vms_fc);
	destructor_Multiarray2_Operator(element->vv0_vmc_fc);

	destructor_Multiarray2_Operator_conditional(element->cv0_vg_vc);
	destructor_Multiarray2_Operator_conditional(element->cv0_vg_vm);

	destructor_Multiarray_Operator_conditional(element->cv0_vgs_vcc);
	destructor_Multiarray_Operator_conditional(element->cv0_vgc_vcs);
	destructor_Multiarray_Operator_conditional(element->vv0_vms_vcc);
	destructor_Multiarray_Operator_conditional(element->vv0_vmc_vcs);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void constructor_derived_Geometry_Element_std (struct Element* element_ptr, const struct Simulation* sim)
{
	struct Geometry_Element* element = (struct Geometry_Element*) element_ptr;

	struct const_Element* b_e = (struct const_Element*)element;

	element->vc0_vgc_vgc = constructor_operators("vc0","vgc","vgc","H_1_P_PM0",b_e,sim); // destructed

	element->cv1_vg_vc[0] = constructor_operators("cv1","vgs","vcs","H_1_P_1P", b_e,sim); // destructed
	element->cv1_vg_vc[1] = constructor_operators("cv1","vgc","vcc","H_1_P_PM0",b_e,sim); // destructed
	element->cv1_vg_vm[0] = constructor_operators("cv1","vgs","vms","H_1_P_1",  b_e,sim); // destructed
	element->cv1_vg_vm[1] = constructor_operators("cv1","vgc","vmc","H_1_P_PM0",b_e,sim); // destructed
	element->vv0_vm_vc[0] = constructor_operators("vv0","vms","vcs","H_1_P_1P", b_e,sim); // destructed
	element->vv0_vm_vc[1] = constructor_operators("vv0","vmc","vcc","H_1_P_PM0",b_e,sim); // destructed

	element->cv0_vgs_fc[0] = constructor_operators("cv0","vgs","fcs","H_1_P_1P", b_e,sim); // destructed
	element->cv0_vgs_fc[1] = constructor_operators("cv0","vgs","fcc","H_1_P_1P", b_e,sim); // destructed
	element->cv0_vgc_fc[0] = constructor_operators("cv0","vgc","fcs","H_1_P_PM1",b_e,sim); // destructed
	element->cv0_vgc_fc[1] = constructor_operators("cv0","vgc","fcc","H_1_P_PM1",b_e,sim); // destructed

	element->vv0_vms_fc[1] = constructor_operators("vv0","vms","fcs","H_1_P_1P", b_e,sim); // destructed
	element->vv0_vms_fc[0] = constructor_operators("vv0","vms","fcc","H_1_P_1P", b_e,sim); // destructed
	element->vv0_vmc_fc[1] = constructor_operators("vv0","vmc","fcs","H_1_P_PM1",b_e,sim); // destructed
	element->vv0_vmc_fc[0] = constructor_operators("vv0","vmc","fcc","H_1_P_PM1",b_e,sim); // destructed
}

static void constructor_derived_Geometry_Element_tp (struct Element* element_ptr, const struct Simulation* sim)
{
	struct Geometry_Element* element = (struct Geometry_Element*) element_ptr;

	const struct const_Element* b_e     = (const struct const_Element*)element;
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
	}

	set_operators_tp(&ops_tp,s_e[0]->vc0_vgc_vgc,NULL,s_e[1]->vc0_vgc_vgc,NULL);
	element->vc0_vgc_vgc = constructor_operators_tp("vc0","vgc","vgc","H_1_P_PM0",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->cv0_vg_vc[0],s_e[0]->cv1_vg_vc[0],s_e[1]->cv0_vg_vc[0],s_e[1]->cv1_vg_vc[0]);
	element->cv1_vg_vc[0] = constructor_operators_tp("cv1","vgs","vcs","H_1_P_1P",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->cv0_vg_vc[1],s_e[0]->cv1_vg_vc[1],s_e[1]->cv0_vg_vc[1],s_e[1]->cv1_vg_vc[1]);
	element->cv1_vg_vc[1] = constructor_operators_tp("cv1","vgs","vcc","H_1_P_PM0",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->cv0_vg_vm[0],s_e[0]->cv1_vg_vm[0],s_e[1]->cv0_vg_vm[0],s_e[1]->cv1_vg_vm[0]);
	element->cv1_vg_vm[0] = constructor_operators_tp("cv1","vgs","vms","H_1_P_1",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->cv0_vg_vm[1],s_e[0]->cv1_vg_vm[1],s_e[1]->cv0_vg_vm[1],s_e[1]->cv1_vg_vm[1]);
	element->cv1_vg_vm[1] = constructor_operators_tp("cv1","vgs","vmc","H_1_P_PM0",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->vv0_vm_vc[0],NULL,s_e[1]->vv0_vm_vc[0],NULL);
	element->vv0_vm_vc[0] = constructor_operators_tp("vv0","vms","vcs","H_1_P_1P", b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->vv0_vm_vc[1],NULL,s_e[1]->vv0_vm_vc[1],NULL);
	element->vv0_vm_vc[1] = constructor_operators_tp("vv0","vmc","vcc","H_1_P_PM0",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->cv0_vg_vc[0],s_e[0]->cv0_vgs_fc[0],s_e[1]->cv0_vg_vc[0],s_e[1]->cv0_vgs_fc[0]);
	element->cv0_vgs_fc[0] = constructor_operators_tp("cv0","vgs","fcs","H_1_P_1P",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->cv0_vgs_vcc,s_e[0]->cv0_vgs_fc[1],s_e[1]->cv0_vgs_vcc,s_e[1]->cv0_vgs_fc[1]);
	element->cv0_vgs_fc[1] = constructor_operators_tp("cv0","vgs","fcc","H_1_P_1P",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->cv0_vgc_vcs,s_e[0]->cv0_vgc_fc[0],s_e[1]->cv0_vgc_vcs,s_e[1]->cv0_vgc_fc[0]);
	element->cv0_vgc_fc[0] = constructor_operators_tp("cv0","vgc","fcs","H_1_P_PM1",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->cv0_vg_vc[1],s_e[0]->cv0_vgc_fc[1],s_e[1]->cv0_vg_vc[1],s_e[1]->cv0_vgc_fc[1]);
	element->cv0_vgc_fc[1] = constructor_operators_tp("cv0","vgc","fcc","H_1_P_PM1",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->vv0_vm_vc[0],s_e[0]->vv0_vms_fc[0],s_e[1]->vv0_vm_vc[0],s_e[1]->vv0_vms_fc[0]);
	element->vv0_vms_fc[0] = constructor_operators_tp("vv0","vms","fcs","H_1_P_1P",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->vv0_vms_vcc,s_e[0]->vv0_vms_fc[1],s_e[1]->vv0_vms_vcc,s_e[1]->vv0_vms_fc[1]);
	element->vv0_vms_fc[1] = constructor_operators_tp("vv0","vms","fcc","H_1_P_1P",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->vv0_vmc_vcs,s_e[0]->vv0_vmc_fc[0],s_e[1]->vv0_vmc_vcs,s_e[1]->vv0_vmc_fc[0]);
	element->vv0_vmc_fc[0] = constructor_operators_tp("vv0","vmc","fcs","H_1_P_PM1",b_e,sim,&ops_tp); // destructed

	set_operators_tp(&ops_tp,s_e[0]->vv0_vm_vc[1],s_e[0]->vv0_vmc_fc[1],s_e[1]->vv0_vm_vc[1],s_e[1]->vv0_vmc_fc[1]);
	element->vv0_vmc_fc[1] = constructor_operators_tp("vv0","vmc","fcc","H_1_P_PM1",b_e,sim,&ops_tp); // destructed
}
