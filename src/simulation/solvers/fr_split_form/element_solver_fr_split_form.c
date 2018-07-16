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

#include "element_solver_fr_split_form.h"

#include "macros.h"
#include "definitions_elements.h"

#include "matrix.h"
#include "multiarray.h"

#include "volume_solver.h"

#include "compute_volume_rlhs.h"
#include "element_operators.h"
#include "element_operators_tp.h"
#include "operator.h"
#include "multiarray_operator.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

/// \brief Constructor for a derived \ref FRSF_Solver_Element using the standard operators.
static void constructor_derived_FRSF_Solver_Element_std
	(struct Element* element_ptr, ///< Defined for \ref constructor_derived_FRSF_Solver_Element.
	 const struct Simulation* sim ///< Defined for \ref constructor_derived_FRSF_Solver_Element.
	);

// Interface functions ********************************************************************************************** //

void constructor_derived_FRSF_Solver_Element (struct Element* element_ptr, const struct Simulation* sim)
{
	if (element_ptr->type == POINT)
		return;

	switch (element_ptr->type) {
	case LINE: case TRI: case TET: case PYR:
		constructor_derived_FRSF_Solver_Element_std(element_ptr,sim);
		break;
	case QUAD: case HEX: case WEDGE:
	default:
		EXIT_UNSUPPORTED;
		break;
	}
}

void destructor_derived_FRSF_Solver_Element (struct Element* element_ptr)
{
	if (element_ptr->type == POINT)
		return;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void constructor_derived_FRSF_Solver_Element_std (struct Element* element_ptr, const struct Simulation* sim)
{
	struct FRSF_Solver_Element* frsf_s_e = (struct FRSF_Solver_Element*) element_ptr;
	struct const_Element* e = (struct const_Element*)element_ptr;
	struct Solver_Volume* s_vol = (struct Solver_Volume*) sim->volumes->first;

	const struct Multiarray_Operator* cv1_vs_vc_MO = constructor_operators("cv1","vsA","vcs","H_CF_P_PM0",e,sim);

	const int p = s_vol->p_ref;
	struct Multiarray_Operator cv1_vs_vc_ma1 = set_MO_from_MO(cv1_vs_vc_MO,1,(ptrdiff_t[]){0,0,p,p});

	const struct const_Matrix_d* cv0_vs_vc = get_operator__cv0_vs_vc(s_vol)->op_std;
	const struct const_Matrix_d* cv1_vs_vc = cv1_vs_vc_ma1.data[0]->op_std;
	const struct const_Matrix_d* cv0_vs_vc_inv = constructor_inverse_const_Matrix_d(cv0_vs_vc);

	frsf_s_e->D_frsf = constructor_mm_const_Matrix_d('N','N',1.0,cv0_vs_vc_inv,cv1_vs_vc,'R'); 
}
