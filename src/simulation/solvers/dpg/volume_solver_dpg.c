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
/** \file
 */

#include "volume_solver_dpg.h"

#include "macros.h"

#include "element_solver_dpg.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

/** \brief Constructor for the H0 norm operator of the input volume.
 *  \return See brief. */
const struct const_Matrix_d* constructor_norm_op_H0
	(const struct DPG_Solver_Volume* dpg_s_vol ///< \ref DPG_Solver_Volume.
	);

// Interface functions ********************************************************************************************** //

void constructor_derived_DPG_Solver_Volume (struct Volume* volume_ptr, const struct Simulation* sim)
{
	struct DPG_Solver_Volume* dpg_s_vol = (struct DPG_Solver_Volume*) volume_ptr;
	UNUSED(sim);

	dpg_s_vol->norm_op_H0 = constructor_norm_op_H0(dpg_s_vol); // destructed
}

void destructor_derived_DPG_Solver_Volume (struct Volume* volume_ptr)
{
	struct DPG_Solver_Volume* dpg_s_vol = (struct DPG_Solver_Volume*) volume_ptr;

	destructor_const_Matrix_d(dpg_s_vol->norm_op_H0);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

const struct const_Matrix_d* constructor_norm_op_H0 (const struct DPG_Solver_Volume* dpg_s_vol)
{
	struct Volume* vol                       = (struct Volume*) dpg_s_vol;
	struct Solver_Volume* s_vol              = (struct Solver_Volume*) dpg_s_vol;
	const struct DPG_Solver_Element* dpg_s_e = (struct DPG_Solver_Element*) vol->element;

	const int p      = s_vol->p_ref,
	          curved = vol->curved;
	const struct Operator* cv0_vt_vc = get_Multiarray_Operator(dpg_s_e->cv0_vt_vc[curved],(ptrdiff_t[]){0,0,p,p});
	const struct const_Vector_d* w_vc = get_operator__w_vc__s_e(s_vol);

	const struct const_Vector_d jacobian_det_vc = interpret_const_Multiarray_as_Vector_d(s_vol->jacobian_det_vc);
	const struct const_Vector_d* wJ_vc = constructor_dot_mult_const_Vector_d(w_vc,&jacobian_det_vc,1); // destructed

	const struct const_Matrix_d* H0_l = cv0_vt_vc->op_std;
	const struct const_Matrix_d* H0_r = constructor_mm_diag_const_Matrix_d(1.0,H0_l,wJ_vc,'L',false); // destructed
	destructor_const_Vector_d(wJ_vc);

	const struct const_Matrix_d* H0 = constructor_mm_const_Matrix_d('T','N',1.0,H0_l,H0_r,'R'); // returned
	destructor_const_Matrix_d(H0_r);

	return H0;
}
