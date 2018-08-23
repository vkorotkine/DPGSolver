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

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "macros.h"
#include "definitions_bc.h"
#include "definitions_opg.h"
#include "definitions_tol.h"

#include "def_templates_face_solver_opg.h"
#include "def_templates_volume_solver_opg.h"

#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"

#include "def_templates_compute_face_rlhs_opg.h"
#include "def_templates_compute_rlhs.h"
#include "def_templates_math_functions.h"
#include "def_templates_numerical_flux.h"
#include "def_templates_penalty_opg.h"
#include "def_templates_face_solver.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void constructor_rlhs_f_test_penalty_advection_upwind_T
	(const struct Flux_Ref_T*const flux_r, const struct Numerical_Flux_T*const num_flux,
	 struct Solver_Face_T*const s_face, struct Solver_Storage_Implicit*const ssi)
{
	UNUSED(flux_r); UNUSED(num_flux);

	const struct Face*const face = (struct Face*) s_face;
	const int bc = face->bc % BC_STEP_SC;
	switch (bc) {
	case BC_INFLOW: case BC_INFLOW_ALT1: case BC_INFLOW_ALT2:
	case BC_UPWIND: case BC_UPWIND_ALT1: case BC_UPWIND_ALT2:
	case BC_UPWIND_ALT3: case BC_UPWIND_ALT4: case BC_UPWIND_ALT5:
		return; // do not impose a boundary condition for the test functions.
		break;
	case BC_OUTFLOW: case BC_OUTFLOW_ALT1: case BC_OUTFLOW_ALT2:
		break; // continue below.
	default:
		EXIT_ERROR("Unsupported: %d\n",face->bc);
		break;
	}

	/** It is assumed that the initial solution for \ref Solver_Volume_T::test_s_coef is equal to zero on all
	 *  boundary faces which require the addition of the penalty term (outflow faces). It is also currently assumed
	 *  that the exact solution test functions are equal to zero (i.e. that the value of `g` described in
	 *  \ref constructor_rlhs_f_b_test_penalty_T is equal to zero).
	 */

	const struct Solver_Volume_T*const s_vol         = (struct Solver_Volume_T*) face->neigh_info[0].volume;
	const struct OPG_Solver_Face_T*const opg_s_face  = (struct OPG_Solver_Face_T*) s_face;

	const struct Lhs_Operators_OPG_T*const ops = constructor_Lhs_Operators_OPG_T(opg_s_face); // destructed

	const int ind_sc = (s_face->cub_type == 'c');
	const int p_t_p = get_set_degree_poly(NULL,"tp")[ind_sc];
	const double scale = pow(face->h,s_vol->p_ref+p_t_p);

	// Note: -ve sign included here. Would need to be moved for g != 0.
	const struct const_Matrix_T*const lhs_r =
		constructor_mm_diag_const_Matrix_R_T(-PENALTY_SCALING_OPG/scale,ops->cv0_vt_fc[0],ops->wJ_fc,'L',false); // d.

	const struct const_Matrix_T*const lhs =
		constructor_mm_RT_const_Matrix_T('T','N',1.0,ops->cv0_vt_fc[0],lhs_r,'R'); // destructed
	destructor_Lhs_Operators_OPG_T(s_face,ops);

#if TYPE_RC == TYPE_REAL
	if (ssi != NULL) {
		const int*const n_vr_eq = get_set_n_var_eq(NULL);
		const struct OPG_Solver_Volume_T*const opg_s_vol = (struct OPG_Solver_Volume_T*) face->neigh_info[0].volume;
		for (int vr = 0; vr < n_vr_eq[0]; ++vr) {
		for (int eq = 0; eq < n_vr_eq[1]; ++eq) {
			if (eq != vr)
				continue;
			set_petsc_Mat_row_col_opg(ssi,opg_s_vol,eq,opg_s_vol,vr);
			add_to_petsc_Mat(ssi,lhs);
		}}
	}
#elif TYPE_RC == TYPE_COMPLEX
	assert(ssi == NULL);
	mm_NNC_Multiarray_TTT(1.0,1.0,lhs,(struct const_Multiarray_T*)s_vol->test_s_coef,s_vol->rhs);
#endif
	destructor_const_Matrix_T(lhs_r);
	destructor_const_Matrix_T(lhs);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

#include "undef_templates_face_solver_opg.h"
#include "undef_templates_volume_solver_opg.h"

#include "undef_templates_matrix.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"

#include "undef_templates_compute_face_rlhs_opg.h"
#include "undef_templates_compute_rlhs.h"
#include "undef_templates_math_functions.h"
#include "undef_templates_numerical_flux.h"
#include "undef_templates_penalty_opg.h"
#include "undef_templates_face_solver.h"
