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

#include "compute_rlhs.h"

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>

#include "matrix.h"
#include "multiarray.h"

#include "volume_solver.h"

#include "computational_elements.h"
#include "intrusive.h"
#include "math_functions.h"
#include "test_case.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "compute_rlhs_T.c"
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "compute_rlhs_T.c"
#include "undef_templates_type.h"

void compute_source_rhs_dg_like (const struct Simulation*const sim)
{
	assert(list_is_derived_from("solver",'v',sim));

	struct Test_Case*const test_case = (struct Test_Case*)sim->test_case_rc->tc;

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		const struct Solver_Volume*const s_vol = (struct Solver_Volume*) curr;
		test_case->compute_source_rhs(sim,s_vol,s_vol->rhs);
	}
}

double compute_max_rhs_dg_like (const struct Simulation*const sim)
{
	double max_rhs = 0.0;
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume*const s_vol = (struct Solver_Volume*) curr;

		struct Multiarray_d* rhs = s_vol->rhs;
		double max_rhs_curr = norm_d(compute_size(rhs->order,rhs->extents),rhs->data,"Inf");
		if (max_rhs_curr > max_rhs)
			max_rhs = max_rhs_curr;
	}
	return max_rhs;
}

void add_to_lhs_p_r
	(const double alpha, const struct const_Matrix_d*const dgc_dsc[DIM], struct Matrix_d*const lhs_p_r,
	 const bool boundary_face_term)
{
	const int n_vr = get_set_n_var_eq(NULL)[0];
	if (!boundary_face_term) {
		for (int d_g = 0; d_g < DIM; ++d_g) {
			const struct const_Matrix_d*const lhs_p_l = dgc_dsc[d_g];
			for (int vr_g = 0; vr_g < n_vr; ++vr_g) {
			for (int vr_s = 0; vr_s < n_vr; ++vr_s) {
				if (vr_g != vr_s)
					continue;
				set_scaled_block_Matrix_d(alpha,lhs_p_r,(vr_g+n_vr*(d_g))*lhs_p_l->ext_0,vr_s*lhs_p_l->ext_1,
							        lhs_p_l,0,0,lhs_p_l->ext_0,lhs_p_l->ext_1,'a');
			}}
		}
	} else {
		for (int d_g = 0; d_g < DIM; ++d_g) {
			const struct const_Matrix_d*const lhs_p_l = dgc_dsc[d_g];
			set_scaled_block_Matrix_d(alpha,lhs_p_r,d_g*lhs_p_l->ext_0,0,lhs_p_l,
			                          0,0,lhs_p_l->ext_0,lhs_p_l->ext_1,'a');
		}
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
