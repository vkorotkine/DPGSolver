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

#include <stdbool.h>

#include "matrix.h"

#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void add_to_lhs_p_r
	(const double alpha, const struct const_Matrix_d*const dgc_dsc[DIM], struct Matrix_d*const lhs_p_r,
	 const bool boundary_face_term, const struct Simulation*const sim)
{
	const struct Test_Case*const test_case = (struct Test_Case*) sim->test_case_rc->tc;
	const int n_vr = test_case->n_var;
	if (!boundary_face_term) {
		for (int d_g = 0; d_g < DIM; ++d_g) {
			const struct const_Matrix_d*const lhs_p_l = dgc_dsc[d_g];
			for (int vr_g = 0; vr_g < n_vr; ++vr_g) {
			for (int vr_s = 0; vr_s < n_vr; ++vr_s) {
				if (vr_g != vr_s)
					continue;
				set_scaled_block_Matrix_d(alpha,lhs_p_r,(vr_g+n_vr*(d_g))*lhs_p_r->ext_0,vr_s*lhs_p_r->ext_1,
							        lhs_p_l,0,0,lhs_p_l->ext_0,lhs_p_l->ext_1,'a');
			}}
		}
	} else {
		for (int d_g = 0; d_g < DIM; ++d_g) {
			const struct const_Matrix_d*const lhs_p_l = dgc_dsc[d_g];
			set_scaled_block_Matrix_d(alpha,lhs_p_r,n_vr*d_g*lhs_p_r->ext_0,0,lhs_p_l,
			                          0,0,lhs_p_l->ext_0,lhs_p_l->ext_1,'a');
		}
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
