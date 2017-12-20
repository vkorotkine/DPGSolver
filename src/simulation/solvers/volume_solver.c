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

#include "volume_solver.h"

#include "element_solver.h"

#include "matrix.h"
#include "multiarray.h"

#include "const_cast.h"
#include "geometry.h"
#include "simulation.h"

#include "test_complex_volume_solver.h"
#include "complex_multiarray.h"

// Templated functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "volume_solver_T.c"
#include "undef_templates_type.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void copy_members_r_to_c_Solver_Volume (struct Solver_Volume_c*const s_vol, const struct Solver_Volume*const s_vol_r)
{
	const_cast_ptrdiff(&s_vol->ind_dof,s_vol_r->ind_dof);
	const_cast_i(&s_vol->p_ref,s_vol_r->p_ref);
	const_cast_i(&s_vol->ml,s_vol_r->ml);
	const_constructor_move_const_Multiarray_d
		(&s_vol->geom_coef,constructor_copy_const_Multiarray_d(s_vol_r->geom_coef)); // destructed

	s_vol->sol_coef  = constructor_copy_Multiarray_c_Multiarray_d(s_vol_r->sol_coef);  // destructed
	s_vol->grad_coef = constructor_copy_Multiarray_c_Multiarray_d(s_vol_r->grad_coef); // destructed

	const_constructor_move_const_Multiarray_d(
		&s_vol->metrics_vm,constructor_copy_const_Multiarray_d(s_vol_r->metrics_vm)); // destructed
	const_constructor_move_const_Multiarray_d(
		&s_vol->metrics_vc,constructor_copy_const_Multiarray_d(s_vol_r->metrics_vc)); // destructed
	const_constructor_move_const_Multiarray_d(
		&s_vol->jacobian_det_vc,constructor_copy_const_Multiarray_d(s_vol_r->jacobian_det_vc)); // destructed
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
