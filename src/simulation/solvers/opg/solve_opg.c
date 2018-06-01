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

#include "solve_opg.h"

#include <assert.h>

#include "macros.h"

#include "face_solver.h"
#include "volume_solver_opg.h"

#include "multiarray.h"
#include "vector.h"

#include "compute_volume_rlhs_opg.h"
#include "const_cast.h"
#include "simulation.h"
#include "solve.h"
#include "solve_implicit.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "solve_opg_T.c"

double compute_rlhs_opg (const struct Simulation*const sim, struct Solver_Storage_Implicit*const ssi)
{
	compute_volume_rlhs_opg(sim,ssi,sim->volumes);
	EXIT_ADD_SUPPORT;

	UNUSED(sim);
	EXIT_UNSUPPORTED; // add_to_petsc_Mat_Vec_opg (should be similar to what is done for dpg, likely in compute_rlhs)

	return compute_max_rhs_from_ssi(ssi);
}

void set_petsc_Mat_row_col_opg
	(struct Solver_Storage_Implicit*const ssi, const struct OPG_Solver_Volume*const v_l, const int eq,
	 const struct OPG_Solver_Volume*const v_r, const int vr)
{
	const struct Solver_Volume*const sv_l = (struct Solver_Volume*) v_l,
	                          *const sv_r = (struct Solver_Volume*) v_r;
	ssi->row = (int)(sv_l->ind_dof_test+v_l->test_s_coef->extents[0]*eq);
	ssi->col = (int)(sv_r->ind_dof_test+v_r->test_s_coef->extents[0]*vr);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
