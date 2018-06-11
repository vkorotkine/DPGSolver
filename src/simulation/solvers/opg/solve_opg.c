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
#include "definitions_test_case.h"

#include "face_solver.h"
#include "volume_solver_opg.h"

#include "multiarray.h"
#include "vector.h"

#include "compute_face_rlhs_opg.h"
#include "compute_rlhs.h"
#include "compute_volume_rlhs_opg.h"
#include "const_cast.h"
#include "simulation.h"
#include "solve.h"
#include "solve_implicit.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Fill \ref Solver_Storage_Implicit::b with the negated rhs values.
 *
 *  See comments in \ref solve_implicit.h for why the negated values are used here.
 */
static void fill_petsc_Vec_b_opg
	(const struct Simulation*const sim,       ///< Standard.
	 struct Solver_Storage_Implicit*const ssi ///< Standard.
	);

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "solve_opg_T.c"
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "solve_opg_T.c"
#include "undef_templates_type.h"

double compute_rlhs_opg (const struct Simulation*const sim, struct Solver_Storage_Implicit*const ssi)
{
	initialize_zero_memory_volumes(sim->volumes);
	compute_volume_rlhs_opg(sim,ssi,sim->volumes);
	compute_face_rlhs_opg(sim,ssi,sim->faces);
	compute_source_rhs_dg_like(sim);
	fill_petsc_Vec_b_opg(sim,ssi);

	struct Test_Case*const test_case = (struct Test_Case*)sim->test_case_rc->tc;
	switch (test_case->lhs_terms) {
	case LHS_FULL_NEWTON:
		break; // Do nothing
	case LHS_CFL_RAMPING:
		printf("*** Warning: CFL ramping not performed for the OPG method. ***\n");
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",test_case->lhs_terms);
		break;
	}

	return compute_max_rhs_dg_like(sim);
}

void set_petsc_Mat_row_col_opg
	(struct Solver_Storage_Implicit*const ssi, const struct OPG_Solver_Volume*const v_l, const int eq,
	 const struct OPG_Solver_Volume*const v_r, const int vr)
{
	const struct Solver_Volume*const sv_l = (struct Solver_Volume*) v_l,
	                          *const sv_r = (struct Solver_Volume*) v_r;
	ssi->row = (int)(sv_l->ind_dof_test+sv_l->test_s_coef->extents[0]*eq);
	ssi->col = (int)(sv_r->ind_dof_test+sv_r->test_s_coef->extents[0]*vr);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void fill_petsc_Vec_b_opg (const struct Simulation*const sim, struct Solver_Storage_Implicit*const ssi)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		const int ind_dof        = (int)((struct Solver_Volume*)curr)->ind_dof_test;
		const struct Multiarray_d*const rhs = ((struct Solver_Volume*)curr)->rhs;

		const int ni = (int)compute_size(rhs->order,rhs->extents);

		PetscInt    ix[ni];
		PetscScalar y[ni];

		for (int i = 0; i < ni; ++i) {
			ix[i] = ind_dof+i;
			y[i]  = -(rhs->data[i]);
		}

		VecSetValues(ssi->b,ni,ix,y,INSERT_VALUES);
	}
}
