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

#include "solve.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"
#include "definitions_intrusive.h"

#include "face_solver.h"
#include "volume_solver.h"

#include "multiarray.h"
#include "vector.h"

#include "geometry.h"
#include "intrusive.h"
#include "simulation.h"
#include "solve_explicit.h"
#include "solve_implicit.h"
#include "solution.h"
#include "test_case.h"

#include "solve_dg.h"
#include "solve_dpg.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "solve_T.c"

void solve_for_solution (struct Simulation* sim)
{
	assert(sim->volumes->name == IL_VOLUME_SOLVER);
	assert(sim->faces->name   == IL_FACE_SOLVER);

	set_up_solver_geometry(sim);
	set_initial_solution(sim);

	struct Test_Case* test_case = (struct Test_Case*)sim->test_case_rc->tc;
	switch (test_case->solver_proc) {
	case SOLVER_E:
		solve_explicit(sim);
		break;
	case SOLVER_I:
		solve_implicit(sim);
		break;
	case SOLVER_EI:
		solve_explicit(sim);
		solve_implicit(sim);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",test_case->solver_proc);
		break;
	}
}

double compute_rhs (const struct Simulation* sim)
{
/// \todo Add assertions relevant to rhs.
	double max_rhs = 0.0;

	switch (sim->method) {
	case METHOD_DG:
		max_rhs = compute_rhs_dg(sim);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",sim->method);
		break;
	}

	return max_rhs;
}

double compute_rlhs (const struct Simulation* sim, struct Solver_Storage_Implicit* s_store_i)
{
	double max_rhs = 0.0;

	switch (sim->method) {
		case METHOD_DG:  max_rhs = compute_rlhs_dg(sim,s_store_i);    break;
		case METHOD_DPG: max_rhs = compute_rlhs_dpg(sim,s_store_i);   break;
		default:         EXIT_ERROR("Unsupported: %d\n",sim->method); break;
	}

	return max_rhs;
}

void enforce_positivity_highorder (struct Solver_Volume* s_vol, const struct Simulation* sim)
{
	UNUSED(s_vol);
	UNUSED(sim);
	EXIT_ADD_SUPPORT;
}

void destructor_Solver_Storage_Implicit (struct Solver_Storage_Implicit* ssi)
{
	MatDestroy(&ssi->A);
	VecDestroy(&ssi->b);

	free(ssi);
}

void increment_nnz (struct Vector_i* nnz, const ptrdiff_t ind_dof, const ptrdiff_t n_row, const ptrdiff_t n_col)
{
	assert(ind_dof >= 0);

	const ptrdiff_t i_max = ind_dof+n_row;
	assert(i_max <= nnz->ext_0);

	for (ptrdiff_t i = ind_dof; i < i_max; ++i)
		nnz->data[i] += (int)n_col;
}

void petsc_mat_vec_assemble (struct Solver_Storage_Implicit* ssi)
{
	MatAssemblyBegin(ssi->A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ssi->A,MAT_FINAL_ASSEMBLY);
	VecAssemblyBegin(ssi->b);
	VecAssemblyEnd(ssi->b);
}

ptrdiff_t compute_dof_sol_1st (const struct Simulation* sim)
{
	ptrdiff_t dof = 0;
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume* s_vol = (struct Solver_Volume*) curr;
		dof += compute_size(s_vol->sol_coef->order,s_vol->sol_coef->extents);
	}
	return dof;
}

ptrdiff_t compute_dof_schur (const char dof_type, const struct Simulation* sim)
{
	ptrdiff_t dof = 0;
	switch (dof_type) {
	case 'f':
		dof += compute_dof_faces(sim);
		break;
	case 'v':
		dof += compute_dof_volumes(sim);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",dof_type);
		break;
	}
	return dof;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
