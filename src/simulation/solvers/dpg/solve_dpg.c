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

#include "solve_dpg.h"

#include <assert.h>

#include "macros.h"

#include "face_solver.h"
#include "volume_solver.h"

#include "multiarray.h"
#include "vector.h"

#include "compute_all_rlhs_dpg.h"
#include "intrusive.h"
#include "simulation.h"
#include "solve.h"
#include "solve_implicit.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void update_ind_dof_dpg (const struct Simulation* sim)
{
	ptrdiff_t dof = 0;
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Solver_Face* s_face = (struct Solver_Face*) curr;

		s_face->ind_dof = dof;

		struct Multiarray_d* nf_coef = s_face->nf_coef;
		dof += compute_size(nf_coef->order,nf_coef->extents);
	}

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume* s_vol = (struct Solver_Volume*) curr;

		s_vol->ind_dof = dof;

		struct Multiarray_d* sol_coef = s_vol->sol_coef;
		dof += compute_size(sol_coef->order,sol_coef->extents);
	}
	assert(dof == compute_dof(sim));
}

struct Vector_i* constructor_nnz_dpg (const struct Simulation* sim)
{
	assert(sim->test_case->has_2nd_order == false); // Add support.

	const ptrdiff_t dof = compute_dof(sim);
	struct Vector_i* nnz = constructor_zero_Vector_i(dof); // returned

	// Volume contribution (Diagonal)
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume* s_vol = (struct Solver_Volume*) curr;

		struct Multiarray_d* sol_coef = s_vol->sol_coef;
		const ptrdiff_t size = compute_size(sol_coef->order,sol_coef->extents);
		increment_nnz(nnz,s_vol->ind_dof,size,size);
	}

	// Face contributions (Diagonal and Off-diagonal)
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Face* face          = (struct Face*) curr;
		struct Solver_Face* s_face = (struct Solver_Face*) curr;

		// Diagonal
		struct Multiarray_d* nf_coef = s_face->nf_coef;
		const ptrdiff_t size_nf = compute_size(nf_coef->order,nf_coef->extents);
		increment_nnz(nnz,s_face->ind_dof,size_nf,size_nf);

		// Off-diagonal
		int ind_n = 0;
		struct Solver_Volume* s_vol = (struct Solver_Volume*) face->neigh_info[ind_n].volume;
		struct Multiarray_d* sol_coef = s_vol->sol_coef;
		ptrdiff_t size_sol = compute_size(sol_coef->order,sol_coef->extents);

		increment_nnz(nnz,s_face->ind_dof,size_nf,size_sol);
		increment_nnz(nnz,s_vol->ind_dof,size_sol,size_nf);

		if (face->boundary)
			continue;

		ind_n = 1;
		s_vol = (struct Solver_Volume*) face->neigh_info[ind_n].volume;
		sol_coef = s_vol->sol_coef;
		size_sol = compute_size(sol_coef->order,sol_coef->extents);

		increment_nnz(nnz,s_face->ind_dof,size_nf,size_sol);
		increment_nnz(nnz,s_vol->ind_dof,size_sol,size_nf);
	}
	return nnz;
}

double compute_rlhs_dpg (const struct Simulation* sim, struct Solver_Storage_Implicit* s_store_i)
{
//	zero_memory_volumes(sim);
//	zero_memory_faces(sim);
	compute_all_rlhs_dpg(sim,s_store_i);
//	compute_face_rlhs_dpg(sim,s_store_i);
//	compute_source_rhs_dpg(sim);

//	fill_petsc_Vec_b_dpg(sim,s_store_i);

EXIT_ADD_SUPPORT;
return 0.0;
	//return compute_max_rhs(sim);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
