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

#include "intrusive.h"
#include "simulation.h"
#include "solve.h"
#include "solve_implicit.h"

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
EXIT_UNSUPPORTED;
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Face* face          = (struct Face*) curr;
		struct Solver_Face* s_face = (struct Solver_Face*) curr;

		// Diagonal
		struct Multiarray_d* nf_coef = s_face->nf_coef;
		const ptrdiff_t size_nf = compute_size(nf_coef->order,nf_coef->extents);
		increment_nnz(nnz,s_face->ind_dof,size_nf,size_nf);

		// Off-diagonal
		struct Solver_Volume* s_vol[2] = { (struct Solver_Volume*) face->neigh_info[0].volume, NULL, };
		struct Multiarray_d* sol_coef[2] = { s_vol[0]->sol_coef, NULL, };
		ptrdiff_t size = compute_size(sol_coef[0]->order,sol_coef[0]->extents);
UNUSED(size);
//		increment_nnz(nnz,s_vol[0]->ind_dof,size[0],size[1]);
EXIT_ADD_SUPPORT;

		if (!face->boundary) {
		}

	}
	return nnz;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
