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

#include "solve_dg_T.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"
#include "definitions_test_case.h"


#include "def_templates_solve_dg.h"

#include "def_templates_compute_face_rlhs.h"
#include "def_templates_face_solver.h"
#include "def_templates_volume_solver.h"
#include "def_templates_volume_solver_dg.h"

#include "def_templates_multiarray.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

struct Vector_i* constructor_nnz_dg_T (const bool diag_only, const struct Simulation* sim)
{
	const ptrdiff_t dof = compute_dof(sim);
	struct Vector_i* nnz = constructor_zero_Vector_i(dof); // returned

	// Diagonal contribution
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) curr;

		struct Multiarray_T* sol_coef = s_vol->sol_coef;
		const ptrdiff_t size = compute_size(sol_coef->order,sol_coef->extents);
		increment_nnz(nnz,s_vol->ind_dof,size,size);
	}

	// Off-diagonal contributions
	if (!diag_only) {
		for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
			struct Face* face = (struct Face*) curr;
			if (face->boundary)
				continue;

			struct Solver_Volume_T* s_vol[2] = { (struct Solver_Volume_T*) face->neigh_info[0].volume,
			                                     (struct Solver_Volume_T*) face->neigh_info[1].volume, };

			struct Multiarray_T* sol_coef[2] = { s_vol[0]->sol_coef, s_vol[1]->sol_coef, };
			const ptrdiff_t size[2] = { compute_size(sol_coef[0]->order,sol_coef[0]->extents),
			                            compute_size(sol_coef[1]->order,sol_coef[1]->extents), };

			increment_nnz(nnz,s_vol[0]->ind_dof,size[0],size[1]);
			increment_nnz(nnz,s_vol[1]->ind_dof,size[1],size[0]);
		}
	}
	return nnz;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

#include "undef_templates_solve_dg.h"

#include "undef_templates_compute_face_rlhs.h"
#include "undef_templates_face_solver.h"
#include "undef_templates_volume_solver.h"
#include "undef_templates_volume_solver_dg.h"

#include "undef_templates_multiarray.h"
