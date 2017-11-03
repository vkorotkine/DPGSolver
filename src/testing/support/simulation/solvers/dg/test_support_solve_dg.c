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

#include "test_support_solve_dg.h"

#include <assert.h>

#include "macros.h"
#include "definitions_elements.h"
#include "definitions_intrusive.h"
#include "definitions_test_integration.h"

#include "test_support_compute_volume_rhs_dg.h"
#include "test_support_complex_multiarray.h"
#include "test_support_computational_elements.h"
#include "test_support_multiarray.h"

#include "face.h"
#include "volume.h"
#include "volume_solver.h"
#include "volume_solver_dg_complex.h"

#include "complex_multiarray.h"
#include "multiarray.h"

#include "intrusive.h"
#include "simulation.h"
#include "solve.h"
#include "solve_implicit.h"

// Static function declarations ************************************************************************************* //

/** \brief Constructor for the list of \ref Volume\*s adjacent to (and including) the current volume.
 *  \return Standard. */
struct Intrusive_List* constructor_Volumes_local
	(const struct Volume* vol,    ///< The centre \ref Volume.
	 const struct Simulation* sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void perturb_solution_dg (const struct Simulation* sim)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next)
		perturb_Multiarray_d(((struct Solver_Volume*)curr)->sol_coef,1e-5);
}

void set_initial_solution_complex_dg (const struct Simulation* sim)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume* s_vol                 = (struct Solver_Volume*) curr;
		struct Complex_DG_Solver_Volume* c_dg_s_vol = (struct Complex_DG_Solver_Volume*) curr;

		set_Multiarray_c_Multiarray_d(c_dg_s_vol->sol_coef,(struct const_Multiarray_d*)s_vol->sol_coef);
	}
}

void compute_lhs_cmplx_step_dg (const struct Simulation* sim, struct Solver_Storage_Implicit* s_store_i)
{
	for (struct Intrusive_Link* curr_c = sim->volumes->first; curr_c; curr_c = curr_c->next) {
		struct Volume* vol = (struct Volume*) curr_c;
		struct Intrusive_List* volumes_local = constructor_Volumes_local(vol,sim);
//		struct Intrusive_List* faces_local   = constructor_Faces_local(vol,sim);
printf("v ind: %d\n",vol->index);

		struct Complex_DG_Solver_Volume* c_dg_s_vol_c = (struct Complex_DG_Solver_Volume*) curr_c;
		struct Multiarray_c* sol_coef_c = c_dg_s_vol_c->sol_coef;
		const ptrdiff_t n_col_l = compute_size(sol_coef_c->order,sol_coef_c->extents);
		for (int col_l = 0; col_l < n_col_l; ++col_l) {
			sol_coef_c->data[col_l] += CX_STEP*I;
			switch (CHECK_LIN) {
			case CHECK_LIN_VOLUME:
				compute_volume_rhs_dg_c(sim,volumes_local);
				break;
			case CHECK_LIN_FACE:
EXIT_ADD_SUPPORT;
				break;
			case CHECK_LIN_ALL:
EXIT_ADD_SUPPORT;
				break;
			default:
				EXIT_ERROR("Unsupported: %d.\n",CHECK_LIN);
				break;
			}

// External function here.
			struct Solver_Volume* s_vol_c = (struct Solver_Volume*) curr_c;
			s_store_i->col = s_vol_c->ind_dof+col_l;

			for (struct Intrusive_Link* curr_r = volumes_local->first; curr_r; curr_r = curr_r->next) {
// Diagonal only.
assert(CHECK_LIN == CHECK_LIN_VOLUME);
				struct Solver_Volume* s_vol_r = (struct Solver_Volume*) curr_r;
				if (s_vol_r->ind_dof != s_vol_c->ind_dof)
					continue;

				s_store_i->row = s_vol_r->ind_dof+0;

				struct Complex_DG_Solver_Volume* c_dg_s_vol_r = (struct Complex_DG_Solver_Volume*) curr_r;
				struct Multiarray_c* sol_coef_r = c_dg_s_vol_r->sol_coef;
				const ptrdiff_t ext_0 = compute_size(sol_coef_r->order,sol_coef_r->extents),
				                ext_1 = 1;

				PetscInt idxm[ext_0],
				         idxn[ext_1];

				for (int i = 0; i < ext_0; ++i)
					idxm[i] = s_store_i->row+i;

				for (int i = 0; i < ext_1; ++i)
					idxn[i] = s_store_i->col+i;

				struct Multiarray_c* rhs_r = c_dg_s_vol_r->rhs;
				PetscScalar vv[ext_0];
				for (int i = 0; i < ext_0; ++i)
					vv[i] = cimag(rhs_r->data[i])/CX_STEP;
//printf("%td %td %d %d %d %d %e %e %e\n",ext_0,ext_1,idxm[0],idxm[1],idxm[2],idxn[0],vv[0],vv[1],vv[2]);

				MatSetValues(s_store_i->A,ext_0,idxm,ext_1,idxn,vv,ADD_VALUES);
			}
			sol_coef_c->data[col_l] -= CX_STEP*I;
		}

		destructor_IL(volumes_local);
UNUSED(s_store_i);
	}
	petsc_mat_vec_assemble(s_store_i);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Check whether the \ref Volume is a neighbour (inclusive) of the current \ref Volume.
 *  \return `true` if is a neighbour (or is the current volume); `false` otherwise. */
static bool is_volume_neighbour
	(const struct Volume* vol,     ///< The volume under investigation.
	 const struct Volume* vol_curr ///< The current volume.
	);

struct Intrusive_List* constructor_Volumes_local (const struct Volume* vol, const struct Simulation* sim)
{
	struct Intrusive_List* volumes = constructor_empty_IL(IL_VOLUME_SOLVER_DG_COMPLEX,NULL);

	const size_t sizeof_base = sizeof(struct Complex_DG_Solver_Volume);
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		if (!is_volume_neighbour((struct Volume*)curr,vol))
			continue;

		push_back_IL(volumes,constructor_copied_Intrusive_Link(curr,sizeof_base,sizeof_base));
	}

	return volumes;
}

// Level 1 ********************************************************************************************************** //

static bool is_volume_neighbour (const struct Volume* vol, const struct Volume* vol_curr)
{
	bool is_neighbour = false;

	for (int f = 0; f < NFMAX; ++f) {
	for (int sf = 0; sf < NSUBFMAX; ++sf) {
		const struct Face* face = vol_curr->faces[f][sf];

		if (!face)
			continue;

		const int n_side = ( face->boundary ? 1 : 2 );
		for (int i = 0; i < n_side; ++i) {
			if (vol == face->neigh_info[i].volume)
				is_neighbour = true;
		}

		if (is_neighbour)
			break;
	}}
	return is_neighbour;
}
