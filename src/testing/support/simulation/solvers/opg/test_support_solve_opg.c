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

#include "test_support_solve_opg.h"
#include "test_support_solve.h"
#include "test_support_multiarray.h"
#include "test_support_computational_elements.h"

#include <assert.h>

#include "macros.h"
#include "definitions_intrusive.h"
#include "definitions_tol.h"
#include "definitions_test_integration.h"

#include "multiarray.h"

#include "face_solver_opg.h"
#include "volume_solver_opg.h"

#include "compute_face_rlhs_opg.h"
#include "compute_volume_rlhs_opg.h"
#include "intrusive.h"
#include "simulation.h"
#include "solve.h"
#include "solve_opg.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Compute the complex rhs terms based on the value of \ref CHECK_LIN.
static void compute_rhs_cmplx_step_opg
	(struct Intrusive_List*const volumes_local_v, ///< The list of volumes over which to iterate for volume terms.
	 struct Intrusive_List*const volumes_local_f, ///< The list of volumes over which to iterate for face terms.
	 struct Intrusive_List*const faces_local,     ///< The list of faces over which to iterate.
	 const struct Simulation*const sim            ///< Standard.
		);

/// \brief Set a column of the lhs matrix using the values of the complex rhs for the opg scheme.
static void set_col_lhs_cmplx_step_opg
	(const int col_l,                            ///< The local (to the volume solution dof) column index.
	 const struct Solver_Volume_c*const s_vol_c, /**< The \ref Solver_Volume_T associated with the current column of
	                                              *   the matrix. */
	 struct Intrusive_List*const volumes_local,  ///< The list of volumes over which to iterate.
	 struct Solver_Storage_Implicit*const ssi    ///< \ref Solver_Storage_Implicit.
		);

// Interface functions ********************************************************************************************** //

void compute_lhs_cmplx_step_opg (const struct Simulation* sim, struct Solver_Storage_Implicit* ssi)
{
	assert(!sim->test_case_rc->is_real);

	struct Test_Case_c*const test_case = (struct Test_Case_c*) sim->test_case_rc->tc;
	assert(test_case->solver_method_curr == 'e'); // Should not be computing Jacobian terms.

	for (struct Intrusive_Link* curr_c = sim->volumes->first; curr_c; curr_c = curr_c->next) {
		struct Volume*const vol = (struct Volume*) curr_c;

		struct Intrusive_List*const volumes_local_v = constructor_Volumes_local_centre_only(vol);
		struct Intrusive_List*const volumes_local_f = constructor_Volumes_local_neigh_only(vol,sim);
		struct Intrusive_List*const faces_local     = constructor_Faces_local_neigh_only(vol,sim);

		struct Solver_Volume_c*const s_vol = (struct Solver_Volume_c*) curr_c;
		struct Multiarray_c* test_s_coef_c = s_vol->test_s_coef;
		const ptrdiff_t n_col_l = compute_size(test_s_coef_c->order,test_s_coef_c->extents);
		for (int col_l = 0; col_l < n_col_l; ++col_l) {
			test_s_coef_c->data[col_l] += CX_STEP*I;
			compute_rhs_cmplx_step_opg(volumes_local_v,volumes_local_f,faces_local,sim);
			test_s_coef_c->data[col_l] -= CX_STEP*I;

			set_col_lhs_cmplx_step_opg(col_l,(struct Solver_Volume_c*)curr_c,volumes_local_f,ssi);
		}
		destructor_IL(volumes_local_v,true);
		destructor_IL(volumes_local_f,true);
		destructor_IL(faces_local,true);
	}
	petsc_mat_vec_assemble(ssi);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void compute_rhs_cmplx_step_opg
	(struct Intrusive_List*const volumes_local_v, struct Intrusive_List*const volumes_local_f,
	 struct Intrusive_List*const faces_local, const struct Simulation*const sim)
{
	initialize_zero_memory_volumes_c(volumes_local_f);
	update_coef_s_v_opg_c(sim,volumes_local_v);
	update_coef_nf_f_opg_c(sim,faces_local);

	assert(get_set_has_1st_2nd_order(NULL)[1] == false); // Add support.
	switch (CHECK_LIN) {
	case CHECK_LIN_VOLUME:
		compute_volume_rlhs_opg_c(sim,NULL,volumes_local_v);
		break;
	case CHECK_LIN_FACE:
		compute_face_rlhs_opg_c(sim,NULL,faces_local);
		break;
	case CHECK_LIN_ALL:
		compute_volume_rlhs_opg_c(sim,NULL,volumes_local_v);
		compute_face_rlhs_opg_c(sim,NULL,faces_local);
		break;
	default:
		EXIT_ERROR("Unsupported: %d.\n",CHECK_LIN);
		break;
	}
}

static void set_col_lhs_cmplx_step_opg
	(const int col_l, const struct Solver_Volume_c*const s_vol_c, struct Intrusive_List*const volumes_local,
	 struct Solver_Storage_Implicit*const ssi)
{
	ssi->col = (int)s_vol_c->ind_dof_test+col_l;

	for (struct Intrusive_Link* curr_r = volumes_local->first; curr_r; curr_r = curr_r->next) {
		struct Solver_Volume_c* s_vol_r = (struct Solver_Volume_c*) curr_r;
		ssi->row = (int)s_vol_r->ind_dof_test+0;

		struct Multiarray_c* test_s_coef_r = s_vol_r->test_s_coef;

		const ptrdiff_t ext_0 = compute_size(test_s_coef_r->order,test_s_coef_r->extents);
		PetscInt idxm[ext_0];
		for (int i = 0; i < ext_0; ++i)
			idxm[i] = ssi->row+i;

		PetscInt idxn[1] = { ssi->col, };

		struct Multiarray_c* rhs_r = s_vol_r->rhs;
		PetscScalar vv[ext_0];
		for (int i = 0; i < ext_0; ++i)
			vv[i] = cimag(rhs_r->data[i])/CX_STEP;

		MatSetValues(ssi->A,(PetscInt)ext_0,idxm,1,idxn,vv,ADD_VALUES);
	}
}
