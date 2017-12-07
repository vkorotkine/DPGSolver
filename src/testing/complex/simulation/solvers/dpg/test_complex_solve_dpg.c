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

#include "test_complex_solve_dpg.h"

#include "complex_multiarray.h"

#include "test_complex_face_solver_dpg.h"
#include "test_complex_volume_solver_dpg.h"

#include "test_complex_boundary.h"
#include "test_complex_boundary_advection.h"
#include "test_complex_boundary_euler.h"
#include "test_complex_compute_all_rhs_dpg.h"
#include "test_support_computational_elements.h"
#include "test_support_multiarray.h"
#include "test_complex_test_case.h"


#include <assert.h>

#include "macros.h"
#include "definitions_intrusive.h"
#include "definitions_test_integration.h"

#include "multiarray.h"

#include "boundary_advection.h"
#include "boundary_euler.h"
#include "intrusive.h"
#include "simulation.h"
#include "solve.h"
#include "solve_implicit.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Compute the complex rhs terms associated with the input volume for the dpg scheme.
static void compute_rhs_cmplx_step_dpg_volume
	(const struct DPG_Solver_Volume_c* dpg_s_vol, ///< The \ref DPG_Solver_Volume_c.
	 struct Solver_Storage_Implicit* ssi,           ///< \ref Solver_Storage_Implicit.
	 const struct Simulation* sim                   ///< \ref Simulation.
	);

/// \brief Compute the complex rhs terms associated with the input face for the dpg scheme.
static void compute_rhs_cmplx_step_dpg_face
	(const struct DPG_Solver_Face_c* dpg_s_face, ///< The \ref DPG_Solver_Face_c.
	 struct Solver_Storage_Implicit* ssi,          ///< \ref Solver_Storage_Implicit.
	 const struct Simulation* sim                  ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void perturb_solution_dpg (const struct Simulation* sim)
{
	struct Test_Case_c* test_case = (struct Test_Case_c*) sim->test_case_rc->tc;
	assert(test_case->has_2nd_order == false); // Add support.

EXIT_ADD_SUPPORT; // Ensure that the real and complex solutions are receiving the same perturbation.
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume_c* s_vol = (struct Solver_Volume_c*) curr;
		perturb_Multiarray_c(s_vol->sol_coef,1e-5);
	}

	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Solver_Face_c* s_face = (struct Solver_Face_c*) curr;
		perturb_Multiarray_c(s_face->nf_coef,1e-5);
	}
}

void compute_lhs_cmplx_step_dpg (const struct Simulation* sim, struct Solver_Storage_Implicit* ssi)
{
	assert(sim->test_case_rc->is_real == false);
	/** As the use of `complex` PETSc Vec containers would require using a different build where **all** containers
	 *  would be complex, it was decided to store the complex portion of the computed rhs term directly in the PETSc
	 *  Mat for this case. */

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		const struct Solver_Volume_c* s_vol   = (struct Solver_Volume_c*) curr;
		struct DPG_Solver_Volume_c* dpg_s_vol = (struct DPG_Solver_Volume_c*) curr;
		struct Multiarray_c* sol_coef_c = s_vol->sol_coef;
		const ptrdiff_t n_col_l = compute_size(sol_coef_c->order,sol_coef_c->extents);
		for (int col_l = 0; col_l < n_col_l; ++col_l) {
			ssi->col = (int)s_vol->ind_dof+col_l;

			sol_coef_c->data[col_l] += CX_STEP*I;
			compute_rhs_cmplx_step_dpg_volume(dpg_s_vol,ssi,sim);
			sol_coef_c->data[col_l] -= CX_STEP*I;
		}
	}

	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		const struct Face* face = (struct Face*) curr;
		if (face->boundary)
			continue;

		const struct Solver_Face_c* s_face   = (struct Solver_Face_c*) curr;
		struct DPG_Solver_Face_c* dpg_s_face = (struct DPG_Solver_Face_c*) curr;
		struct Multiarray_c* nf_coef_c = s_face->nf_coef;
		const ptrdiff_t n_col_l = compute_size(nf_coef_c->order,nf_coef_c->extents);
		for (int col_l = 0; col_l < n_col_l; ++col_l) {
			ssi->col = (int)s_face->ind_dof+col_l;

			nf_coef_c->data[col_l] += CX_STEP*I;
			compute_rhs_cmplx_step_dpg_face(dpg_s_face,ssi,sim);
			nf_coef_c->data[col_l] -= CX_STEP*I;
		}
	}
	petsc_mat_vec_assemble(ssi);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void compute_rhs_cmplx_step_dpg_volume
	(const struct DPG_Solver_Volume_c* dpg_s_vol, struct Solver_Storage_Implicit* ssi,
	 const struct Simulation* sim)
{
	struct Intrusive_List* volumes_local = constructor_empty_IL(IL_VOLUME_SOLVER_DPG,NULL); // destructed
	struct Intrusive_Link* curr = (struct Intrusive_Link*) dpg_s_vol;

	const size_t sizeof_base = sizeof(struct DPG_Solver_Volume_c);
	push_back_IL(volumes_local,constructor_copied_Intrusive_Link(curr,sizeof_base,sizeof_base));

	compute_all_rhs_dpg_c(sim,ssi,volumes_local);
	destructor_IL(volumes_local);
}

static void compute_rhs_cmplx_step_dpg_face
	(const struct DPG_Solver_Face_c* dpg_s_face, struct Solver_Storage_Implicit* ssi,
	 const struct Simulation* sim)
{
	const struct Face* face = (struct Face*) dpg_s_face;
	assert(!face->boundary);

	for (int i = 0; i < 2; ++i) {
		const struct DPG_Solver_Volume_c* dpg_s_vol = (struct DPG_Solver_Volume_c*) face->neigh_info[i].volume;
		compute_rhs_cmplx_step_dpg_volume(dpg_s_vol,ssi,sim);
	}
}
