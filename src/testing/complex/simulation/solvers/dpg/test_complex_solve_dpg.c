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
#include "vector.h"

#include "face_solver_dpg.h"
#include "volume_solver_dpg.h"

#include "const_cast.h"
#include "intrusive.h"
#include "simulation.h"
#include "solve.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Constructor for the list of \ref Volume\*s including only the current volume.
 *  \return Standard. */
static struct Intrusive_List* constructor_Volumes_local_v
	(const struct Volume* vol ///< The \ref Volume.
	);

/** \brief Constructor for the list of \ref Volume\*s including only the volumes neighbouring the current face.
 *  \return Standard. */
static struct Intrusive_List* constructor_Volumes_local_f
	(const struct Face* face ///< The \ref Face.
	);

// Interface functions ********************************************************************************************** //

#include "def_templates_type_dc.h"
#include "solve_dpg_T.c"

void perturb_solution_dpg (const struct Simulation* sim)
{
/// \todo template this but as part of test_complex_solve_dpg_T such that functions are not accessible to the main code.
	if (sim->test_case_rc->is_real) {
		struct Test_Case* test_case = (struct Test_Case*) sim->test_case_rc->tc;
		assert(test_case->has_2nd_order == false); // Add support.

		for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
			struct Solver_Volume* s_vol = (struct Solver_Volume*) curr;
			perturb_Multiarray_d(s_vol->sol_coef,MAX_PERTURB);
		}

		for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
			struct Solver_Face* s_face = (struct Solver_Face*) curr;
			perturb_Multiarray_d(s_face->nf_coef,MAX_PERTURB);
		}
	} else {
		struct Test_Case_c* test_case = (struct Test_Case_c*) sim->test_case_rc->tc;
		assert(test_case->has_2nd_order == false); // Add support.

		for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
			struct Solver_Volume_c* s_vol = (struct Solver_Volume_c*) curr;
			perturb_Multiarray_c(s_vol->sol_coef,MAX_PERTURB);
		}

		for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
			struct Solver_Face_c* s_face = (struct Solver_Face_c*) curr;
			perturb_Multiarray_c(s_face->nf_coef,MAX_PERTURB);
		}
	}
}

void compute_lhs_cmplx_step_dpg (const struct Simulation* sim, struct Solver_Storage_Implicit* ssi)
{
	assert(sim->test_case_rc->is_real == false);
	/** As the use of `complex` PETSc Vec containers would require using a different build where **all** containers
	 *  would be complex, it was decided to store the complex portion of the computed rhs term directly in the PETSc
	 *  Mat for this case. */

#define TESTING
#ifdef TESTING
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Volume* vol = (struct Volume*) curr;
		struct Intrusive_List* volumes_local = constructor_Volumes_local_v(vol); // destructed
//volumes_local = sim->volumes;

		const struct Solver_Volume_c* s_vol = (struct Solver_Volume_c*) curr;
		struct Multiarray_c* sol_coef_c = s_vol->sol_coef;
		const ptrdiff_t n_col_l = compute_size(sol_coef_c->order,sol_coef_c->extents);
		for (int col_l = 0; col_l < n_col_l; ++col_l) {
			ssi->col = (int)s_vol->ind_dof+col_l;
printf(" vol: %d %d\n",col_l,ssi->col);
if (col_l < 1)
	continue;

//			sol_coef_c->data[col_l] += CX_STEP*I;
			sol_coef_c->data[col_l] += 1e-8;
			compute_all_rhs_dpg_c(sim,ssi,volumes_local);
EXIT_UNSUPPORTED;
			sol_coef_c->data[col_l] -= CX_STEP*I;
		}
		destructor_IL(volumes_local);
	}
#endif

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Volume* vol = (struct Volume*) curr;
		struct Intrusive_List* volumes_local = constructor_Volumes_local_v(vol); // destructed

		const struct Solver_Volume_c* s_vol = (struct Solver_Volume_c*) curr;
		struct Multiarray_c* sol_coef_c = s_vol->sol_coef;
		const ptrdiff_t n_col_l = compute_size(sol_coef_c->order,sol_coef_c->extents);
		for (int col_l = 0; col_l < n_col_l; ++col_l) {
			ssi->col = (int)s_vol->ind_dof+col_l;
printf(" vol: %d %d %d\n",col_l,ssi->col,vol->index);

			sol_coef_c->data[col_l] += CX_STEP*I;
			compute_all_rhs_dpg_c(sim,ssi,volumes_local);
			sol_coef_c->data[col_l] -= CX_STEP*I;
		}
		destructor_IL(volumes_local);
	}

	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		const struct Face* face = (struct Face*) curr;
		if (face->boundary)
			continue;

		struct Intrusive_List* volumes_local = constructor_Volumes_local_f(face); // destructed

		const struct Solver_Face_c* s_face = (struct Solver_Face_c*) curr;
		struct Multiarray_c* nf_coef_c = s_face->nf_coef;
		const ptrdiff_t n_col_l = compute_size(nf_coef_c->order,nf_coef_c->extents);
		for (int col_l = 0; col_l < n_col_l; ++col_l) {
			ssi->col = (int)s_face->ind_dof+col_l;
printf("face: %d %d\n",col_l,ssi->col);

			nf_coef_c->data[col_l] += CX_STEP*I;
			compute_all_rhs_dpg_c(sim,ssi,volumes_local);
			nf_coef_c->data[col_l] -= CX_STEP*I;
		}
		destructor_IL(volumes_local);
	}
	petsc_mat_vec_assemble(ssi);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct Intrusive_List* constructor_Volumes_local_v (const struct Volume* vol)
{
	struct Intrusive_List* volumes = constructor_empty_IL(IL_VOLUME_SOLVER_DPG,NULL); // returned

	const size_t sizeof_base = sizeof(struct DPG_Solver_Volume_c);
	struct Intrusive_Link* curr = (struct Intrusive_Link*) vol;

	// A copy is required such that the link in the global list is not modified.
	push_back_IL(volumes,constructor_copied_Intrusive_Link(curr,sizeof_base,sizeof_base));

	return volumes;
}

static struct Intrusive_List* constructor_Volumes_local_f (const struct Face* face)
{
	assert(!face->boundary);

	struct Intrusive_List* volumes = constructor_empty_IL(IL_VOLUME_SOLVER_DPG,NULL); // returned

	const size_t sizeof_base = sizeof(struct DPG_Solver_Volume_c);
	for (int i = 0; i < 2; ++i) {
		struct Intrusive_Link* curr = (struct Intrusive_Link*) face->neigh_info[i].volume;

		// A copy is required such that the link in the global list is not modified.
		push_back_IL(volumes,constructor_copied_Intrusive_Link(curr,sizeof_base,sizeof_base));
	}

	return volumes;
}
