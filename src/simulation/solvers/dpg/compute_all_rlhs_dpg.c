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

#include "compute_all_rlhs_dpg.h"

#include "face_solver_dpg.h"
#include "volume_solver_dpg.h"
#include "element_solver_dpg.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "compute_rlhs.h"
#include "compute_face_rlhs.h"
#include "compute_volume_rlhs.h"
#include "flux.h"
#include "intrusive.h"
#include "multiarray_operator.h"
#include "numerical_flux.h"
#include "operator.h"
#include "simulation.h"
#include "solve.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Add entries from the current volume to Solver_Storage_Implicit::A and Solver_Storage_Implicit::b.
 *  \attention **When using the schur complement method to solve the global system, the block diagonal contributions are
 *             inverted before being added to global system matrix.**
 */
static void add_to_petsc_Mat_Vec_dpg
	(const struct Solver_Volume* s_vol,    ///< The current volume.
	 const struct const_Vector_d* rhs_neg, ///< The 'neg'ated local 'r'ight-'h'and 's'ide vector.
	 const struct const_Matrix_d* lhs,     ///< The local 'l'eft-'h'and 's'ide matrix.
	 struct Solver_Storage_Implicit* ssi,  ///< \ref Solver_Storage_Implicit.
	 const struct Simulation* sim          ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "compute_all_rlhs_dpg_T.c"

void compute_all_rlhs_dpg (const struct Simulation* sim, struct Solver_Storage_Implicit* ssi)
{
	assert(sim->volumes->name == IL_VOLUME_SOLVER_DPG);
	assert(sim->faces->name == IL_FACE_SOLVER_DPG);
	assert(sim->elements->name == IL_ELEMENT_SOLVER_DPG);

	/** \note Linearized terms are required for the computation of optimal test functions, even if only computing rhs
	 *        terms. */
	struct Test_Case* test_case = (struct Test_Case*)sim->test_case_rc->tc;
	assert(test_case->solver_method_curr == 'i');

	struct S_Params_DPG s_params = set_s_params_dpg(sim);
	struct Flux_Input* flux_i = constructor_Flux_Input(sim); // destructed

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
//struct Volume*        vol   = (struct Volume*) curr;
//printf("v_ind: %d\n",vol->index);
		struct Solver_Volume* s_vol         = (struct Solver_Volume*) curr;
		struct DPG_Solver_Volume* dpg_s_vol = (struct DPG_Solver_Volume*) curr;

		struct Flux_Ref* flux_r = constructor_Flux_Ref_vol(&s_params.spvs,flux_i,s_vol,sim); // destructed

		const struct const_Matrix_d* norm_op = s_params.constructor_norm_op(dpg_s_vol,flux_r,sim); // destructed

		s_params.compute_rlhs(norm_op,flux_r,dpg_s_vol,ssi,sim);
		destructor_const_Matrix_d(norm_op);
		destructor_Flux_Ref(flux_r);
	}
	destructor_Flux_Input(flux_i);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void add_to_petsc_Mat_Vec_dpg
	(const struct Solver_Volume* s_vol, const struct const_Vector_d* rhs_neg, const struct const_Matrix_d* lhs,
	 struct Solver_Storage_Implicit* ssi, const struct Simulation* sim)
{
	assert(sizeof(int) == sizeof(PetscInt));
	assert(sizeof(double) == sizeof(PetscScalar));

	const ptrdiff_t ext_0 = rhs_neg->ext_0;

	const struct const_Vector_i* idxm = constructor_petsc_idxm_dpg(ext_0,s_vol); // destructed.

	struct Test_Case* test_case = (struct Test_Case*)sim->test_case_rc->tc;
	if (test_case->use_schur_complement) {
		const ptrdiff_t dof_s = compute_size(s_vol->sol_coef->order,s_vol->sol_coef->extents),
		                dof_g = compute_size(s_vol->grad_coef->order,s_vol->grad_coef->extents);
		invert_sub_block_Matrix_d((struct Matrix_d*)lhs,0,0,dof_s);
		assert(dof_g == 0); // Add support.
	}

	MatSetValues(ssi->A,(PetscInt)ext_0,idxm->data,(PetscInt)ext_0,idxm->data,lhs->data,ADD_VALUES);
	VecSetValues(ssi->b,(PetscInt)ext_0,idxm->data,rhs_neg->data,ADD_VALUES);

	destructor_const_Vector_i(idxm);
}
