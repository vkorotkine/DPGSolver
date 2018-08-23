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

#include "test_support_solve_dpg.h"
#include "test_support_solve.h"
#include "test_support_computational_elements.h"

#include <assert.h>

#include "macros.h"
#include "definitions_intrusive.h"
#include "definitions_tol.h"
#include "definitions_test_integration.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "face_solver_dpg.h"
#include "volume_solver_dpg.h"

#include "computational_elements.h"
#include "compute_all_rlhs_dpg.h"
#include "const_cast.h"
#include "intrusive.h"
#include "simulation.h"
#include "solve.h"
#include "solve_dpg.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Constructor for the list of \ref Volume\*s including only the volumes neighbouring the current face.
 *  \return Standard. */
static struct Intrusive_List* constructor_Volumes_local_f
	(const struct Face* face ///< The \ref Face.
	);

/** \brief Constructor for the 'S'olution norm operator required for the generalized eigenvalue problem to compute the
 *         inf-sup constant.
 *  \return A petsc Mat holding the operator. */
static Mat constructor_Mat_S_dpg
	(const struct Simulation*const sim ///< Standard.
	 );

// Interface functions ********************************************************************************************** //

void compute_lhs_cmplx_step_dpg (const struct Simulation* sim, struct Solver_Storage_Implicit* ssi)
{
	assert(sim->test_case_rc->is_real == false);
	/** As the use of `complex` PETSc Vec containers would require using a different build where **all** containers
	 *  would be complex, it was decided to store the complex portion of the computed rhs term directly in the PETSc
	 *  Mat for this case. */

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Volume* vol = (struct Volume*) curr;
		struct Intrusive_List* volumes_local = constructor_Volumes_local_centre_only(vol); // destructed

		const struct Solver_Volume_c* s_vol = (struct Solver_Volume_c*) curr;
		struct Multiarray_c* sol_coef_c = s_vol->sol_coef;
		const ptrdiff_t n_col_l = compute_size(sol_coef_c->order,sol_coef_c->extents);
		for (int col_l = 0; col_l < n_col_l; ++col_l) {
			ssi->col = (int)s_vol->ind_dof+col_l;
/* printf(" vol: %d %d %d\n",col_l,ssi->col,vol->index); */

			sol_coef_c->data[col_l] += CX_STEP*I;
			compute_all_rhs_dpg_c(sim,ssi,volumes_local);
			sol_coef_c->data[col_l] -= CX_STEP*I;
		}
		destructor_IL(volumes_local,true);
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
/* printf("face: %d %d\n",col_l,ssi->col); */

			nf_coef_c->data[col_l] += CX_STEP*I;
			compute_all_rhs_dpg_c(sim,ssi,volumes_local);
			nf_coef_c->data[col_l] -= CX_STEP*I;
		}
		destructor_IL(volumes_local,true);
	}

	if (test_case_explicitly_enforces_conservation(sim)) {
		for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
			struct Volume* vol = (struct Volume*) curr;
			struct Intrusive_List* volumes_local = constructor_Volumes_local_centre_only(vol); // destructed

			const struct Solver_Volume_c* s_vol = (struct Solver_Volume_c*) curr;
			struct Multiarray_c* l_mult_c = s_vol->l_mult;
			const ptrdiff_t n_col_l = compute_size(l_mult_c->order,l_mult_c->extents);
			for (int col_l = 0; col_l < n_col_l; ++col_l) {
				ssi->col = (int)s_vol->ind_dof_constraint+col_l;
/* printf(" vol (l_mult): %d %d %d\n",col_l,ssi->col,vol->index); */

				l_mult_c->data[col_l] += CX_STEP*I;
				compute_all_rhs_dpg_c(sim,ssi,volumes_local);
				l_mult_c->data[col_l] -= CX_STEP*I;
			}
			destructor_IL(volumes_local,true);
		}
	}
	petsc_mat_vec_assemble(ssi);
}

const struct Gen_Eig_Data* constructor_Gen_Eig_Data_inf_sup_dpg (struct Simulation*const sim)
{
	struct Test_Case* test_case = (struct Test_Case*) sim->test_case_rc->tc;
	test_case->solver_method_curr = 'i';

	constructor_derived_Elements(sim,IL_ELEMENT_SOLVER_DPG);       // destructed
	constructor_derived_computational_elements(sim,IL_SOLVER_DPG); // destructed
	initialize_zero_memory_volumes(sim->volumes);

	// Compute A' inv(T) A
	struct Solver_Storage_Implicit*const ssi = constructor_Solver_Storage_Implicit(sim); // destructed/moved
	compute_all_rlhs_dpg(sim,ssi,sim->volumes);
	Mat A_T_inv_A = ssi->A;
	MatAssemblyBegin(A_T_inv_A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A_T_inv_A,MAT_FINAL_ASSEMBLY);

	// Compute S.
	Mat S = constructor_Mat_S_dpg(sim); // moved

	ssi->do_not_destruct_A = true;
	destructor_Solver_Storage_Implicit(ssi);

	destructor_derived_Elements(sim,IL_ELEMENT_SOLVER);
	destructor_derived_computational_elements(sim,IL_SOLVER);

	struct Gen_Eig_Data*const ged = calloc(1,sizeof* ged); // returned
	ged->A = A_T_inv_A;
	ged->B = S;
	return ged;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Constructor for a \ref Solver_Storage_Implicit container for the dpg method's 'S'olution norm.
 *  \return See brief. */
struct Solver_Storage_Implicit* constructor_Solver_Storage_Implicit_dpg_S
	(const struct Simulation*const sim ///< Standard.
	 );

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

static Mat constructor_Mat_S_dpg (const struct Simulation*const sim)
{
	if (test_case_explicitly_enforces_conservation(sim))
		EXIT_ERROR("Add required contributions for the Lagrange multiplier terms.\n");
	if (get_set_has_1st_2nd_order(NULL)[1])
		EXIT_ERROR("Add support.\n");

	struct Solver_Storage_Implicit*const ssi = constructor_Solver_Storage_Implicit_dpg_S(sim); // destructed/returned

	const int eq = 0;
	const int vr = 0;
	const int n_vr = get_set_n_var_eq(NULL)[0];
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		const struct Face*const face = (struct Face*) curr;
		if (face->boundary)
			continue;

		const struct Solver_Face*const s_face = (struct Solver_Face*) curr;
		const struct const_Matrix_d*const m = constructor_mass_face_d(s_face); // destructed
		const struct const_Matrix_d*const m_b = constructor_block_diagonal_const_Matrix_d(m,n_vr); // destructed
		destructor_const_Matrix_d(m);

		ssi->row = (int)(s_face->ind_dof+s_face->nf_coef->extents[0]*eq);
		ssi->col = (int)(s_face->ind_dof+s_face->nf_coef->extents[0]*vr);
		add_to_petsc_Mat(ssi,m_b);
		destructor_const_Matrix_d(m_b);
	}

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		const struct Solver_Volume*const s_vol = (struct Solver_Volume*) curr;
		const struct const_Matrix_d*const m = constructor_mass(s_vol); // destructed
		const struct const_Matrix_d*const m_b = constructor_block_diagonal_const_Matrix_d(m,n_vr); // destructed
		destructor_const_Matrix_d(m);

		ssi->row = (int)(s_vol->ind_dof+s_vol->sol_coef->extents[0]*eq);
		ssi->col = (int)(s_vol->ind_dof+s_vol->sol_coef->extents[0]*vr);
		add_to_petsc_Mat(ssi,m_b);
		destructor_const_Matrix_d(m_b);
	}

	Mat S = ssi->A;

	ssi->do_not_destruct_A = true;
	destructor_Solver_Storage_Implicit(ssi);

	MatAssemblyBegin(S,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(S,MAT_FINAL_ASSEMBLY);
	return S;
}

// Level 1 ********************************************************************************************************** //

struct Solver_Storage_Implicit* constructor_Solver_Storage_Implicit_dpg_S (const struct Simulation*const sim)
{
	update_ind_dof_d(sim);
	struct Vector_i*const nnz = constructor_nnz_dpg(true,sim); // destructed
	const ptrdiff_t dof_solve = nnz->ext_0;

	struct Solver_Storage_Implicit*const ssi = calloc(1,sizeof *ssi); // free

	MatCreateSeqAIJ(MPI_COMM_WORLD,(PetscInt)dof_solve,(PetscInt)dof_solve,0,nnz->data,&ssi->A); // destructed
	MatSetFromOptions(ssi->A);
	MatSetUp(ssi->A);

	destructor_Vector_i(nnz);
	return ssi;
}
