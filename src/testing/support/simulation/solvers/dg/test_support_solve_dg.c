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
#include "test_support_solve.h"
#include "test_support_computational_elements.h"

#include <assert.h>

#include "macros.h"
#include "definitions_elements.h"
#include "definitions_intrusive.h"
#include "definitions_tol.h"
#include "definitions_test_integration.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "face_solver_dg.h"
#include "volume_solver_dg.h"

#include "computational_elements.h"
#include "compute_face_rlhs.h"
#include "compute_face_rlhs_dg.h"
#include "compute_grad_coef_dg.h"
#include "compute_volume_rlhs_dg.h"
#include "const_cast.h"
#include "intrusive.h"
#include "simulation.h"
#include "solve.h"
#include "solve_dg.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Compute the complex rhs terms based on the value of \ref CHECK_LIN.
static void compute_rhs_cmplx_step_dg
	(struct Intrusive_List* volumes_local, ///< The list of volumes over which to iterate.
	 struct Intrusive_List* faces_local,   ///< The list of faces over which to iterate.
	 const struct Simulation* sim          ///< Standard.
	);

/// \brief Set a column of the lhs matrix using the values of the complex rhs for the dg scheme.
static void set_col_lhs_cmplx_step_dg
	(const int col_l,                       ///< The local (to the volume solution dof) column index.
	 const struct Solver_Volume_c* s_vol_c, /**< The \ref Solver_Volume_T associated with the current column of the
	                                         *   matrix. */
	 struct Intrusive_List* volumes_local,  ///< The list of volumes over which to iterate.
	 struct Solver_Storage_Implicit* ssi    ///< \ref Solver_Storage_Implicit.
	);

/** \brief Constructor for the 'S'olution norm operator required for the generalized eigenvalue problem to compute the
 *         inf-sup constant.
 *  \return A petsc Mat holding the operator. */
static Mat constructor_Mat_S_dg
	(const struct Simulation*const sim ///< Standard.
	 );

/** \brief Constructor for the inverse of the 'T'est norm operator required for the generalized eigenvalue problem to
 *         compute the inf-sup constant.
 *  \return A petsc Mat holding the operator. */
static Mat constructor_Mat_T_inv_dg
	(const struct Simulation*const sim ///< Standard.
	 );

// Interface functions ********************************************************************************************** //

void compute_lhs_cmplx_step_dg (const struct Simulation* sim, struct Solver_Storage_Implicit* ssi)
{
	assert(!sim->test_case_rc->is_real);

	struct Test_Case_c* test_case = (struct Test_Case_c*) sim->test_case_rc->tc;
	assert(test_case->solver_method_curr == 'e'); // Should not be computing Jacobian terms.

	for (struct Intrusive_Link* curr_c = sim->volumes->first; curr_c; curr_c = curr_c->next) {
		struct Volume* vol = (struct Volume*) curr_c;
/// \todo Add a "volume_central" for volume/face terms, use volumes_local (including all neigh) for grad_coef.
		struct Intrusive_List* volumes_local = constructor_Volumes_local_neigh_only(vol,sim);
		struct Intrusive_List* faces_local   = constructor_Faces_local_neigh_only(vol,sim);

		struct Solver_Volume_c* s_vol = (struct Solver_Volume_c*) curr_c;
		struct Multiarray_c* sol_coef_c = s_vol->sol_coef;
		const ptrdiff_t n_col_l = compute_size(sol_coef_c->order,sol_coef_c->extents);
		for (int col_l = 0; col_l < n_col_l; ++col_l) {
			sol_coef_c->data[col_l] += CX_STEP*I;
			compute_rhs_cmplx_step_dg(volumes_local,faces_local,sim);
			sol_coef_c->data[col_l] -= CX_STEP*I;

			set_col_lhs_cmplx_step_dg(col_l,(struct Solver_Volume_c*)curr_c,volumes_local,ssi);
		}
		destructor_IL(volumes_local,true);
		destructor_IL(faces_local,true);
	}
	petsc_mat_vec_assemble(ssi);
}

const struct Gen_Eig_Data* constructor_Gen_Eig_Data_inf_sup_dg (struct Simulation*const sim)
{
	struct Test_Case* test_case = (struct Test_Case*) sim->test_case_rc->tc;
	test_case->solver_method_curr = 'i';

	constructor_derived_Elements(sim,IL_ELEMENT_SOLVER_DG);       // destructed
	constructor_derived_computational_elements(sim,IL_SOLVER_DG); // destructed
	initialize_zero_memory_volumes(sim->volumes);

	// Compute A
	struct Solver_Storage_Implicit*const ssi = constructor_Solver_Storage_Implicit(sim); // destructed
	compute_volume_rlhs_dg(sim,ssi,sim->volumes);
	compute_face_rlhs_dg(sim,ssi,sim->faces);
	Mat A = ssi->A;
	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

	/* Compute S and T_inv.
	 * The inverse of the test norm operator is computed directly as it only required local block inversion while
	 * the available Petsc functions to solve for inv(T)*A (i.e. using MatMatSolve to solve T*X = A) necessarily
	 * returns a dense matrix for X which is not necessary. */
	Mat S     = constructor_Mat_S_dg(sim);     // moved
	Mat T_inv = constructor_Mat_T_inv_dg(sim); // destroyed

	Mat T_inv_A;
	Mat A_T_inv_A;
	MatMatMult(T_inv,A,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&T_inv_A); // destroyed
	MatDestroy(&T_inv);

	MatTransposeMatMult(A,T_inv_A,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&A_T_inv_A); // moved
	MatDestroy(&T_inv_A);

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

/** \brief Constructor for a \ref Solver_Storage_Implicit container for the dg method's 'S'olution norm.
 *  \return See brief. */
struct Solver_Storage_Implicit* constructor_Solver_Storage_Implicit_dg_S
	(const struct Simulation*const sim ///< Standard.
	 );

static void compute_rhs_cmplx_step_dg
	(struct Intrusive_List* volumes_local, struct Intrusive_List* faces_local, const struct Simulation* sim)
{
	initialize_zero_memory_volumes_c(volumes_local);
	compute_grad_coef_dg_c(sim,volumes_local,faces_local);
	switch (CHECK_LIN) {
	case CHECK_LIN_VOLUME:
		compute_volume_rlhs_dg_c(sim,NULL,volumes_local);
		break;
	case CHECK_LIN_FACE:
		compute_face_rlhs_dg_c(sim,NULL,faces_local);
		break;
	case CHECK_LIN_ALL:
		compute_volume_rlhs_dg_c(sim,NULL,volumes_local);
		compute_face_rlhs_dg_c(sim,NULL,faces_local);
		break;
	default:
		EXIT_ERROR("Unsupported: %d.\n",CHECK_LIN);
		break;
	}
}

static void set_col_lhs_cmplx_step_dg
	(const int col_l, const struct Solver_Volume_c* s_vol_c, struct Intrusive_List* volumes_local,
	 struct Solver_Storage_Implicit* ssi)
{
	ssi->col = (int)s_vol_c->ind_dof+col_l;

	for (struct Intrusive_Link* curr_r = volumes_local->first; curr_r; curr_r = curr_r->next) {
		struct Solver_Volume_c* s_vol_r = (struct Solver_Volume_c*) curr_r;
		ssi->row = (int)s_vol_r->ind_dof+0;

		struct Multiarray_c* sol_coef_r = s_vol_r->sol_coef;

		const ptrdiff_t ext_0 = compute_size(sol_coef_r->order,sol_coef_r->extents);
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

static Mat constructor_Mat_S_dg (const struct Simulation*const sim)
{
	struct Solver_Storage_Implicit*const ssi = constructor_Solver_Storage_Implicit_dg_S(sim); // destructed/returned

	const int n_vr = get_set_n_var_eq(NULL)[0];
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		const struct Solver_Volume*const s_vol = (struct Solver_Volume*) curr;
		const struct const_Matrix_d*const m = constructor_mass(s_vol); // destructed
		const struct const_Matrix_d*const m_b = constructor_block_diagonal_const_Matrix_d(m,n_vr); // destructed
		destructor_const_Matrix_d(m);

		set_petsc_Mat_row_col_dg(ssi,s_vol,0,s_vol,0);
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

static Mat constructor_Mat_T_inv_dg (const struct Simulation*const sim)
{
	struct Solver_Storage_Implicit*const ssi = constructor_Solver_Storage_Implicit_dg_S(sim); // destructed/returned

	const int n_vr = get_set_n_var_eq(NULL)[0];
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		const struct Solver_Volume*const s_vol = (struct Solver_Volume*) curr;
		const struct const_Matrix_d*const m_inv = constructor_inverse_mass(s_vol,NULL); // destructed
		const struct const_Matrix_d*const m_inv_b = constructor_block_diagonal_const_Matrix_d(m_inv,n_vr); // d.
		destructor_const_Matrix_d(m_inv);

		set_petsc_Mat_row_col_dg(ssi,s_vol,0,s_vol,0);
		add_to_petsc_Mat(ssi,m_inv_b);
		destructor_const_Matrix_d(m_inv_b);
	}

	Mat T_inv = ssi->A;

	ssi->do_not_destruct_A = true;
	destructor_Solver_Storage_Implicit(ssi);

	MatAssemblyBegin(T_inv,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(T_inv,MAT_FINAL_ASSEMBLY);
	return T_inv;
}

// Level 1 ********************************************************************************************************** //

struct Solver_Storage_Implicit* constructor_Solver_Storage_Implicit_dg_S (const struct Simulation*const sim)
{
	update_ind_dof_d(sim);
	struct Vector_i*const nnz = constructor_nnz_dg(true,sim); // destructed
	const ptrdiff_t dof_solve = nnz->ext_0;

	struct Solver_Storage_Implicit*const ssi = calloc(1,sizeof *ssi); // free

	MatCreateSeqAIJ(MPI_COMM_WORLD,(PetscInt)dof_solve,(PetscInt)dof_solve,0,nnz->data,&ssi->A); // destructed
	MatSetFromOptions(ssi->A);
	MatSetUp(ssi->A);

	destructor_Vector_i(nnz);
	return ssi;
}
