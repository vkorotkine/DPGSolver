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

/** \brief Constructor for the Jacobian of the lhs volume term with respect to the solution coefficients of 1st order
 *         equations only.
 *  \return See additional comments.
 *
 *  The term returned from this function, henceforth denoted by dLds, will eventually be transposed and multiplied on
 *  the right by the solution coefficient related rows of the RHS vector. To perform the multiplication using a single
 *  BLAS 2 call (i.e. in the most efficient manner), the dLds term is stored in memory with the following constraints:
 *  - ext_0 = n_dof_t*n_eq (size of RHS(s_coef), where n_dof_t == number of test basis functions);
 *  - ext_1 = (n_dof_s*n_vr)^2 (number of columns of standard LHS matrix times number of components to linearize with
 *                              respect to);
 *  - layout = 'C' (Required to ensure that memory of each slice of the 3-tensor is contiguous).
 *
 *  This function is nearly identical to \ref constructor_lhs_v_1_T, except that the \ref Solver_Element::cv0_vs_vc
 *  operator is replaced by the \ref DPG_Solver_Element::cvcv0_vs_vc operator.
 */
static struct Matrix_d* constructor_dlhs_ds_v_1
	(const struct Flux_Ref* flux_r,             ///< \ref Flux_Ref_T.
	 const struct DPG_Solver_Volume* dpg_s_vol, ///< \ref DPG_Solver_Volume_T.
	 const struct Simulation* sim               ///< \ref Simulation.
	);

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

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Get the pointer to the appropriate \ref DPG_Solver_Element::cvcv0_vs_vc operator.
 *  \return See brief. */
const struct Operator* get_operator__cvcv0_vs_vc
	(const struct DPG_Solver_Volume* dpg_s_vol ///< The current volume.
	);

static struct Matrix_d* constructor_dlhs_ds_v_1
	(const struct Flux_Ref* flux_r, const struct DPG_Solver_Volume* dpg_s_vol, const struct Simulation* sim)
{
	const struct Solver_Volume* s_vol = (struct Solver_Volume*) dpg_s_vol;

	const struct Multiarray_Operator tw1_vt_vc = get_operator__tw1_vt_vc(s_vol);
	const struct Operator* cvcv0_vs_vc         = get_operator__cvcv0_vs_vc(dpg_s_vol);

	struct Test_Case* test_case = (struct Test_Case*)sim->test_case_rc->tc;
	const ptrdiff_t n_dof_t  = tw1_vt_vc.data[0]->op_std->ext_0,
	                n_dof_s  = s_vol->sol_coef->extents[0],
	                n_dof_s2 = n_dof_s*n_dof_s,
	                n_eq     = test_case->n_eq,
	                n_vr     = test_case->n_var,
	                n_vr2    = n_vr*n_vr;

	const ptrdiff_t ext_1 = tw1_vt_vc.data[0]->op_std->ext_1;
	struct Matrix_d* tw1_r  = constructor_empty_Matrix_d('R',n_dof_t,ext_1);               // destructed
	struct Matrix_d* lhs_ll = constructor_empty_Matrix_d('C',n_dof_t,n_dof_s);             // destructed
	struct Matrix_d* lhs_l  = constructor_empty_Matrix_d('C',n_dof_t,n_dof_s2);            // destructed
	struct Matrix_d* lhs    = constructor_empty_Matrix_d('C',n_eq*n_dof_t,n_vr2*n_dof_s2); // returned
set_to_value_Matrix_d(lhs,0.0);

	const struct const_Multiarray_d* d2fr_ds2_Ma = flux_r->d2fr_ds2;
	struct Vector_d d2fr_ds2 = { .ext_0 = d2fr_ds2_Ma->extents[0], .owns_data = false, .data = NULL, };

	for (int vr2 = 0; vr2 < n_vr; ++vr2) {
	for (int vr = 0; vr < n_vr; ++vr) {
	for (int eq = 0; eq < n_eq; ++eq) {
printf("eq, vr, vr2: %d %d %d\n",eq,vr,vr2);
		set_to_value_Matrix_d(tw1_r,0.0);
		for (int dim = 0; dim < DIM; ++dim) {
			const ptrdiff_t ind = compute_index_sub_container(
				d2fr_ds2_Ma->order,1,d2fr_ds2_Ma->extents,(ptrdiff_t[]){eq,vr,vr2,dim});
			d2fr_ds2.data = (double*)&d2fr_ds2_Ma->data[ind];
printf("%d\n",dim);
print_Vector_d(&d2fr_ds2);
			mm_diag_d('R',1.0,1.0,tw1_vt_vc.data[dim]->op_std,(struct const_Vector_d*)&d2fr_ds2,tw1_r,false);
		}

		mm_d('N','N',1.0,0.0,(struct const_Matrix_d*)tw1_r,cvcv0_vs_vc->op_std,lhs_l);
//print_Matrix_d(lhs_l);

		for (int dof_s = 0; dof_s < n_dof_s; ++dof_s) {
			const ptrdiff_t row = eq*n_dof_t,
			                col = (n_dof_s*n_vr*vr2+dof_s*n_vr+vr)*n_dof_s;
			lhs_ll->data = get_col_Matrix_d(dof_s*n_dof_s,lhs_l);
			set_block_Matrix_d(lhs,(struct const_Matrix_d*)lhs_ll,row,col,'i');
		}
if (eq == 3 && vr == 1 && vr2 == 1) {
print_Matrix_d(lhs);
EXIT_UNSUPPORTED;
}
	}}}
	destructor_Matrix_d(tw1_r);
	destructor_Matrix_d(lhs_l);
	destructor_Matrix_d(lhs_ll);
EXIT_UNSUPPORTED;

	return lhs;
}

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

// Level 1 ********************************************************************************************************** //

const struct Operator* get_operator__cvcv0_vs_vc (const struct DPG_Solver_Volume* dpg_s_vol)
{
	const struct Volume* vol           = (struct Volume*) dpg_s_vol;
	const struct Solver_Volume* s_vol  = (struct Solver_Volume*) dpg_s_vol;
	const struct DPG_Solver_Element* e = (struct DPG_Solver_Element*) vol->element;

	const int p = s_vol->p_ref,
	          curved = vol->curved;
	return get_Multiarray_Operator(e->cvcv0_vs_vc[curved],(ptrdiff_t[]){0,0,p,p});
}
