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

#include "macros.h"
#include "definitions_tol.h"

#include "face_solver_dpg.h"
#include "volume_solver_dpg.h"
#include "element_solver_dpg.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "computational_elements.h"
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


#include "test_complex_volume_solver_dpg.h"

#include "complex_multiarray.h"
#include "complex_matrix.h"

#include "test_complex_computational_elements.h"
#include "test_complex_compute_all_rhs_dpg.h"
#include "test_complex_test_case.h"

// Static function declarations ************************************************************************************* //

struct Norm_DPG;

/** \brief Add the contribution of the linearization of the optimal test functions wrt to the solution coefficients to
 *         the lhs term for the DPG scheme.
 *
 *  When the pde is not linear, the "form" is no longer linear in the test functions and the associated linearization
 *  contribution must also be included.
 */
static void add_to_lhs_opt__d_opt_t_ds
	(const struct Flux_Ref*const flux_r,             ///< \ref Flux_Ref_T.
	 const struct DPG_Solver_Volume*const dpg_s_vol, ///< \ref DPG_Solver_Volume_T.
	 const struct Norm_DPG*const norm,               ///< \ref Norm_DPG.
	 const struct const_Matrix_d*const opt_t,        ///< Computed optimal test functions.
	 const struct const_Vector_d*const rhs_std,      ///< The standard (DG-like) contribution to the rhs.
	 const struct const_Matrix_d*const lhs_std,      ///< The standard (DG-like) contribution to the lhs.
	 struct Matrix_d*const lhs_opt,                  ///< Container used to store the optimal lhs contribution.
	 const struct Simulation*const sim,              ///< \ref Simulation.
	 struct Simulation*const sim_c                   ///< The complex \ref Simulation.
	);

/** \brief Add the contribution from the lagrange multiplier enforcing conservation (if applicable).
 *
 *  \note When the numerical flux used to enforce the boundary conditions is non-linear, an additional contribution to
 *        the standard Lagrange multiplier terms (eq. (4.47) \cite Zienkiewicz2013_ch4) is added to the
 *        sol_coef-sol_coef portion of lhs_opt. */
static void add_to_lhs_opt__l_mult
	(struct Matrix_d**const lhs_opt_ptr,             ///< Pointer to the optimal lhs matrix contribution.
	 const struct const_Matrix_d*const lhs_std,      ///< The standard (DG-like) contribution to the lhs.
	 const struct DPG_Solver_Volume*const dpg_s_vol, ///< \ref DPG_Solver_Volume_T.
	 const struct Simulation*const sim               ///< \ref Simulation.
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

/** \brief Portion of \ref constructor_norm_DPG__h1_upwind constructing the Jacobian of the norm wrt the solution
 *         coefficients.
 *  \return See brief. */
static const struct const_Matrix_d* constructor_norm_DPG_dN_ds__h1_upwind
	(const struct DPG_Solver_Volume* dpg_s_vol, ///< See brief.
	 const struct Flux_Ref* flux_r,             ///< See brief.
	 const struct const_Matrix_d* n1_lt,        /**< The transpose of the left component of the 1st order term of the
	                                             *   h1 upwind norm. */
	 const struct Simulation* sim               ///< See brief.
	);

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "compute_all_rlhs_dpg_T.c"

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Get the pointer to the appropriate \ref DPG_Solver_Element::cvcv0_vs_vc operator.
 *  \return See brief. */
static const struct Operator* get_operator__cvcv0_vs_vc
	(const struct DPG_Solver_Volume* dpg_s_vol ///< The current volume.
	);

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

/** \brief Adds the linearization of the face boundary terms to the Hessian of the residual (rhs) wrt the solution
 *         coefficients. */
static void add_to_dlhs_ds__face_boundary
	(struct Matrix_d*const dlhs_ds,             ///< Defined for \ref add_to_dlhs_ds__norm.
	 const struct DPG_Solver_Volume* dpg_s_vol, ///< Defined for \ref add_to_rlhs__face_boundary.
	 const struct Simulation*const sim,         ///< Defined for \ref add_to_rlhs__face_boundary.
	 struct Simulation*const sim_c,             ///< Defined for \ref add_to_rlhs__face_boundary.
	 const char lin_method                      ///< Linearization method. Options: 'a'nalytical, 'c'omplex step.
	);

/** \brief Adds the linearization of \ref Norm_DPG::N to the Hessian of the residual (rhs) wrt the solution coefficients
 *         if the norm has a dependence on the solution.
 *
 *  Columns corresponding to the face dof are added to `dlhs_ds` here if there is a contribution to the linearization
 *  from the DPG norm.
 */
static void add_to_dlhs_ds__norm
	(struct Matrix_d**const dlhs_ds_ptr, /**< Pointer to the term storing the linearization of the LHS wrt to the
	                                      *   solution coefficients. */
	 const struct Norm_DPG* norm,        ///< \ref Norm_DPG.
	 const struct const_Matrix_d* opt_t  ///< The optimal test functions.
	);

/** \brief Get the appropriate sub-range of the \ref DPG_Solver_Element::cvcv1_vt_vc operators.
 *  \return See brief. */
static struct Multiarray_Operator get_operator__cvcv1_vt_vc__rlhs
	(const struct DPG_Solver_Volume* dpg_s_vol ///< The current volume.
	);

static void add_to_lhs_opt__d_opt_t_ds
	(const struct Flux_Ref*const flux_r, const struct DPG_Solver_Volume*const dpg_s_vol,
	 const struct Norm_DPG*const norm, const struct const_Matrix_d*const opt_t,
	 const struct const_Vector_d*const rhs_std, const struct const_Matrix_d*const lhs_std,
	 struct Matrix_d*const lhs_opt, const struct Simulation*const sim, struct Simulation*const sim_c)
{
	struct Test_Case* test_case = (struct Test_Case*) sim->test_case_rc->tc;
	if (test_case->is_linear)
		return;

	struct Matrix_d* dlhs_ds = constructor_dlhs_ds_v_1(flux_r,dpg_s_vol,sim); // destructed
	add_to_dlhs_ds__face_boundary(dlhs_ds,dpg_s_vol,sim,sim_c,'c');
	add_to_dlhs_ds__norm(&dlhs_ds,norm,opt_t);

	const struct const_Matrix_d* dopt_t_ds =
		constructor_sysv_const_Matrix_d(norm->N,(struct const_Matrix_d*)dlhs_ds); // destructed
	destructor_Matrix_d(dlhs_ds);

	const struct const_Vector_d* lhs_opt_t = constructor_mv_const_Vector_d('T',1.0,dopt_t_ds,rhs_std); // destructed
	destructor_const_Matrix_d(dopt_t_ds);

	const struct Solver_Volume*const s_vol = (struct Solver_Volume*) dpg_s_vol;

	// Note: If the `norm->dN_ds == NULL`, i.e. the norm does not depend on the solution, the redundant 0.0 entries
	//       are not included in `dopt_t_ds` and hence also not in `lhs_opt_t`, resulting in the condition below to
	//       determine `ext_1`.
	const ptrdiff_t size_s_coef = compute_size(s_vol->sol_coef->order,s_vol->sol_coef->extents),
	                ext_0       = ( norm->dN_ds ? lhs_std->ext_1 : size_s_coef ),
	                ext_1       = (lhs_opt_t->ext_0)/ext_0;
	struct Matrix_d lhs_opt_t_M =
		{ .layout = 'C', .ext_0 = ext_0, .ext_1 = ext_1, .owns_data = false, .data = (Type*)lhs_opt_t->data, };
	transpose_Matrix_d(&lhs_opt_t_M,true);
	set_block_Matrix_d(lhs_opt,0,0,(struct const_Matrix_d*)&lhs_opt_t_M,0,0,
	                   lhs_opt_t_M.ext_0,lhs_opt_t_M.ext_1,'a');

	destructor_const_Vector_d(lhs_opt_t);
}

static void add_to_lhs_opt__l_mult
	(struct Matrix_d**const lhs_opt_ptr, const struct const_Matrix_d*const lhs_std,
	 const struct DPG_Solver_Volume*const dpg_s_vol, const struct Simulation*const sim)
{
	if (!test_case_explicitly_enforces_conservation(sim))
		return;

	struct Test_Case* test_case = (struct Test_Case*)sim->test_case_rc->tc;

	const struct const_Matrix_d*const lhs_l_mult_M = constructor_l_mult_M(lhs_std,dpg_s_vol,sim); // destructed
	assert(lhs_l_mult_M->ext_1 == test_case->n_eq);

	struct Matrix_d* lhs_opt = *lhs_opt_ptr;
	assert(lhs_opt->ext_0 == lhs_opt->ext_1);

	const ptrdiff_t ext_1_i = lhs_opt->ext_1,
	                ext_1   = ext_1_i + lhs_l_mult_M->ext_1;

	struct Matrix_d* lhs_add = constructor_zero_Matrix_d('R',ext_1,ext_1); // moved
	set_block_Matrix_d(lhs_add,0,0,(struct const_Matrix_d*)lhs_opt,0,0,lhs_opt->ext_0,lhs_opt->ext_1,'i');
	set_block_Matrix_d(lhs_add,0,ext_1_i,lhs_l_mult_M,0,0,lhs_l_mult_M->ext_0,lhs_l_mult_M->ext_1,'i');
	transpose_Matrix_d((struct Matrix_d*)lhs_l_mult_M,false);
	set_block_Matrix_d(lhs_add,ext_1_i,0,lhs_l_mult_M,0,0,lhs_l_mult_M->ext_0,lhs_l_mult_M->ext_1,'i');

	destructor_Matrix_d(lhs_opt);
	lhs_opt = lhs_add;

	destructor_const_Matrix_d(lhs_l_mult_M);

	*lhs_opt_ptr = lhs_add;

/// \todo Separate function after implementing
	if (test_case->is_linear)
		return;
	EXIT_ADD_SUPPORT;
}

static void add_to_petsc_Mat_Vec_dpg
	(const struct Solver_Volume* s_vol, const struct const_Vector_d* rhs_neg, const struct const_Matrix_d* lhs,
	 struct Solver_Storage_Implicit* ssi, const struct Simulation* sim)
{
	assert(sizeof(int) == sizeof(PetscInt));
	assert(sizeof(double) == sizeof(PetscScalar));

	const ptrdiff_t ext_0 = rhs_neg->ext_0;

	const struct const_Vector_i* idxm = constructor_petsc_idxm_dpg(ext_0,s_vol,sim); // destructed.

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

static const struct const_Matrix_d* constructor_norm_DPG_dN_ds__h1_upwind
	(const struct DPG_Solver_Volume* dpg_s_vol, const struct Flux_Ref* flux_r, const struct const_Matrix_d* n1_lt,
	 const struct Simulation* sim)
{
	struct Test_Case* test_case = (struct Test_Case*)sim->test_case_rc->tc;
	if (test_case->is_linear)
		return NULL;

/// \todo Combine common code from constructor_norm_DPG__h1_upwind into external function.
	const int n_eq = test_case->n_eq,
	          n_vr = test_case->n_var;

	const struct Multiarray_Operator cvcv1_vt_vc = get_operator__cvcv1_vt_vc__rlhs(dpg_s_vol);

	const ptrdiff_t ext_0 = cvcv1_vt_vc.data[0]->op_std->ext_0,
	                ext_1 = cvcv1_vt_vc.data[0]->op_std->ext_1;
	const ptrdiff_t n_dof_t = (n1_lt->ext_1)/n_eq,
	                n_dof_s = ext_1/n_dof_t;
	assert(n_dof_s == ((struct Solver_Volume*)dpg_s_vol)->sol_coef->extents[0]);

	struct Matrix_d* cvcv1r = constructor_empty_Matrix_d('C',n_vr*ext_0,n_eq*ext_1); // destructed

	struct Matrix_d* cvcv1r_l = constructor_empty_Matrix_d('C',ext_0,ext_1); // destructed
	const struct const_Multiarray_d* d2fr_ds2_Ma = flux_r->d2fr_ds2;
	struct Vector_d d2fr_ds2 = { .ext_0 = d2fr_ds2_Ma->extents[0], .owns_data = false, .data = NULL, };

	// dN_ds should be interpreted as a (neq*n_dof_t)x(neq*n_dof_t)x(n_vr*n_dof_s) 3-tensor.
	struct Matrix_d* dN_ds   = constructor_empty_Matrix_d('R',(n_eq*n_dof_t)*(n_vr*n_dof_s),n_eq*n_dof_t); // returned
	struct Matrix_d* dN_ds_l = constructor_empty_Matrix_d('R',n_dof_t,n_dof_t); // destructed
	for (int vr2 = 0; vr2 < n_vr; ++vr2) {
		for (int vr = 0; vr < n_vr; ++vr) {
		for (int eq = 0; eq < n_eq; ++eq) {
			set_to_value_Matrix_d(cvcv1r_l,0.0);
			for (int dim = 0; dim < DIM; ++dim) {
				const ptrdiff_t ind = compute_index_sub_container(
					d2fr_ds2_Ma->order,1,d2fr_ds2_Ma->extents,(ptrdiff_t[]){eq,vr,vr2,dim});
				d2fr_ds2.data = (Type*)&d2fr_ds2_Ma->data[ind];
				mm_diag_d('L',1.0,1.0,cvcv1_vt_vc.data[dim]->op_std,
				          (struct const_Vector_d*)&d2fr_ds2,cvcv1r_l,false);
			}
			set_block_Matrix_d(cvcv1r,eq*ext_0,vr*ext_1,
			                   (struct const_Matrix_d*)cvcv1r_l,0,0,cvcv1r_l->ext_0,cvcv1r_l->ext_1,'i');
		}}

/// \todo Add link to pdf documentation writing this out clearly.
		/* The slice of dN_ds corresponding to the 2nd contribution from the chain-rule differentiation of the
		 * df_ds'*df_ds term in the norm, where both f and s are vectors and where "'" denotes transposition.
		 * Explicitly:
		 *
		 * Given:
		 *   d/ds_{vr2} (df_ds'*df_ds) = d^2f_{ds ds_{vr2}}'*df_ds + df_ds'*d^2f_{ds ds_{vr2}}
		 *
		 * `dn1_ds` is the contribution from the df_ds'*d^2f_{ds ds_{vr2}} portion.
		 *
		 *  The 1st contribution can then be obtained from the transpose of the opposite term. Once again,
		 *  explicitly:
		 *
		 *  Considering the eq=0, vr=1 contribution of d/ds_{vr2} (df_ds'*df_ds), the 1st contribution is given by
		 *  the transpose of the eq=1, vr=0 contribution of d/ds_{vr2} (df_ds'*df_ds).
		 */
		const struct const_Matrix_d* dn1_ds =
			constructor_mm_const_Matrix_T('T','N',1.0,n1_lt,(struct const_Matrix_T*)cvcv1r,'R'); // destructed
		for (int vr = 0; vr < n_vr; ++vr) {
		for (int eq = 0; eq < n_eq; ++eq) {
		for (int dof_s = 0; dof_s < n_dof_s; ++dof_s) {
			const ptrdiff_t row_sub[] = { n_dof_t*(vr), n_dof_t*(eq), },
			                col_sub[] = { n_dof_t*(dof_s+n_dof_s*(eq)), n_dof_t*(dof_s+n_dof_s*(vr)), };

			set_block_Matrix_d(dN_ds_l,0,0,dn1_ds,row_sub[0],col_sub[0],dN_ds_l->ext_0,dN_ds_l->ext_1,'i');
			transpose_Matrix_d(dN_ds_l,false);
			set_block_Matrix_d(dN_ds_l,0,0,dn1_ds,row_sub[1],col_sub[1],dN_ds_l->ext_0,dN_ds_l->ext_1,'a');

			const ptrdiff_t row = n_dof_t*(eq+n_eq*(dof_s+n_dof_s*(vr2))),
			                col = n_dof_t*(vr);
			set_block_Matrix_d(dN_ds,row,col,
			                   (struct const_Matrix_d*)dN_ds_l,0,0,dN_ds_l->ext_0,dN_ds_l->ext_1,'i');
		}}}
		destructor_const_Matrix_d(dn1_ds);
	}
	destructor_Matrix_d(cvcv1r);
	destructor_Matrix_d(cvcv1r_l);
	destructor_Matrix_d(dN_ds_l);

	return (struct const_Matrix_d*) dN_ds;
}

// Level 1 ********************************************************************************************************** //

/// \brief Version of \ref add_to_dlhs_ds__face_boundary using the complex step method.
static void add_to_dlhs_ds__face_boundary_cmplx_step
	(struct Matrix_d*const dlhs_ds,             ///< See brief.
	 const struct DPG_Solver_Volume* dpg_s_vol, ///< See brief.
	 struct Simulation*const sim_c              ///< See brief.
	);

static const struct Operator* get_operator__cvcv0_vs_vc (const struct DPG_Solver_Volume* dpg_s_vol)
{
	const struct Volume* vol           = (struct Volume*) dpg_s_vol;
	const struct Solver_Volume* s_vol  = (struct Solver_Volume*) dpg_s_vol;
	const struct DPG_Solver_Element* e = (struct DPG_Solver_Element*) vol->element;

	const int p = s_vol->p_ref,
	          curved = vol->curved;
	return get_Multiarray_Operator(e->cvcv0_vs_vc[curved],(ptrdiff_t[]){0,0,p,p});
}

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
	struct Matrix_d* lhs_l  = constructor_empty_Matrix_d('C',n_dof_t,n_dof_s2);            // destructed
	struct Matrix_d* lhs    = constructor_empty_Matrix_d('C',n_eq*n_dof_t,n_vr2*n_dof_s2); // returned

	const struct const_Multiarray_d* d2fr_ds2_Ma = flux_r->d2fr_ds2;
	struct Vector_d d2fr_ds2 = { .ext_0 = d2fr_ds2_Ma->extents[0], .owns_data = false, .data = NULL, };

	/** \note: Many entries in the Hessian are 0 and need not be included in the loop below.
	 *  \todo Include a list of non-zero entries in the Hessian term (as part of \ref Flux_Ref_T) and skip
	 *        unnecessary terms below. Change to constructor_zero_\*. Also do this for the Jacobian terms. */
	for (int vr2 = 0; vr2 < n_vr; ++vr2) {
	for (int vr = 0; vr < n_vr; ++vr) {
	for (int eq = 0; eq < n_eq; ++eq) {
		set_to_value_Matrix_d(tw1_r,0.0);
		for (int dim = 0; dim < DIM; ++dim) {
			const ptrdiff_t ind = compute_index_sub_container(
				d2fr_ds2_Ma->order,1,d2fr_ds2_Ma->extents,(ptrdiff_t[]){eq,vr,vr2,dim});
			d2fr_ds2.data = (double*)&d2fr_ds2_Ma->data[ind];
			mm_diag_d('R',1.0,1.0,tw1_vt_vc.data[dim]->op_std,(struct const_Vector_d*)&d2fr_ds2,tw1_r,false);
		}
		mm_d('N','N',1.0,0.0,(struct const_Matrix_d*)tw1_r,cvcv0_vs_vc->op_std,lhs_l);

		for (int dof_s = 0; dof_s < n_dof_s; ++dof_s) {
			const ptrdiff_t row = eq*n_dof_t,
			                col = (n_dof_s*n_vr*vr2+dof_s*n_vr+vr)*n_dof_s;
			const double*const data = (const double*) get_col_Matrix_d(dof_s*n_dof_s,lhs_l);
			const struct const_Matrix_d lhs_ll =
				{ .layout = lhs_l->layout, .ext_0 = n_dof_t, .ext_1 = n_dof_s, .data = data, };
			set_block_Matrix_d(lhs,row,col,&lhs_ll,0,0,lhs_ll.ext_0,lhs_ll.ext_1,'i');
		}
	}}}
	destructor_Matrix_d(tw1_r);
	destructor_Matrix_d(lhs_l);

	return lhs;
}

static void add_to_dlhs_ds__norm
	(struct Matrix_d**const dlhs_ds_ptr, const struct Norm_DPG* norm, const struct const_Matrix_d* opt_t)
{
	if (!norm->dN_ds)
		return;

	struct Matrix_d* dlhs_ds = *dlhs_ds_ptr;

	const char layout = dlhs_ds->layout;
	assert(layout == 'C');
	struct Matrix_d* dlhs_ds_total = constructor_mm_Matrix_T('N','N',-1.0,norm->dN_ds,opt_t,'R'); // moved

	const ptrdiff_t n_t  = norm->N->ext_0,
	                n_s  = (norm->dN_ds->ext_0)/n_t,
	                n_sf = opt_t->ext_1;

	// Transpose sub-blocks and reinterpret `dlhs_ds_total` as having different extents.
	for (int s = 0; s < n_s; ++s) {
		const ptrdiff_t size_block = n_t*n_sf;
		struct Matrix_d sub_block = { .layout = dlhs_ds_total->layout, .ext_0 = n_t, .ext_1 = n_sf,
		                              .owns_data = false, .data = &dlhs_ds_total->data[s*size_block], };
		transpose_Matrix_d(&sub_block,true);
	}
	dlhs_ds_total->layout = 'C';
	dlhs_ds_total->ext_0 = n_t;
	dlhs_ds_total->ext_1 = n_sf*n_s;

	// As the `dlhs_ds` input does not currently have any contribution from the face unknowns, its contributions are
	// added to `dlhs_ds_total` in the locations corresponding to the solution dof.
	for (int s = 0; s < n_s; ++s) {
		const ptrdiff_t col_sf = s*n_sf,
		                col_s  = s*n_s;

		set_block_Matrix_d(dlhs_ds_total,0,col_sf,(struct const_Matrix_d*)dlhs_ds,0,col_s,n_t,n_s,'a');
	}

	destructor_Matrix_d(dlhs_ds);
	*dlhs_ds_ptr = dlhs_ds_total;

	assert((*dlhs_ds_ptr)->layout == layout);
}

static void add_to_dlhs_ds__face_boundary
	(struct Matrix_d*const dlhs_ds, const struct DPG_Solver_Volume* dpg_s_vol, const struct Simulation*const sim,
	 struct Simulation*const sim_c, const char lin_method)
{
UNUSED(sim);
	if (lin_method == 'a')
		EXIT_ADD_SUPPORT; // Will require linearization of boundary condition/numerical flux Jacobians.
	else if (lin_method == 'c')
		add_to_dlhs_ds__face_boundary_cmplx_step(dlhs_ds,dpg_s_vol,sim_c);
	else
		EXIT_ERROR("Unsupported: %c\n",lin_method);
}

static struct Multiarray_Operator get_operator__cvcv1_vt_vc__rlhs (const struct DPG_Solver_Volume* dpg_s_vol)
{
	struct Volume* vol          = (struct Volume*) dpg_s_vol;
	struct Solver_Volume* s_vol = (struct Solver_Volume*) vol;

	const struct DPG_Solver_Element* dpg_s_e = (struct DPG_Solver_Element*) vol->element;

	const int p      = s_vol->p_ref,
	          curved = vol->curved;

	return set_MO_from_MO(dpg_s_e->cvcv1_vt_vc[curved],1,(ptrdiff_t[]){0,0,p,p});
}

// Level 2 ********************************************************************************************************** //

/// \brief Constructor for the complex \ref Simulation volumes and faces DPG computational element lists.
static void constructor_Simulation_c_comp_elems
	(struct Simulation*const sim_c,            ///< The complex \ref Simulation.
	 const struct DPG_Solver_Volume* dpg_s_vol ///< The current \ref DPG_Solver_Volume_T.
	);

/// \brief Destructor for the complex \ref Simulation volumes and faces DPG computational element lists.
static void destructor_Simulation_c_comp_elems
	(struct Simulation*const sim_c ///< The complex \ref Simulation.
	);

static void add_to_dlhs_ds__face_boundary_cmplx_step
	(struct Matrix_d*const dlhs_ds, const struct DPG_Solver_Volume* dpg_s_vol, struct Simulation*const sim_c)
{
	constructor_Simulation_c_comp_elems(sim_c,dpg_s_vol); // destructed

	const struct Solver_Volume_c* s_vol_c =  (struct Solver_Volume_c*) sim_c->volumes->first;
	struct Multiarray_c*const s_coef_c = s_vol_c->sol_coef;

	const ptrdiff_t size_s = compute_size(s_coef_c->order,s_coef_c->extents);

	struct Matrix_c* lhs = constructor_empty_Matrix_c('R',dlhs_ds->ext_0,size_s); // destructed
	for (int col_l = 0; col_l < size_s; ++col_l) {
		set_to_value_Matrix_c(lhs,0.0);

		s_coef_c->data[col_l] += CX_STEP*I;
		add_to_rlhs__face_c(NULL,&lhs,(struct DPG_Solver_Volume_c*)s_vol_c,sim_c,false);
		s_coef_c->data[col_l] -= CX_STEP*I;

		transpose_Matrix_c(lhs,true);
		set_block_Matrix_d_cmplx_step(dlhs_ds,(struct const_Matrix_c*)lhs,0,col_l*size_s,'a');

		// Not transposed back as memory is zero'ed => Only need to change `layout`.
		lhs->layout = 'R';
	}
	destructor_Matrix_c(lhs);

	destructor_Simulation_c_comp_elems(sim_c);
}

// Level 3 ********************************************************************************************************** //

/** \brief Constructor for a list of volumes including only a copy of the current volume.
 *  \return See brief. */
static struct Intrusive_List* constructor_Volumes_dpg_local
	(const struct DPG_Solver_Volume*const dpg_s_vol, ///< The current \ref DPG_Solver_Volume_T.
	 struct Simulation*const sim                     ///< The complex \ref Simulation.
	);

/** \brief Constructor for a list of faces including only copies of boundary faces neighbouring the current volume.
 *  \return See brief. */
static struct Intrusive_List* constructor_Faces_dpg_local
	(const struct DPG_Solver_Volume*const dpg_s_vol, ///< The current \ref DPG_Solver_Volume_T.
	 struct Simulation*const sim                     ///< The complex \ref Simulation.
	);

/** \brief Copy members of real to complex computational elements of the current type as determined by the name of the
 *         \ref Simulation::volumes list. */
static void copy_members_computational_elements_dpg
	(const struct DPG_Solver_Volume*const dpg_s_vol, ///< The current \ref DPG_Solver_Volume_T.
	 const struct Simulation*const sim               ///< The complex \ref Simulation.
	);

static void constructor_Simulation_c_comp_elems
	(struct Simulation*const sim_c, const struct DPG_Solver_Volume* dpg_s_vol)
{
	assert(sim_c->volumes == NULL && sim_c->faces == NULL);
	sim_c->volumes = constructor_Volumes_dpg_local(dpg_s_vol,sim_c);
	sim_c->faces   = constructor_Faces_dpg_local(dpg_s_vol,sim_c);

	constructor_derived_computational_elements_c(sim_c,IL_SOLVER);
	copy_members_computational_elements_dpg(dpg_s_vol,sim_c);

	constructor_derived_computational_elements_c(sim_c,IL_SOLVER_DPG);
	copy_members_computational_elements_dpg(dpg_s_vol,sim_c);
}

static void destructor_Simulation_c_comp_elems (struct Simulation*const sim_c)
{
	destructor_derived_computational_elements_c(sim_c,IL_SOLVER);
	destructor_derived_computational_elements_c(sim_c,IL_BASE);

	destructor_Volumes(sim_c->volumes);
	sim_c->volumes = NULL;

	destructor_Faces(sim_c->faces);
	sim_c->faces = NULL;
}

// Level 4 ********************************************************************************************************** //

/** \brief Return the pointer to the \ref DPG_Solver_Face_T of the input \ref DPG_Solver_Volume_T with the given index.
 *  \return See brief. */
struct DPG_Solver_Face* get_dpg_face_ptr
	(const struct DPG_Solver_Volume*const dpg_s_vol, ///< The current volume.
	 const int index_f                               ///< The index of the face.
	);

static struct Intrusive_List* constructor_Volumes_dpg_local
	(const struct DPG_Solver_Volume*const dpg_s_vol, struct Simulation*const sim)
{
	struct Intrusive_List* volumes = constructor_empty_IL(IL_VOLUME,NULL); // returned

	push_back_IL(volumes,(struct Intrusive_Link*)constructor_copy_Volume((struct Volume*)dpg_s_vol,sim,false));
	sim->n_v = 1;

	return volumes;
}

static struct Intrusive_List* constructor_Faces_dpg_local
	(const struct DPG_Solver_Volume*const dpg_s_vol, struct Simulation*const sim)
{
	struct Intrusive_List* faces = constructor_empty_IL(IL_FACE,NULL); // returned

	const struct Volume*const vol_master = (struct Volume*) dpg_s_vol;

	struct Volume* vol = (struct Volume*) sim->volumes->first;
	assert(sim->volumes->first->next == NULL);

	ptrdiff_t n_f = 0;
	for (int i = 0; i < NFMAX;    ++i) {
	for (int j = 0; j < NSUBFMAX; ++j) {
		const struct Face* face_master = vol_master->faces[i][j];
		if (!face_master || !face_master->boundary)
			continue;

		push_back_IL(faces,(struct Intrusive_Link*)constructor_copy_Face(face_master,sim,false));
		++n_f;

		struct Face* face = (struct Face*) faces->last;

		const_cast_Face(&vol->faces[i][j],face);
		face->neigh_info[0].volume = vol;
	}}
	sim->n_f = n_f;
	assert(n_f < NFMAX);

	return faces;
}

static void copy_members_computational_elements_dpg
	(const struct DPG_Solver_Volume*const dpg_s_vol, const struct Simulation*const sim)
{
	int category = -1;
	if (sim->volumes->name == IL_VOLUME_SOLVER) {
		assert(sim->faces->name == IL_FACE_SOLVER);
		category = IL_SOLVER;
	} else if (sim->volumes->name == IL_VOLUME_SOLVER_DPG) {
		assert(sim->faces->name == IL_FACE_SOLVER_DPG);
		category = IL_SOLVER_DPG;
	} else {
		EXIT_ERROR("Unsupported: %d\n",sim->volumes->name);
	}

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		if (category == IL_SOLVER)
			copy_members_r_to_c_Solver_Volume((struct Solver_Volume_c*)curr,(struct Solver_Volume*)dpg_s_vol,sim);
		else if (category == IL_SOLVER_DPG)
			copy_members_r_to_c_DPG_Solver_Volume((struct DPG_Solver_Volume_c*)curr,dpg_s_vol,sim);
	}

	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		const struct Face*const face = (struct Face*) curr;
		const struct DPG_Solver_Face*const dpg_s_face_r = get_dpg_face_ptr(dpg_s_vol,face->index);
		if (category == IL_SOLVER)
			copy_members_r_to_c_Solver_Face((struct Solver_Face_c*)curr,(struct Solver_Face*)dpg_s_face_r,sim);
		else if (category == IL_SOLVER_DPG)
			copy_members_r_to_c_DPG_Solver_Face((struct DPG_Solver_Face_c*)curr,dpg_s_face_r,sim);
	}
}

// Level 5 ********************************************************************************************************** //

struct DPG_Solver_Face* get_dpg_face_ptr (const struct DPG_Solver_Volume*const dpg_s_vol, const int index_f)
{
	const struct Volume*const vol = (struct Volume*) dpg_s_vol;
	for (int i = 0; i < NFMAX;    ++i) {
	for (int j = 0; j < NSUBFMAX; ++j) {
		const struct Face* face = vol->faces[i][j];
		if (!face || (face->index != index_f))
			continue;

		return (struct DPG_Solver_Face*) face;
	}}
	EXIT_ERROR("Did not find the face.\n");
}
