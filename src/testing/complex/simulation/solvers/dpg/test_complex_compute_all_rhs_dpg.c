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

#include "test_complex_compute_all_rhs_dpg.h"

#include <assert.h>
#include "petscmat.h"

#include "macros.h"
#include "definitions_dpg.h"
#include "definitions_intrusive.h"
#include "definitions_test_integration.h"

#include "test_complex_flux.h"
#include "test_complex_compute_face_rhs.h"
#include "test_complex_compute_volume_rhs.h"
#include "test_complex_numerical_flux.h"
#include "test_complex_operators.h"
#include "test_complex_test_case.h"

#include "face_solver_dpg_complex.h"
#include "volume_solver_dpg_complex.h"

#include "complex_matrix.h"
#include "complex_multiarray.h"
#include "complex_vector.h"
#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "compute_all_rlhs_dpg.h"
#include "compute_face_rlhs.h"
#include "compute_volume_rlhs.h"
#include "intrusive.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "solve.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief `complex` version of \ref constructor_Numerical_Flux_Input_data for the dpg scheme.
static void constructor_Numerical_Flux_Input_c_data_dpg
	(struct Numerical_Flux_Input_c* num_flux_i,         ///< See brief.
	 const struct Complex_DPG_Solver_Face* c_dg_s_face, ///< See brief.
	 const struct Simulation* sim                       ///< See brief.
	);

#include "def_templates_type_dc.h"

#include "def_templates_matrix_c.h"
#include "def_templates_multiarray_c.h"
#include "def_templates_vector_c.h"

#include "def_templates_compute_all_rlhs_dpg_c.h"
#include "def_templates_compute_face_rlhs_c.h"
#include "def_templates_multiarray_operator_c.h"
#include "def_templates_numerical_flux_c.h"

#include "compute_all_rlhs_dpg_T.c"

/** \brief See \ref constructor_norm_op_fptr.
 *
 *  \param c_dpg_s_vol See brief.
 *  \param flux_r      See brief.
 *  \param sim         See brief.
 */
typedef const struct const_Matrix_c* (*constructor_norm_op_c_fptr)
	(const struct Complex_DPG_Solver_Volume* c_dpg_s_vol,
	 const struct Flux_Ref_c* flux_r,
	 const struct Simulation* sim
	);

/// \brief Container for solver-related parameters.
struct S_Params_DPG_c {
	constructor_norm_op_c_fptr constructor_norm_op; ///< Pointer to the appropriate function.
//	compute_rhs_c_fptr compute_rlhs;                ///< Pointer to the appropriate function.
};

/** \brief Set the parameters of \ref S_Params_DPG_c.
 *  \return A statically allocated \ref S_Params_DPG_c container. */
static struct S_Params_DPG_c set_s_params_dpg_c
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief `complex` version of \ref compute_rlhs_1 computing rhs terms only.
static void compute_rhs_1_c
	(const struct const_Matrix_c* norm_op,                ///< See brief.
	 const struct Flux_Ref_c* flux_r,                     ///< See brief.
	 const struct Complex_DPG_Solver_Volume* c_dpg_s_vol, ///< See brief.
	 struct Solver_Storage_Implicit* ssi,                 ///< See brief.
	 const struct Simulation* sim                         ///< See brief.
	);

// Interface functions ********************************************************************************************** //

void compute_all_rhs_dpg_c
	(const struct Complex_DPG_Solver_Volume* c_dpg_s_vol, struct Solver_Storage_Implicit* ssi,
	 const struct Simulation* sim)
{
	assert(sim->volumes->name == IL_VOLUME_SOLVER_DPG_COMPLEX);
	assert(sim->elements->name == IL_ELEMENT_SOLVER_DPG);

	/** \note Linearized terms are required for the computation of optimal test functions, even if only computing rhs
	 *        terms. */
	assert(sim->test_case->solver_method_curr == 'i');

	struct Flux_Input_c* flux_i = constructor_Flux_Input_c(sim); // destructed

	struct S_Params_Volume_Structor_c spvs = { .constructor_sol_vc = constructor_sol_vc_dpg_c, };
	struct Solver_Volume* s_vol = (struct Solver_Volume*) c_dpg_s_vol;
	struct Flux_Ref_c* flux_r = constructor_Flux_Ref_vol_c(&spvs,flux_i,s_vol,sim); // destructed

	struct S_Params_DPG_c s_params = set_s_params_dpg_c(sim);
	const struct const_Matrix_c* norm_op = s_params.constructor_norm_op(c_dpg_s_vol,flux_r,sim); // destructed

	assert(sim->test_case->has_2nd_order == false); // add support.
	compute_rhs_1_c(norm_op,flux_r,c_dpg_s_vol,ssi,sim);

	destructor_const_Matrix_c(norm_op);
	destructor_Flux_Ref_c(flux_r);

	destructor_Flux_Input_c(flux_i);
}

const struct const_Multiarray_c* constructor_sol_vc_dpg_c
	(const struct Solver_Volume* s_vol, const struct Simulation* sim)
{
	UNUSED(sim);
	const struct Operator* cv0_vs_vc = get_operator__cv0_vs_vc(s_vol);

	struct Complex_DPG_Solver_Volume* s_vol_c = (struct Complex_DPG_Solver_Volume*) s_vol;
	const struct const_Multiarray_c* s_coef = (const struct const_Multiarray_c*) s_vol_c->sol_coef;

	return constructor_mm_NN1_Operator_const_Multiarray_c(cv0_vs_vc,s_coef,'C','d',s_coef->order,NULL);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief `complex` version of \ref constructor_norm_op__h1_upwind.
 *  \return See brief. */
static const struct const_Matrix_c* constructor_norm_op__h1_upwind_c
	(const struct Complex_DPG_Solver_Volume* c_dpg_s_vol, ///< See brief.
	 const struct Flux_Ref_c* flux_r,                     ///< See brief.
	 const struct Simulation* sim                         ///< See brief.
	);

/** \brief `complex` version of \ref constructor_rhs_v_1.
 *  \return See brief. */
static struct Vector_c* constructor_rhs_v_1_c
	(const struct Flux_Ref_c* flux_r,   ///< See brief.
	 const struct Solver_Volume* s_vol, ///< See brief.
	 const struct Simulation* sim       ///< See brief.
	);

/** \brief `complex` version of \ref constructor_lhs_v_1.
 *  \return See brief. */
struct Matrix_c* constructor_lhs_v_1_c
	(const struct Flux_Ref_c* flux_r,   ///< See brief.
	 const struct Solver_Volume* s_vol, ///< See brief.
	 const struct Simulation* sim       ///< See brief.
	);

/** \brief Version of \ref add_to_petsc_Mat_Vec_dpg setting a single column of Solver_Storage_Implicit::A to the values
 *         computed using the complex step rhs. */
static void add_to_petsc_Mat_dpg_c
	(const struct Solver_Volume* s_vol,    ///< See brief.
	 const struct const_Vector_c* rhs_neg, ///< See brief.
	 struct Solver_Storage_Implicit* ssi   ///< See brief.
	);

static struct S_Params_DPG_c set_s_params_dpg_c (const struct Simulation* sim)
{
	struct S_Params_DPG_c s_params;

	struct Test_Case* test_case = sim->test_case;
	switch (test_case->ind_test_norm) {
		case TEST_NORM_H1_UPWIND: s_params.constructor_norm_op = constructor_norm_op__h1_upwind_c; break;
		default:                  EXIT_ERROR("Unsupported: %d\n",test_case->ind_test_norm);        break;
	}

	return s_params;
}

static void compute_rhs_1_c
	(const struct const_Matrix_c* norm_op, const struct Flux_Ref_c* flux_r,
	 const struct Complex_DPG_Solver_Volume* c_dpg_s_vol, struct Solver_Storage_Implicit* ssi,
	 const struct Simulation* sim)
{
	const struct Solver_Volume* s_vol = (struct Solver_Volume*) c_dpg_s_vol;
	struct Vector_c* rhs = constructor_rhs_v_1_c(flux_r,s_vol,sim); // destructed
	struct Matrix_c* lhs = constructor_lhs_v_1_c(flux_r,s_vol,sim); // destructed

	const struct DPG_Solver_Volume* dpg_s_vol = (struct DPG_Solver_Volume*) c_dpg_s_vol;
	increment_and_add_dof_rlhs_f_1_c(rhs,&lhs,dpg_s_vol,sim);

	const struct const_Matrix_c* optimal_test =
		constructor_sysv_const_Matrix_c(norm_op,(struct const_Matrix_c*)lhs); // destructed
	destructor_Matrix_c(lhs);

	const struct const_Vector_c* rhs_opt =
		constructor_mv_const_Vector_c('T',-1.0,optimal_test,(struct const_Vector_c*)rhs); // destructed
	destructor_Vector_c(rhs);
	destructor_const_Matrix_c(optimal_test);

	add_to_petsc_Mat_dpg_c(s_vol,rhs_opt,ssi);

	destructor_const_Vector_c(rhs_opt);
}

// Level 1 ********************************************************************************************************** //

static const struct const_Matrix_c* constructor_norm_op__h1_upwind_c
	(const struct Complex_DPG_Solver_Volume* c_dpg_s_vol, const struct Flux_Ref_c* flux_r,
	 const struct Simulation* sim)
{
	const int n_eq = sim->test_case->n_eq,
	          n_vr = sim->test_case->n_var;

	struct DPG_Solver_Volume* dpg_s_vol = (struct DPG_Solver_Volume*) c_dpg_s_vol;
	const struct Multiarray_Operator cvt1_vt_vc = get_operator__cvt1_vt_vc__rlhs(dpg_s_vol);

	const ptrdiff_t ext_0 = cvt1_vt_vc.data[0]->op_std->ext_0,
	                ext_1 = cvt1_vt_vc.data[0]->op_std->ext_1;

	struct Matrix_c* cvt1r = constructor_empty_Matrix_c('R',n_eq*ext_0,n_vr*ext_1); // destructed

	struct Matrix_c* cvt1r_l = constructor_empty_Matrix_c('R',ext_0,ext_1); // destructed
	const struct const_Multiarray_c* dfr_ds_Ma = flux_r->dfr_ds;
	struct Vector_c dfr_ds = { .ext_0 = dfr_ds_Ma->extents[0], .owns_data = false, .data = NULL, };

	for (int vr = 0; vr < n_vr; ++vr) {
	for (int eq = 0; eq < n_eq; ++eq) {
		set_to_value_Matrix_c(cvt1r_l,0.0);
		for (int dim = 0; dim < DIM; ++dim) {
			const ptrdiff_t ind =
				compute_index_sub_container(dfr_ds_Ma->order,1,dfr_ds_Ma->extents,(ptrdiff_t[]){eq,vr,dim});
			dfr_ds.data = (double complex*)&dfr_ds_Ma->data[ind];
			mm_diag_c('R',1.0,1.0,cvt1_vt_vc.data[dim]->op_std,(struct const_Vector_c*)&dfr_ds,cvt1r_l,false);
		}
		set_block_Matrix_c(cvt1r,(struct const_Matrix_c*)cvt1r_l,eq*ext_0,vr*ext_1,'i');
	}}
	destructor_Matrix_c(cvt1r_l);

	struct Solver_Volume* s_vol = (struct Solver_Volume*) c_dpg_s_vol;
	const struct const_Vector_d* w_vc = get_operator__w_vc__s_e(s_vol);
	const struct const_Vector_d J_vc  = interpret_const_Multiarray_as_Vector_d(s_vol->jacobian_det_vc);

	const struct const_Vector_d* J_inv_vc = constructor_inverse_const_Vector_d(&J_vc);               // destructed
	const struct const_Vector_d* wJ_vc    = constructor_dot_mult_const_Vector_d(w_vc,J_inv_vc,n_vr); // destructed
	destructor_const_Vector_d(J_inv_vc);

	const struct const_Matrix_c* n1_l =
		constructor_mm_diag_const_Matrix_c_d(1.0,(struct const_Matrix_c*)cvt1r,wJ_vc,'R',false); // destructed
	destructor_const_Vector_d(wJ_vc);

	const struct const_Matrix_c* n1 =
		constructor_mm_const_Matrix_c('N','T',1.0,n1_l,(struct const_Matrix_c*)cvt1r,'R'); // destructed
	destructor_const_Matrix_c(n1_l);
	destructor_Matrix_c(cvt1r);

	const struct const_Matrix_d* norm_op_H0 = dpg_s_vol->norm_op_H0;
	assert(norm_op_H0->ext_0 == ext_0);

	struct Matrix_c* norm_op = constructor_empty_Matrix_c('R',n_eq*ext_0,n_eq*ext_0); // returned

	set_block_Matrix_c(norm_op,n1,0,0,'i');
	for (int eq = 0; eq < n_eq; ++eq)
		set_block_Matrix_c_d(norm_op,norm_op_H0,eq*ext_0,eq*ext_0,'a');
	destructor_const_Matrix_c(n1);

	return (struct const_Matrix_c*) norm_op;
}

static struct Vector_c* constructor_rhs_v_1_c
	(const struct Flux_Ref_c* flux_r, const struct Solver_Volume* s_vol, const struct Simulation* sim)
{
	const int n_eq = sim->test_case->n_eq;

	const struct Multiarray_Operator tw1_vt_vc = get_operator__tw1_vt_vc(s_vol);

	const ptrdiff_t ext_0 = tw1_vt_vc.data[0]->op_std->ext_0;

	struct Vector_c* rhs = constructor_zero_Vector_c(ext_0*n_eq); // returned

	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	const char op_format = 'd';

	ptrdiff_t extents[2] = { ext_0, n_eq, };
	struct Multiarray_c rhs_Ma =
		{ .layout = 'C', .order = 2, .extents = extents, .owns_data = false, .data = rhs->data, };
	for (ptrdiff_t dim = 0; dim < DIM; ++dim)
		mm_NNC_Operator_Multiarray_c(1.0,1.0,tw1_vt_vc.data[dim],flux_r->fr,&rhs_Ma,op_format,2,&dim,NULL);

	return rhs;
}

struct Matrix_c* constructor_lhs_v_1_c
	(const struct Flux_Ref_c* flux_r, const struct Solver_Volume* s_vol, const struct Simulation* sim)
{
	const struct Multiarray_Operator tw1_vt_vc = get_operator__tw1_vt_vc(s_vol);
	const struct Operator* cv0_vs_vc = get_operator__cv0_vs_vc(s_vol);

	const ptrdiff_t ext_0 = tw1_vt_vc.data[0]->op_std->ext_0,
	                ext_1 = tw1_vt_vc.data[0]->op_std->ext_1;
	const int n_eq = sim->test_case->n_eq,
	          n_vr = sim->test_case->n_var;

	struct Matrix_c* tw1_r = constructor_empty_Matrix_c('R',ext_0,ext_1);                         // destructed
	struct Matrix_c* lhs_l = constructor_empty_Matrix_c('R',ext_0,cv0_vs_vc->op_std->ext_1);      // destructed
	struct Matrix_c* lhs   = constructor_empty_Matrix_c('R',n_eq*lhs_l->ext_0,n_vr*lhs_l->ext_1); // returned

	const struct const_Multiarray_c* dfr_ds_Ma = flux_r->dfr_ds;
	struct Vector_c dfr_ds = { .ext_0 = dfr_ds_Ma->extents[0], .owns_data = false, .data = NULL, };

	for (int vr = 0; vr < n_vr; ++vr) {
	for (int eq = 0; eq < n_eq; ++eq) {
		set_to_value_Matrix_c(tw1_r,0.0);
		for (int dim = 0; dim < DIM; ++dim) {
			const ptrdiff_t ind =
				compute_index_sub_container(dfr_ds_Ma->order,1,dfr_ds_Ma->extents,(ptrdiff_t[]){eq,vr,dim});
			dfr_ds.data = (double complex*)&dfr_ds_Ma->data[ind];
			mm_diag_c('R',1.0,1.0,tw1_vt_vc.data[dim]->op_std,(struct const_Vector_c*)&dfr_ds,tw1_r,false);
		}

		mm_cdc('N','N',1.0,0.0,(struct const_Matrix_c*)tw1_r,cv0_vs_vc->op_std,lhs_l);
		set_block_Matrix_c(lhs,(struct const_Matrix_c*)lhs_l,eq*lhs_l->ext_0,vr*lhs_l->ext_1,'i');
	}}
	destructor_Matrix_c(tw1_r);
	destructor_Matrix_c(lhs_l);

	return lhs;
}

static void add_to_petsc_Mat_dpg_c
	(const struct Solver_Volume* s_vol, const struct const_Vector_c* rhs_neg, struct Solver_Storage_Implicit* ssi)
{
	const ptrdiff_t ext_0 = rhs_neg->ext_0;

	const struct const_Vector_i* idxm = constructor_petsc_idxm_dpg(ext_0,s_vol); // destructed.

	PetscScalar rhs_c_data[ext_0];
	for (int i = 0; i < ext_0; ++i)
		rhs_c_data[i] = cimag((-rhs_neg->data[i])/CX_STEP);

	MatSetValues(ssi->A,(PetscInt)ext_0,idxm->data,1,&ssi->col,rhs_c_data,ADD_VALUES);

	destructor_const_Vector_i(idxm);
}

// Level 2 ********************************************************************************************************** //

static void constructor_Numerical_Flux_Input_c_data_dpg
	(struct Numerical_Flux_Input_c* num_flux_i, const struct Complex_DPG_Solver_Face* c_dpg_s_face,
	const struct Simulation* sim)
{
	struct Complex_Test_Case* test_case = (struct Complex_Test_Case*) sim->test_case;
	struct Solver_Face* s_face = (struct Solver_Face*) c_dpg_s_face;

	test_case->constructor_Boundary_Value_Input_c_face_fcl(&num_flux_i->bv_l,s_face,sim);           // destructed
	c_dpg_s_face->constructor_Boundary_Value_c_fcl(&num_flux_i->bv_r,&num_flux_i->bv_l,s_face,sim); // destructed
}
