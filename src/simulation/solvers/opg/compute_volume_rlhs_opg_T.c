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

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"
#include "definitions_intrusive.h"


#include "def_templates_compute_volume_rlhs_opg.h"

#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"

#include "def_templates_volume_solver_opg.h"

#include "def_templates_compute_volume_rlhs.h"
#include "def_templates_flux.h"
#include "def_templates_test_case.h"
#include "def_templates_operators.h"

// Static function declarations ************************************************************************************* //

/// \brief Container for solver-related parameters.
struct S_Params_T {
	struct S_Params_Volume_Structor_T spvs; ///< \ref S_Params_Volume_Structor_T.

	compute_rlhs_v_fptr_T compute_rlhs; ///< Pointer to the appropriate function.
};

/** \brief Set the parameters of \ref S_Params_T.
 *  \return A statically allocated \ref S_Params_T container. */
static struct S_Params_T set_s_params_T
	(const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Constructor for the 'test' function 'diff'erentiation 'op'erator for the '1'st order 'v'olume term for the
 *         OPG scheme.
 *  \return See brief.
 *
 *  This operator is the rhs portion of the operator returned from \ref constructor_operator__test_s_coef_to_sol_coef_T.
 */
static const struct const_Matrix_T* constructor_test_diff_op_1v_opg_T
	(const struct Flux_Ref_T*const flux_r,     ///< Standard.
	 const struct Solver_Volume_T*const s_vol, ///< Standard.
	 const bool include_det_j                  /**< Flag for whether the inverse Jacobian determinant should be
	                                            *   included. */
		);

// Interface functions ********************************************************************************************** //

void compute_volume_rlhs_opg_T
	(const struct Simulation* sim, struct Solver_Storage_Implicit* ssi, struct Intrusive_List* volumes)
{
	assert(sim->volumes->name == IL_VOLUME_SOLVER_OPG);
	assert(sim->elements->name == IL_ELEMENT_SOLVER_OPG);

	struct S_Params_T s_params = set_s_params_T(sim);
	struct Flux_Input_T* flux_i = constructor_Flux_Input_T(sim); // destructed

	for (struct Intrusive_Link* curr = volumes->first; curr; curr = curr->next) {
		struct Solver_Volume_T*const s_vol = (struct Solver_Volume_T*) curr;

		struct Flux_Ref_T* flux_r = constructor_Flux_Ref_vol_T(&s_params.spvs,flux_i,s_vol); // destructed

		// Compute the rhs and the lhs terms.
		s_params.compute_rlhs(flux_r,s_vol,ssi);
		destructor_Flux_Ref_T(flux_r);
	}
	destructor_Flux_Input_T(flux_i);
}

const struct const_Matrix_T* constructor_operator__test_s_coef_to_sol_coef_T
	(const struct Flux_Ref_T*const flux_r, const struct OPG_Solver_Volume_T*const opg_s_vol)
{
	struct Solver_Volume_T*const s_vol = (struct Solver_Volume_T*) opg_s_vol;

	const struct Operator*const cv0_vs_vc = get_operator__cv0_vs_vc_T(s_vol);
	const struct const_Matrix_T*const op_proj_L2 =
		constructor_mm_TR_const_Matrix_T('N','T',1.0,opg_s_vol->m_inv,cv0_vs_vc->op_std,'R'); // destructed

	const struct const_Vector_R*const w_vc = get_operator__w_vc__s_e_T(s_vol);
	scale_Matrix_T_by_Vector_R('R',1.0,(struct Matrix_T*)op_proj_L2,w_vc,false);

	// Note: Inverse Jacobian determinant from this operator cancels with that from the integral (both omitted).
	const struct const_Matrix_T*const op_v1_opg = constructor_test_diff_op_1v_opg_T(flux_r,s_vol,false); // dest.

	const int n_eq = get_set_n_var_eq(NULL)[1];

	const char layout = op_v1_opg->layout;
	const ptrdiff_t sub_ext[2] = { op_v1_opg->ext_0/n_eq, op_proj_L2->ext_0, };

	assert(layout == 'R');

	struct Matrix_T*const op_final = constructor_empty_Matrix_T(layout,n_eq*op_proj_L2->ext_0,op_v1_opg->ext_1); // r.
	struct Matrix_T sub_op_v1_opg = { .layout = layout, .ext_0 = sub_ext[0], .ext_1 = op_v1_opg->ext_1, .data = NULL, };
	struct Matrix_T sub_op_final  = { .layout = layout, .ext_0 = sub_ext[1], .ext_1 = op_v1_opg->ext_1, .data = NULL, };

	for (int eq = 0; eq < n_eq; ++eq) {
		sub_op_v1_opg.data = (Type*) get_row_const_Matrix_T(sub_ext[0]*eq,op_v1_opg);
		sub_op_final.data  = get_row_Matrix_T(sub_ext[1]*eq,op_final);
		mm_T('N','N',1.0,0.0,op_proj_L2,(struct const_Matrix_T*)&sub_op_v1_opg,&sub_op_final);
	}
	destructor_const_Matrix_T(op_proj_L2);
	destructor_const_Matrix_T(op_v1_opg);

	return (struct const_Matrix_T*) op_final;
}

void update_coef_s_v_opg_T (const struct Simulation*const sim, struct Intrusive_List*const volumes)
{
	/* The L2 projection is used to compute \ref Solver_Volume_T::sol_coef from \ref Solver_Volume_T::test_s_coef.
	 * For collocated schemes when the solution polynomial and test function polynomial degrees are equal, the
	 * projection operator reduces to identity. */

	struct Test_Case_T*const test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const char smc = test_case->solver_method_curr;
	test_case->solver_method_curr = 'i';

	struct Flux_Input_T*const flux_i = constructor_Flux_Input_T(sim); // destructed

	struct S_Params_T s_params;
	set_S_Params_Volume_Structor_T(&s_params.spvs,sim);

	for (struct Intrusive_Link* curr = volumes->first; curr; curr = curr->next) {
		struct Solver_Volume_T*const s_vol         = (struct Solver_Volume_T*) curr;
		struct OPG_Solver_Volume_T*const opg_s_vol = (struct OPG_Solver_Volume_T*) curr;

		struct Flux_Ref_T*const flux_r = constructor_Flux_Ref_vol_T(&s_params.spvs,flux_i,s_vol);
		const struct const_Matrix_T*const op__t_to_s =
			constructor_operator__test_s_coef_to_sol_coef_T(flux_r,opg_s_vol); // destructed
		destructor_Flux_Ref_T(flux_r);

		struct Vector_T test_s_coef_V = interpret_Multiarray_as_Vector_T(s_vol->test_s_coef);
		struct Vector_T s_coef_V      = interpret_Multiarray_as_Vector_T(s_vol->sol_coef);
		mv_T('N',1.0,1.0,op__t_to_s,(struct const_Vector_T*)&test_s_coef_V,&s_coef_V);

		destructor_const_Matrix_T(op__t_to_s);
	}
	destructor_Flux_Input_T(flux_i);

	test_case->solver_method_curr = smc;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct S_Params_T set_s_params_T (const struct Simulation*const sim)
{
	struct S_Params_T s_params;

	set_S_Params_Volume_Structor_T(&s_params.spvs,sim);

	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	switch (test_case->solver_method_curr) {
#if TYPE_RC == TYPE_COMPLEX
	case 'e':
		s_params.compute_rlhs = compute_rhs_v_dg_like_T;
		break;
#elif TYPE_RC == TYPE_REAL
	case 'i':
		if (test_case->has_1st_order && !test_case->has_2nd_order)
			s_params.compute_rlhs = compute_rlhs_1;
		else if (!test_case->has_1st_order && test_case->has_2nd_order)
			EXIT_ADD_SUPPORT; //s_params.compute_rlhs = compute_rlhs_2;
		else if (test_case->has_1st_order && test_case->has_2nd_order)
			EXIT_ADD_SUPPORT; //s_params.compute_rlhs = compute_rlhs_12;
		else
			EXIT_ADD_SUPPORT;
		break;
#endif
	default:
		EXIT_ERROR("Unsupported: %c (type_rc: %d)\n",test_case->solver_method_curr,TYPE_RC);
		break;
	}

	return s_params;
}

static const struct const_Matrix_T* constructor_test_diff_op_1v_opg_T
	(const struct Flux_Ref_T*const flux_r, const struct Solver_Volume_T*const s_vol, const bool include_det_j)
{
	const int*const n_var_eq = get_set_n_var_eq(NULL);
	const int n_vr = n_var_eq[0],
	          n_eq = n_var_eq[1];

	const struct Multiarray_Operator cv1_vt_vc = get_operator__cv1_vt_vc_T(s_vol);

	const ptrdiff_t ext_0 = cv1_vt_vc.data[0]->op_std->ext_0,
	                ext_1 = cv1_vt_vc.data[0]->op_std->ext_1;

	struct Matrix_T* cv1r = constructor_empty_Matrix_T('R',n_vr*ext_0,n_eq*ext_1); // returned

	struct Matrix_T* cv1r_l = constructor_empty_Matrix_T('R',ext_0,ext_1); // destructed
	const struct const_Multiarray_T* dfr_ds_Ma = flux_r->dfr_ds;

	assert(dfr_ds_Ma != NULL);
	struct Vector_T dfr_ds = { .ext_0 = dfr_ds_Ma->extents[0], .owns_data = false, .data = NULL, };

	for (int vr = 0; vr < n_vr; ++vr) {
	for (int eq = 0; eq < n_eq; ++eq) {
		set_to_value_Matrix_T(cv1r_l,0.0);
		for (int dim = 0; dim < DIM; ++dim) {
			const ptrdiff_t ind =
				compute_index_sub_container(dfr_ds_Ma->order,1,dfr_ds_Ma->extents,(ptrdiff_t[]){eq,vr,dim});
			dfr_ds.data = (Type*)&dfr_ds_Ma->data[ind];
			mm_diag_T('L',-1.0,1.0,cv1_vt_vc.data[dim]->op_std,(struct const_Vector_T*)&dfr_ds,cv1r_l,false);
		}
		/** \note Building the right operator here (as opposed to the left in \ref constructor_lhs_v_1_T). The
		 *  flux Jacobian terms are thus transposed by swapping the 'eq' and 'vr' indices. */
		set_block_Matrix_T(cv1r,vr*ext_0,eq*ext_1,
		                   (struct const_Matrix_T*)cv1r_l,0,0,cv1r_l->ext_0,cv1r_l->ext_1,'i');
	}}
	destructor_Matrix_T(cv1r_l);

	if (include_det_j) {
		const struct const_Vector_T J_vc = interpret_const_Multiarray_as_Vector_T(s_vol->jacobian_det_vc);
		const struct const_Vector_T*const Jr_vc = constructor_repeated_const_Vector_T(1.0,&J_vc,n_vr); // dest.
		scale_Matrix_by_Vector_T('L',1.0,cv1r,Jr_vc,true);
		destructor_const_Vector_T(Jr_vc);
	}
	return (struct const_Matrix_T*) cv1r;
}

#include "undef_templates_compute_volume_rlhs_opg.h"

#include "undef_templates_matrix.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"

#include "undef_templates_volume_solver_opg.h"

#include "undef_templates_compute_volume_rlhs.h"
#include "undef_templates_flux.h"
#include "undef_templates_test_case.h"
#include "undef_templates_operators.h"
