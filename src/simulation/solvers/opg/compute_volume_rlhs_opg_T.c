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

#include "def_templates_compute_rlhs.h"
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

/** \brief Constructor for the solution at the 'v'olume 's'olution nodes.
 *  \return See brief. */
static const struct const_Multiarray_T* constructor_sol_vs_T
	(const struct Solver_Volume_T*const s_vol ///< Standard.
	 );

/** \brief Constructor for the gradients at the 'v'olume 's'olution nodes.
 *  \return See brief. */
static const struct const_Multiarray_T* constructor_grad_vs_T
	(const struct Solver_Volume_T*const s_vol ///< Standard.
	 );

/** \brief Constructor for the xyz coordinates at the 'v'olume 's'olution nodes.
 *  \return See brief. */
static const struct const_Multiarray_T* constructor_xyz_vs_T
	(const struct Solver_Volume_T*const s_vol ///< Standard.
	 );

/** \brief Constructor for the 'test' function 'diff'erentiation 'op'erator for the '1'st order 'v'olume term for the
 *         OPG scheme.
 *  \return See brief.
 *
 *  This operator is the rhs portion of the operator returned from \ref constructor_operator__test_s_coef_to_sol_coef_T.
 */
static const struct const_Matrix_T* constructor_test_diff_op_1v_opg_T
	(const struct Flux_Ref_T*const flux_r,            ///< Standard.
	 const struct OPG_Solver_Volume_T*const opg_s_vol ///< Standard.
	 );

/// \brief Version of \ref compute_rlhs_v_fptr_T computing only the rhs term for the opg scheme.
static void compute_rhs_v_opg_T
	(const struct Flux_Ref_T*const flux_r,    ///< See brief.
	 struct Solver_Volume_T*const s_vol,      ///< See brief.
	 struct Solver_Storage_Implicit*const ssi ///< See brief.
		);

// Interface functions ********************************************************************************************** //

void compute_volume_rlhs_opg_T
	(const struct Simulation* sim, struct Solver_Storage_Implicit* ssi, struct Intrusive_List* volumes)
{
	assert(sim->volumes->name == IL_VOLUME_SOLVER_OPG);
	assert(sim->elements->name == IL_ELEMENT_SOLVER_OPG);

	struct S_Params_T s_params = set_s_params_T(sim);

	struct Test_Case_T*const test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const char smc = test_case->solver_method_curr;
	test_case->solver_method_curr = 'i'; // Jacobians needed for rhs volume term.
	struct Flux_Input_T* flux_i = constructor_Flux_Input_T(sim); // destructed
	test_case->solver_method_curr = smc;

	for (struct Intrusive_Link* curr = volumes->first; curr; curr = curr->next) {
		struct Solver_Volume_T*const s_vol = (struct Solver_Volume_T*) curr;

		struct Flux_Ref_T* flux_r = constructor_Flux_Ref_vol_opg_T(flux_i,s_vol); // destructed

		// Compute the rhs and the lhs terms.
		s_params.compute_rlhs(flux_r,s_vol,ssi);
		destructor_Flux_Ref_T(flux_r);
	}
	destructor_Flux_Input_T(flux_i);
}

struct Flux_Ref_T* constructor_Flux_Ref_vol_opg_T (struct Flux_Input_T* flux_i, const struct Solver_Volume_T* s_vol)
{
	flux_i->s   = constructor_sol_vs_T(s_vol);
	flux_i->g   = constructor_grad_vs_T(s_vol);
	flux_i->xyz = constructor_xyz_vs_T(s_vol);

	struct Flux_T* flux = constructor_Flux_T(flux_i); // destructed
	destructor_const_Multiarray_T(flux_i->s);
	destructor_conditional_const_Multiarray_T(flux_i->g);
	destructor_conditional_const_Multiarray_T(flux_i->xyz);

	struct Flux_Ref_T* flux_r = constructor_Flux_Ref_T(s_vol->metrics_vs,flux);
	destructor_Flux_T(flux);

	return flux_r;
}

const struct const_Matrix_T* constructor_operator__test_s_coef_to_sol_coef_T
	(const struct Flux_Ref_T*const flux_r, const struct OPG_Solver_Volume_T*const opg_s_vol)
{
	const struct const_Matrix_R*const vc0_vs_vs = get_operator__vc0_vs_vs_T(opg_s_vol)->op_std;
	const struct const_Matrix_T*const op_v1_opg = constructor_test_diff_op_1v_opg_T(flux_r,opg_s_vol); // destructed

	const int n_eq = get_set_n_var_eq(NULL)[1];
	const char layout = 'R';
	const ptrdiff_t sub_ext = vc0_vs_vs->ext_0;

	struct Matrix_T*const op_final = constructor_empty_Matrix_T(layout,n_eq*vc0_vs_vs->ext_0,op_v1_opg->ext_1); // r.
	struct Matrix_T sub_op_v1_opg = { .layout = layout, .ext_0 = sub_ext, .ext_1 = op_v1_opg->ext_1, .data = NULL, };
	struct Matrix_T sub_op_final  = { .layout = layout, .ext_0 = sub_ext, .ext_1 = op_v1_opg->ext_1, .data = NULL, };

	for (int eq = 0; eq < n_eq; ++eq) {
		sub_op_v1_opg.data = (Type*) get_row_const_Matrix_T(sub_ext*eq,op_v1_opg);
		sub_op_final.data  = get_row_Matrix_T(sub_ext*eq,op_final);
		mm_RTT('N','N',1.0,0.0,vc0_vs_vs,(struct const_Matrix_T*)&sub_op_v1_opg,&sub_op_final);
	}
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

	for (struct Intrusive_Link* curr = volumes->first; curr; curr = curr->next) {
		struct Solver_Volume_T*const s_vol         = (struct Solver_Volume_T*) curr;
		struct OPG_Solver_Volume_T*const opg_s_vol = (struct OPG_Solver_Volume_T*) curr;

		struct Flux_Ref_T*const flux_r = constructor_Flux_Ref_vol_opg_T(flux_i,s_vol);
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
		s_params.compute_rlhs = compute_rhs_v_opg_T;
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

static const struct const_Multiarray_T* constructor_sol_vs_T (const struct Solver_Volume_T*const s_vol)
{
	const struct Operator*const cv0_vs_vs = get_operator__cv0_vs_vs_T(s_vol);
	const struct const_Multiarray_T*const s_coef = (const struct const_Multiarray_T*) s_vol->sol_coef;
	const char op_format = get_set_op_format(0);
	return constructor_mm_NN1_Operator_const_Multiarray_T(cv0_vs_vs,s_coef,'C',op_format,s_coef->order,NULL);
}

static const struct const_Multiarray_T* constructor_grad_vs_T (const struct Solver_Volume_T*const s_vol)
{
	if (!get_set_has_1st_2nd_order(NULL)[1])
		return NULL;

	const struct Operator*const cv0_vr_vs = get_operator__cv0_vr_vs_T(s_vol);
	const struct const_Multiarray_T*const g_coef = (const struct const_Multiarray_T*) s_vol->grad_coef;
	const char op_format = get_set_op_format(0);
	return constructor_mm_NN1_Operator_const_Multiarray_T(cv0_vr_vs,g_coef,'C',op_format,g_coef->order,NULL);
}

static const struct const_Multiarray_T* constructor_xyz_vs_T (const struct Solver_Volume_T*const s_vol)
{
	switch (get_set_pde_index(NULL)) {
	case PDE_ADVECTION:
		break; // Do nothing (continue below).
	case PDE_DIFFUSION: case PDE_EULER: case PDE_NAVIER_STOKES:
		return NULL;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",get_set_pde_index(NULL));
		break;
	}

	const struct Operator*const cv0_vg_vs = get_operator__cv0_vg_vs_T(s_vol);
	const struct const_Multiarray_T*const g_coef = s_vol->geom_coef;
	const char op_format = get_set_op_format(0);
	return constructor_mm_NN1_Operator_const_Multiarray_T(cv0_vg_vs,g_coef,'C',op_format,g_coef->order,NULL);
}

static const struct const_Matrix_T* constructor_test_diff_op_1v_opg_T
	(const struct Flux_Ref_T*const flux_r, const struct OPG_Solver_Volume_T*const opg_s_vol)
{
	const int*const n_var_eq = get_set_n_var_eq(NULL);
	const int n_vr = n_var_eq[0],
	          n_eq = n_var_eq[1];

	const struct Multiarray_Operator cv1_vt_vs = get_operator__cv1_vt_vs_T(opg_s_vol);

	const ptrdiff_t ext_0 = cv1_vt_vs.data[0]->op_std->ext_0,
	                ext_1 = cv1_vt_vs.data[0]->op_std->ext_1;

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
			mm_diag_T('L',-1.0,1.0,cv1_vt_vs.data[dim]->op_std,(struct const_Vector_T*)&dfr_ds,cv1r_l,false);
		}
		/** \note Building the right operator here (as opposed to the left in \ref constructor_lhs_v_1_T). The
		 *  flux Jacobian terms are thus transposed by swapping the 'eq' and 'vr' indices. */
		set_block_Matrix_T(cv1r,vr*ext_0,eq*ext_1,
		                   (struct const_Matrix_T*)cv1r_l,0,0,cv1r_l->ext_0,cv1r_l->ext_1,'i');
	}}
	destructor_Matrix_T(cv1r_l);

	const struct Solver_Volume_T*const s_vol = (struct Solver_Volume_T*) opg_s_vol;
	const struct const_Vector_T J_vs = interpret_const_Multiarray_as_Vector_T(s_vol->jacobian_det_vs);
	const struct const_Vector_T*const Jr_vs = constructor_repeated_const_Vector_T(1.0,&J_vs,n_vr); // dest.
	scale_Matrix_by_Vector_T('L',1.0,cv1r,Jr_vs,true);
	destructor_const_Vector_T(Jr_vs);

	return (struct const_Matrix_T*) cv1r;
}

static void compute_rhs_v_opg_T
	(const struct Flux_Ref_T*const flux_r, struct Solver_Volume_T*const s_vol,
	 struct Solver_Storage_Implicit*const ssi)
{
	UNUSED(ssi);
	const int n_vr = get_set_n_var_eq(NULL)[0];

	const struct OPG_Solver_Volume_T*const opg_s_vol = (struct OPG_Solver_Volume_T*) s_vol;
	const struct const_Matrix_T*const op_t_to_s =
		constructor_operator__test_s_coef_to_sol_coef_T(flux_r,opg_s_vol); // destructed
	const struct const_Matrix_T*const M = constructor_block_diagonal_const_Matrix_T(opg_s_vol->m,n_vr); // dest.
	const struct const_Matrix_T*const lhs_l = constructor_mm_const_Matrix_T('T','N',1.0,op_t_to_s,M,'R'); // dest.
	destructor_const_Matrix_T(op_t_to_s);
	destructor_const_Matrix_T(M);

	// -ve sign included to cancel that coming from op_t_to_s for the v'*df_ds term.
	mm_NNC_Multiarray_TTT(-1.0,1.0,lhs_l,(struct const_Multiarray_T*)s_vol->sol_coef,s_vol->rhs);
	destructor_const_Matrix_T(lhs_l);
}

#include "undef_templates_compute_volume_rlhs_opg.h"

#include "undef_templates_matrix.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"

#include "undef_templates_volume_solver_opg.h"

#include "undef_templates_compute_rlhs.h"
#include "undef_templates_compute_volume_rlhs.h"
#include "undef_templates_flux.h"
#include "undef_templates_test_case.h"
#include "undef_templates_operators.h"
