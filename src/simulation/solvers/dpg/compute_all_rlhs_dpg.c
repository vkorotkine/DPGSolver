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

#include <assert.h>
#include <stdlib.h>

#include "macros.h"
#include "definitions_dpg.h"
#include "definitions_elements.h"
#include "definitions_intrusive.h"

#include "face_solver_dpg.h"
#include "volume_solver_dpg.h"
#include "element_solver_dpg.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "compute_volume_rlhs.h"
#include "flux.h"
#include "intrusive.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Function pointer to functions constructing the norm operator used to evaluate the optimal test functions.
 *
 *  \param dpg_s_vol The current volume.
 *  \param flux_r    \ref Flux_Ref.
 *  \param sim       \ref Simulation.
 */
typedef const struct const_Matrix_d* (*constructor_norm_op_fptr)
	(const struct DPG_Solver_Volume* dpg_s_vol,
	 const struct Flux_Ref* flux_r,
	 const struct Simulation* sim
	);

/** \brief Pointer to functions computing rhs and lhs terms for the dpg solver.
 *
 *  \param norm_op   The dpg norm operator.
 *  \param flux_r    \ref Flux_Ref.
 *  \param dpg_s_vol Pointer to the current volume.
 *  \param sim       \ref Simulation.
 */
typedef void (*compute_rlhs_fptr)
	(const struct const_Matrix_d* norm_op,
	 const struct Flux_Ref* flux_r,
	 struct DPG_Solver_Volume* dpg_s_vol,
	 const struct Simulation* sim
	);

/// \brief Container for solver-related parameters.
struct S_Params {
	struct S_Params_Volume_Structor spvs; ///< \ref S_Params_Volume_Structor.

	constructor_norm_op_fptr constructor_norm_op; ///< Pointer to the appropriate function.
	compute_rlhs_fptr compute_rlhs;               ///< Pointer to the appropriate function.
};

/** \brief Set the parameters of \ref S_Params.
 *  \return A statically allocated \ref S_Params container. */
static struct S_Params set_s_params
	(const struct Simulation* sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void compute_all_rlhs_dpg (const struct Simulation* sim, struct Solver_Storage_Implicit* ssi)
{
	assert(sim->test_case->solver_method_curr == 'i');
	assert(sim->volumes->name == IL_VOLUME_SOLVER_DPG);
	assert(sim->faces->name == IL_FACE_SOLVER_DPG);
	assert(sim->elements->name == IL_ELEMENT_SOLVER_DPG);

	struct S_Params s_params = set_s_params(sim);
	struct Flux_Input* flux_i = constructor_Flux_Input(sim); // destructed

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
//struct Volume*        vol   = (struct Volume*) curr;
//printf("v_ind: %d\n",vol->index);
		struct Solver_Volume* s_vol         = (struct Solver_Volume*) curr;
		struct DPG_Solver_Volume* dpg_s_vol = (struct DPG_Solver_Volume*) curr;

		struct Flux_Ref* flux_r = constructor_Flux_Ref_vol(&s_params.spvs,flux_i,s_vol,sim);

		const struct const_Matrix_d* norm_op = s_params.constructor_norm_op(dpg_s_vol,flux_r,sim); // destructed

		s_params.compute_rlhs(norm_op,flux_r,dpg_s_vol,sim);
		destructor_const_Matrix_d(norm_op);
		destructor_Flux_Ref(flux_r);
EXIT_UNSUPPORTED;

UNUSED(ssi);
EXIT_UNSUPPORTED;
	}
	destructor_Flux_Input(flux_i);
}

struct Multiarray_Operator get_operator__cvt1_vt_vc__rlhs (const struct DPG_Solver_Volume* dpg_s_vol)
{
	struct Volume* vol          = (struct Volume*) dpg_s_vol;
	struct Solver_Volume* s_vol = (struct Solver_Volume*) vol;

	const struct DPG_Solver_Element* dpg_s_e = (struct DPG_Solver_Element*) vol->element;

	const int p      = s_vol->p_ref,
	          curved = vol->curved;

	return set_MO_from_MO(dpg_s_e->cvt1_vt_vc[curved],1,(ptrdiff_t[]){0,0,p,p});
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Version of \ref constructor_norm_op_fptr; see comments for \ref TEST_NORM_H1_UPWIND.
 *  \return See brief. */
static const struct const_Matrix_d* constructor_norm_op__h1_upwind
	(const struct DPG_Solver_Volume* dpg_s_vol, ///< See brief.
	 const struct Flux_Ref* flux_r,             ///< See brief.
	 const struct Simulation* sim               ///< See brief.
	);

/// \brief Version of \ref compute_rlhs_fptr computing rhs and lhs terms for 1st order equations only.
static void compute_rlhs_1
	(const struct const_Matrix_d* norm_op, ///< See brief.
	 const struct Flux_Ref* flux_r,        ///< See brief.
	 struct DPG_Solver_Volume* dpg_s_vol,  ///< See brief.
	 const struct Simulation* sim          ///< See brief.
	);

static struct S_Params set_s_params (const struct Simulation* sim)
{
	struct S_Params s_params;

	set_S_Params_Volume_Structor(&s_params.spvs,sim);

	struct Test_Case* test_case = sim->test_case;
	switch (test_case->solver_method_curr) {
	case 'i':
		if (test_case->has_1st_order && !test_case->has_2nd_order)
			s_params.compute_rlhs = compute_rlhs_1;
		else if (!test_case->has_1st_order && test_case->has_2nd_order)
			EXIT_ADD_SUPPORT; //s_params.compute_rlhs = compute_rlhs_2;
		else if (test_case->has_1st_order && test_case->has_2nd_order)
			EXIT_ADD_SUPPORT; //s_params.compute_rlhs = compute_rlhs_12;
		else
			EXIT_ERROR("Unsupported: %d %d\n",test_case->has_1st_order,test_case->has_2nd_order);
		break;
	case 'e': // fallthrough
	default:
		EXIT_ERROR("Unsupported: %c\n",test_case->solver_method_curr);
		break;
	}

	switch (test_case->ind_test_norm) {
		case TEST_NORM_H1_UPWIND: s_params.constructor_norm_op = constructor_norm_op__h1_upwind; break;
		default:                  EXIT_ERROR("Unsupported: %d\n",test_case->ind_test_norm);      break;
	}

	return s_params;
}

// Level 1 ********************************************************************************************************** //

/** \brief Constructor for the rhs \ref Vector_d with volume contributions from 1st order equations included.
 *  \return See brief. */
static struct Vector_d* constructor_rhs_vol_1
	(const struct Flux_Ref* flux_r,     ///< Defined for \ref compute_rlhs_fptr.
	 const struct Solver_Volume* s_vol, ///< Defined for \ref compute_rlhs_fptr.
	 const struct Simulation* sim       ///< Defined for \ref compute_rlhs_fptr.
	);

static const struct const_Matrix_d* constructor_norm_op__h1_upwind
	(const struct DPG_Solver_Volume* dpg_s_vol, const struct Flux_Ref* flux_r, const struct Simulation* sim)
{
	const int d    = sim->d,
	          n_eq = sim->test_case->n_eq,
	          n_vr = sim->test_case->n_var;

	const struct Multiarray_Operator cvt1_vt_vc = get_operator__cvt1_vt_vc__rlhs(dpg_s_vol);

	const ptrdiff_t ext_0 = cvt1_vt_vc.data[0]->op_std->ext_0,
	                ext_1 = cvt1_vt_vc.data[0]->op_std->ext_1;

	struct Matrix_d* cvt1r = constructor_empty_Matrix_d('R',n_eq*ext_0,n_vr*ext_1); // destructed

	struct Matrix_d* cvt1r_l = constructor_empty_Matrix_d('R',ext_0,ext_1); // destructed
	const struct const_Multiarray_d* dfr_ds_Ma = flux_r->dfr_ds;
	struct Vector_d dfr_ds = { .ext_0 = dfr_ds_Ma->extents[0], .owns_data = false, .data = NULL, };

	for (int vr = 0; vr < n_vr; ++vr) {
	for (int eq = 0; eq < n_eq; ++eq) {
		set_to_value_Matrix_d(cvt1r_l,0.0);
		for (int dim = 0; dim < d; ++dim) {
			const ptrdiff_t ind =
				compute_index_sub_container(dfr_ds_Ma->order,1,dfr_ds_Ma->extents,(ptrdiff_t[]){eq,vr,dim});
			dfr_ds.data = (double*)&dfr_ds_Ma->data[ind];
			mm_diag_d('R',1.0,1.0,cvt1_vt_vc.data[dim]->op_std,(struct const_Vector_d*)&dfr_ds,cvt1r_l,false);
		}
		set_block_Matrix_d(cvt1r,(struct const_Matrix_d*)cvt1r_l,eq*ext_0,vr*ext_1,'i');
	}}
	destructor_Matrix_d(cvt1r_l);

	struct Solver_Volume* s_vol = (struct Solver_Volume*) dpg_s_vol;
	const struct const_Vector_d* w_vc = get_operator__w_vc__s_e(s_vol);
	const struct const_Vector_d J_vc  = interpret_const_Multiarray_as_Vector_d(s_vol->jacobian_det_vc);

	const struct const_Vector_d* J_inv_vc = constructor_inverse_const_Vector_d(&J_vc);               // destructed
	const struct const_Vector_d* wJ_vc    = constructor_dot_mult_const_Vector_d(w_vc,J_inv_vc,n_vr); // destructed
	destructor_const_Vector_d(J_inv_vc);

	const struct const_Matrix_d* n1_l =
		constructor_mm_diag_const_Matrix_d(1.0,(struct const_Matrix_d*)cvt1r,wJ_vc,'R',false); // destructed
	destructor_const_Vector_d(wJ_vc);

	const struct const_Matrix_d* n1 =
		constructor_mm_const_Matrix_d('N','T',1.0,n1_l,(struct const_Matrix_d*)cvt1r,'R'); // destructed
	destructor_const_Matrix_d(n1_l);

	const struct const_Matrix_d* norm_op_H0 = dpg_s_vol->norm_op_H0;
	assert(norm_op_H0->ext_0 == ext_0);

	struct Matrix_d* norm_op = constructor_empty_Matrix_d('R',n_eq*ext_0,n_eq*ext_0); // returned

	set_block_Matrix_d(norm_op,n1,0,0,'i');
	for (int eq = 0; eq < n_eq; ++eq)
		set_block_Matrix_d(norm_op,norm_op_H0,eq*ext_0,eq*ext_0,'a');
	destructor_const_Matrix_d(n1);

	return (struct const_Matrix_d*) norm_op;
}

static void compute_rlhs_1
	(const struct const_Matrix_d* norm_op, const struct Flux_Ref* flux_r, struct DPG_Solver_Volume* dpg_s_vol,
	 const struct Simulation* sim)
{
// Likely:
// 2) Add to rhs with face terms (increment dof).
// 3) Constructor lhs using dof.
// 4) Add to lhs volume.
// 5) Add to lhs face.
	const struct Solver_Volume* s_vol = (struct Solver_Volume*) dpg_s_vol;

	struct Vector_d* rhs = constructor_rhs_vol_1(flux_r,s_vol,sim); // destructed

	// rhs (face)
// loop over volume faces. Check if volume is 'l' or 'r', appropriately reorder the operator.

	destructor_Vector_d(rhs);
UNUSED(norm_op);
EXIT_UNSUPPORTED;
}

// Level 2 ********************************************************************************************************** //

static struct Vector_d* constructor_rhs_vol_1
	(const struct Flux_Ref* flux_r, const struct Solver_Volume* s_vol, const struct Simulation* sim)
{
	const int d    = sim->d,
	          n_eq = sim->test_case->n_eq;

	const struct Multiarray_Operator tw1_vt_vc = get_operator__tw1_vt_vc__rlhs(s_vol);

	const ptrdiff_t ext_0 = tw1_vt_vc.data[0]->op_std->ext_0;

	struct Vector_d* rhs = constructor_zero_Vector_d(ext_0*n_eq); // returned

	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	const char op_format = 'd';

	ptrdiff_t extents[2] = { ext_0, n_eq, };
	struct Multiarray_d rhs_Ma =
		{ .layout = 'C', .order = 2, .extents = extents, .owns_data = false, .data = rhs->data, };
	for (ptrdiff_t dim = 0; dim < d; ++dim)
		mm_NNC_Operator_Multiarray_d(1.0,1.0,tw1_vt_vc.data[dim],flux_r->fr,&rhs_Ma,op_format,2,&dim,NULL);

	return rhs;
}

static void increment_rhs_face_1 (struct Vector_d* rhs, const struct Solver_Volume* s_vol, const struct Simulation* sim)
{
	const struct Volume* vol = (struct Volume*) s_vol;

	const int n_eq = sim->test_case->n_eq;
	const ptrdiff_t ext_0 = (rhs->ext_0)/n_eq;

	struct Matrix_d rhs_M = { .layout = 'C', .ext_0 = ext_0, .ext_1 = n_eq, .owns_data = false, .data = rhs->data, };
	for (int i = 0; i < NFMAX;    ++i) {
	for (int j = 0; j < NSUBFMAX; ++j) {
		struct Face* vol->face[i][j];
		if (!face)
			continue;

		const int side_index = compute_side_index_face(face,vol);
// get_operator_* needs name change.
		const struct Operator* tw0_vt_fc_op = get_operator__tw0_vt_fc__rlhs_dg(side_index,face),
		                     * cv0_ff_fc_op = get_operator__cv0_ff_fc__dpg_s_e(side_index,dpg_s_face);

		struct Solver_Face* s_face = (struct Solver_Face*) face;

		struct Matrix_d* tw0_vt_fc =
			constructor_copy_Matrix_d((struct Matrix_d*)tw0_vt_fc_op->op_std); // destructed.

		const struct const_Vector_d j_det_V = interpret_const_Multiarray_as_Vector_d(s_face->jacobian_det_fc);
		scale_Matrix_by_Vector_d(side,1.0,a,&j_det_V,false);
EXIT_UNSUPPORTED;
		if (side_index == 0) {
			mm_NNC_Operator_Multiarray_d(-1.0,1.0,tw0_vt_fc,nf_coef,&rhs_Ma,op_format,2,NULL,NULL);
		} else {
			permute_Multiarray_d_fc(nf_coef,'R',1,face);
//scale_Multiarray_d((struct Multiarray_d*)num_flux->nnf,-1.0); // Use "-ve" normal.
// Can remove both of the scalings (here and and finalize_face_rhs_dg).
			mm_NNC_Operator_Multiarray_d(-1.0,1.0,tw0_vt_fc,nf_coef,&rhs_Ma,op_format,2,NULL,NULL);
			EXIT_ADD_SUPPORT;
		}
		struct Matrix_d nf_coef = interpret_Multiarray_as_Matrix_d(s_face->nf_coef);

// Remove redundant "-ve" signs and add comment.
		mm_d('N','N',-1.0,1.0,(struct const_Matrix_d*)lhs_l,(struct const_Matrix_d*)&nf_coef,&rhs_M);
	}}
// Note that the lhs is **always** linear wrt the trace unknowns and the lhs terms should thus be computed at the same
// time as the rhs.
EXIT_UNSUPPORTED;
}

// Level 3 ********************************************************************************************************** //

const struct Operator* get_operator__cv0_ff_fc__dpg_s_e (const int side_index, const struct DPG_Solver_Face* dpg_s_face)
{
	const struct Face* face           = (struct Face*) dpg_s_face;
	const struct Solver_Face* s_face  = (struct Solver_Face*) dpg_s_face;
	const struct Volume* vol          = face->neigh_info[side_index].volume;
	const struct Solver_Volume* s_vol = (struct Solver_Volume*) vol;

	const struct DPG_Solver_Element* dpg_s_e = (struct DPG_Solver_Element*) vol->element;

	const int ind_lf = face->neigh_info[side_index].ind_lf,
	          ind_e  = get_face_element_index(face);
	const int p_f = s_face->p_ref;

	const int curved = ( (s_face->cub_type == 's') ? 0 : 1 );
	return get_Multiarray_Operator(dpg_s_e->cv0_ff_fc[curved],(ptrdiff_t[]){ind_lf,ind_e,ind_e,0,0,p_f,p_f});
}

