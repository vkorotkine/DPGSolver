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
#include "petscmat.h"

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
 *  \param ssi       \ref Solver_Storage_Implicit.
 *  \param sim       \ref Simulation.
 */
typedef void (*compute_rlhs_fptr)
	(const struct const_Matrix_d* norm_op,
	 const struct Flux_Ref* flux_r,
	 const struct DPG_Solver_Volume* dpg_s_vol,
	 struct Solver_Storage_Implicit* ssi,
	 const struct Simulation* sim
	);

/// \brief Container for solver-related parameters.
struct S_Params_DPG {
	struct S_Params_Volume_Structor spvs; ///< \ref S_Params_Volume_Structor.

	constructor_norm_op_fptr constructor_norm_op; ///< Pointer to the appropriate function.
	compute_rlhs_fptr compute_rlhs;               ///< Pointer to the appropriate function.
};

/** \brief Set the parameters of \ref S_Params_DPG.
 *  \return A statically allocated \ref S_Params_DPG container. */
static struct S_Params_DPG set_s_params_dpg
	(const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Get the pointer to the appropriate \ref DPG_Solver_Element::cv0_ff_fc operator.
 *  \return See brief. */
static const struct Operator* get_operator__cv0_ff_fc
	(const int side_index,                    ///< The index of the side of the face under consideration.
	 const struct DPG_Solver_Face* dpg_s_face ///< The current face.
	);

/// \brief Set the values of the global indices corresponding to the current unknown.
static void set_idxm
	(int* ind_idxm,                  ///< Pointer to the index in `idxm`.
	 struct Vector_i* idxm,          ///< Vector of global indices.
	 const int ind_dof,              ///< Index of the first dof corresponding to the `coef`.
	 const struct Multiarray_d* coef ///< Pointer to the coefficients under consideration.
	);

// Interface functions ********************************************************************************************** //

void compute_all_rlhs_dpg (const struct Simulation* sim, struct Solver_Storage_Implicit* ssi)
{
	assert(sim->test_case->solver_method_curr == 'i');
	assert(sim->volumes->name == IL_VOLUME_SOLVER_DPG);
	assert(sim->faces->name == IL_FACE_SOLVER_DPG);
	assert(sim->elements->name == IL_ELEMENT_SOLVER_DPG);

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

struct Multiarray_Operator get_operator__cvt1_vt_vc__rlhs (const struct DPG_Solver_Volume* dpg_s_vol)
{
	struct Volume* vol          = (struct Volume*) dpg_s_vol;
	struct Solver_Volume* s_vol = (struct Solver_Volume*) vol;

	const struct DPG_Solver_Element* dpg_s_e = (struct DPG_Solver_Element*) vol->element;

	const int p      = s_vol->p_ref,
	          curved = vol->curved;

	return set_MO_from_MO(dpg_s_e->cvt1_vt_vc[curved],1,(ptrdiff_t[]){0,0,p,p});
}

const struct const_Matrix_d* constructor_lhs_l_internal_face_dpg
	(const struct DPG_Solver_Volume* dpg_s_vol, const struct DPG_Solver_Face* dpg_s_face)
{
	const struct Volume* vol         = (struct Volume*) dpg_s_vol;
	const struct Face* face          = (struct Face*) dpg_s_face;
	const struct Solver_Face* s_face = (struct Solver_Face*) dpg_s_face;

	const int side_index = compute_side_index_face(face,vol);
//printf("%d %d %d\n",vol->index,face->index,side_index);
	const struct Operator* tw0_vt_fc_op = get_operator__tw0_vt_fc(side_index,s_face),
	                     * cv0_ff_fc_op = get_operator__cv0_ff_fc(side_index,dpg_s_face);

	struct Matrix_d* cv0_ff_fc = constructor_copy_Matrix_d((struct Matrix_d*)cv0_ff_fc_op->op_std); // destructed

	const struct const_Multiarray_d* j_det = s_face->jacobian_det_fc;
	const struct const_Vector_d* j_det_V =
		constructor_copy_const_Vector_d_d(compute_size(j_det->order,j_det->extents),j_det->data); // destructed
	scale_Matrix_by_Vector_d('L',1.0,cv0_ff_fc,j_det_V,false);
	destructor_const_Vector_d(j_det_V);

	// Use "-ve" sign when looking from volume[0] as the face term was moved to the rhs of the equation. When looking
	// from volume[1], the "-ve" sign is cancelled by the "-ve" sign of the inverted normal vector.
	double alpha = -1.0;
	if (side_index == 1) {
		permute_Matrix_d_fc(cv0_ff_fc,'R',side_index,s_face);
		alpha = 1.0;
	}

	const struct const_Matrix_d* lhs_l = constructor_mm_const_Matrix_d
		('N','N',alpha,tw0_vt_fc_op->op_std,(struct const_Matrix_d*)cv0_ff_fc,'R'); // returned
	destructor_Matrix_d(cv0_ff_fc);

	return lhs_l;
}

ptrdiff_t compute_n_dof_nf (const struct Solver_Volume* s_vol)
{
	ptrdiff_t dof = 0;

	const struct Volume* vol = (struct Volume*) s_vol;
	for (int i = 0; i < NFMAX;    ++i) {
	for (int j = 0; j < NSUBFMAX; ++j) {
		const struct Face* face = vol->faces[i][j];
		if (!face)
			continue;

		const struct Solver_Face* s_face = (struct Solver_Face*) face;
		const ptrdiff_t size = compute_size(s_face->nf_coef->order,s_face->nf_coef->extents);
		dof += size;

		assert((size > 0) || (face->boundary));
	}}
	return dof;
}

const struct const_Vector_i* constructor_petsc_idxm_dpg
	(const ptrdiff_t n_dof, const struct Solver_Volume* s_vol)
{
	struct Vector_i* idxm = constructor_empty_Vector_i(n_dof); // returned
	int ind_idxm = 0;

	// volume (sol_coef, grad_coef)
	set_idxm(&ind_idxm,idxm,s_vol->ind_dof,s_vol->sol_coef);
	// grad_coef: To be done.

	// face (nf_coef, sol_coef)
	const struct Volume* vol = (struct Volume*) s_vol;
	for (int i = 0; i < NFMAX;    ++i) {
	for (int j = 0; j < NSUBFMAX; ++j) {
		const struct Face* face = vol->faces[i][j];
		if (!face)
			continue;

		const struct Solver_Face* s_face = (struct Solver_Face*) face;
		set_idxm(&ind_idxm,idxm,s_face->ind_dof,s_face->nf_coef);
		// sol_coef: To be done.
	}}
printf("n_dof %td\n",n_dof);

	assert(ind_idxm == n_dof);

	return (struct const_Vector_i*) idxm;
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
	(const struct const_Matrix_d* norm_op,      ///< See brief.
	 const struct Flux_Ref* flux_r,             ///< See brief.
	 const struct DPG_Solver_Volume* dpg_s_vol, ///< See brief.
	 struct Solver_Storage_Implicit* ssi,       ///< See brief.
	 const struct Simulation* sim               ///< See brief.
	);

static struct S_Params_DPG set_s_params_dpg (const struct Simulation* sim)
{
	struct S_Params_DPG s_params;

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

static const struct Operator* get_operator__cv0_ff_fc (const int side_index, const struct DPG_Solver_Face* dpg_s_face)
{
	const struct Face* face            = (struct Face*) dpg_s_face;
	const struct Solver_Face* s_face   = (struct Solver_Face*) dpg_s_face;
	const struct Volume* vol           = face->neigh_info[side_index].volume;
	const struct DPG_Solver_Element* e = (struct DPG_Solver_Element*) vol->element;

	const int ind_e  = get_face_element_index(face),
	          p_f    = s_face->p_ref;
	const int curved = ( (s_face->cub_type == 's') ? 0 : 1 );

	return get_Multiarray_Operator(e->cv0_ff_fc[curved],(ptrdiff_t[]){ind_e,ind_e,0,0,p_f,p_f});
}

static void set_idxm (int* ind_idxm, struct Vector_i* idxm, const int ind_dof, const struct Multiarray_d* coef)
{
	if (coef == NULL)
		return;

	const ptrdiff_t size = compute_size(coef->order,coef->extents);
printf("size: %td\n",size);
	for (int i = 0; i < size; ++i) {
		idxm->data[*ind_idxm] = ind_dof+i;
		++*ind_idxm;
	}
}

// Level 1 ********************************************************************************************************** //

/** \brief Constructor for the rhs \ref Vector_d with volume contributions from 1st order equations included.
 *  \return See brief. */
static struct Vector_d* constructor_rhs_v_1
	(const struct Flux_Ref* flux_r,     ///< Defined for \ref compute_rlhs_fptr.
	 const struct Solver_Volume* s_vol, ///< Defined for \ref compute_rlhs_fptr.
	 const struct Simulation* sim       ///< Defined for \ref compute_rlhs_fptr.
	);

/// \brief Increment and add dof for the rhs and lhs with the face contributions from 1st order equations.
static void increment_and_add_dof_rlhs_f_1
	(struct Vector_d* rhs,                      ///< Holds the values of the rhs.
	 struct Matrix_d** lhs_ptr,                 ///< Pointer to the matrix holding the values of the lhs.
	 const struct DPG_Solver_Volume* dpg_s_vol, ///< \ref DPG_Solver_Volume.
	 const struct Simulation* sim               ///< \ref Simulation.
	);

/// \brief Increment the rhs terms with the source contribution.
static void increment_rhs_source
	(struct Vector_d* rhs,              ///< Holds the values of the rhs.
	 const struct Solver_Volume* s_vol, ///< Defined for \ref compute_rlhs_fptr.
	 const struct Simulation* sim       ///< Defined for \ref compute_rlhs_fptr.
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

static const struct const_Matrix_d* constructor_norm_op__h1_upwind
	(const struct DPG_Solver_Volume* dpg_s_vol, const struct Flux_Ref* flux_r, const struct Simulation* sim)
{
	const int n_eq = sim->test_case->n_eq,
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
		for (int dim = 0; dim < DIM; ++dim) {
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
	destructor_Matrix_d(cvt1r);

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
	(const struct const_Matrix_d* norm_op, const struct Flux_Ref* flux_r, const struct DPG_Solver_Volume* dpg_s_vol,
	 struct Solver_Storage_Implicit* ssi, const struct Simulation* sim)
{
	const struct Solver_Volume* s_vol = (struct Solver_Volume*) dpg_s_vol;

	struct Vector_d* rhs = constructor_rhs_v_1(flux_r,s_vol,sim); // destructed
	struct Matrix_d* lhs = constructor_lhs_v_1(flux_r,s_vol,sim); // destructed

	increment_and_add_dof_rlhs_f_1(rhs,&lhs,dpg_s_vol,sim);
	increment_rhs_source(rhs,s_vol,sim);

//print_const_Matrix_d(norm_op);
//print_Matrix_d(lhs);
//print_Multiarray_d(s_vol->sol_coef);
//print_Vector_d(rhs);
	const struct const_Matrix_d* optimal_test =
		constructor_sysv_const_Matrix_d(norm_op,(struct const_Matrix_d*)lhs); // destructed

//print_const_Matrix_d(optimal_test);

	const struct const_Matrix_d* lhs_opt =
		constructor_mm_const_Matrix_d('T','N',1.0,optimal_test,(struct const_Matrix_d*)lhs,'R'); // destructed
	destructor_Matrix_d(lhs);

	const struct const_Vector_d* rhs_opt =
		constructor_mv_const_Vector_d('T',-1.0,optimal_test,(struct const_Vector_d*)rhs); // destructed
	destructor_Vector_d(rhs);
	destructor_const_Matrix_d(optimal_test);

//print_const_Matrix_d(lhs_opt);
//print_const_Vector_d(rhs_opt);
	add_to_petsc_Mat_Vec_dpg(s_vol,rhs_opt,lhs_opt,ssi,sim);

	destructor_const_Matrix_d(lhs_opt);
	destructor_const_Vector_d(rhs_opt);
}

// Level 2 ********************************************************************************************************** //

/// \brief Increment the rhs and lhs entries corresponding to internal faces.
static void increment_rlhs_internal_face
	(const struct DPG_Solver_Volume* dpg_s_vol, ///< The current \ref DPG_Solver_Volume.
	 const struct DPG_Solver_Face* dpg_s_face,  ///< The current \ref DPG_Solver_Face.
	 struct Matrix_d* lhs,                      ///< The lhs matrix contribution for the current volume/faces.
	 struct Matrix_d* rhs,                      ///< The rhs matrix contribution for the current volume/faces.
	 int* ind_dof,                              ///< The index of the current dof under consideration.
	 const struct Simulation* sim               ///< \ref Simulation.
	);

/// \brief Increment the rhs and lhs entries corresponding to boundary faces.
static void increment_rlhs_boundary_face
	(const struct DPG_Solver_Volume* dpg_s_vol, ///< The current \ref DPG_Solver_Volume.
	 const struct DPG_Solver_Face* dpg_s_face,  ///< The current \ref DPG_Solver_Face.
	 struct Matrix_d* lhs,                      ///< The lhs matrix contribution for the current volume/faces.
	 struct Matrix_d* rhs,                      ///< The rhs matrix contribution for the current volume/faces.
	 const struct Simulation* sim               ///< \ref Simulation.
	);

static struct Vector_d* constructor_rhs_v_1
	(const struct Flux_Ref* flux_r, const struct Solver_Volume* s_vol, const struct Simulation* sim)
{
	const int n_eq = sim->test_case->n_eq;

	const struct Multiarray_Operator tw1_vt_vc = get_operator__tw1_vt_vc(s_vol);

	const ptrdiff_t ext_0 = tw1_vt_vc.data[0]->op_std->ext_0;

	struct Vector_d* rhs = constructor_zero_Vector_d(ext_0*n_eq); // returned

	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	const char op_format = 'd';

	ptrdiff_t extents[2] = { ext_0, n_eq, };
	struct Multiarray_d rhs_Ma =
		{ .layout = 'C', .order = 2, .extents = extents, .owns_data = false, .data = rhs->data, };
	for (ptrdiff_t dim = 0; dim < DIM; ++dim)
		mm_NNC_Operator_Multiarray_d(1.0,1.0,tw1_vt_vc.data[dim],flux_r->fr,&rhs_Ma,op_format,2,&dim,NULL);

	return rhs;
}

static void increment_and_add_dof_rlhs_f_1
	(struct Vector_d* rhs, struct Matrix_d** lhs_ptr, const struct DPG_Solver_Volume* dpg_s_vol,
	 const struct Simulation* sim)
{
	struct Matrix_d* lhs = *lhs_ptr;

	const struct Volume* vol          = (struct Volume*) dpg_s_vol;
	const struct Solver_Volume* s_vol = (struct Solver_Volume*) dpg_s_vol;

	const int n_eq = sim->test_case->n_eq,
	          n_vr = sim->test_case->n_var;
	const ptrdiff_t ext_0 = (rhs->ext_0)/n_eq;

	const ptrdiff_t n_dof_s  = (lhs->ext_1)/n_vr,
	                n_dof_nf = compute_n_dof_nf(s_vol);
	struct Matrix_d* lhs_add = constructor_empty_Matrix_d('R',lhs->ext_0,(n_dof_s+n_dof_nf)*n_vr); // moved
//set_to_value_Matrix_d(lhs_add,0.0);
	set_block_Matrix_d(lhs_add,(struct const_Matrix_d*)lhs,0,0,'i');
	destructor_Matrix_d(lhs);
	lhs = lhs_add;
//print_Matrix_d(lhs);

	struct Matrix_d rhs_M = { .layout = 'C', .ext_0 = ext_0, .ext_1 = n_eq, .owns_data = false, .data = rhs->data, };

	int ind_dof = n_vr*n_dof_s;
	for (int i = 0; i < NFMAX;    ++i) {
	for (int j = 0; j < NSUBFMAX; ++j) {
		const struct Face* face = vol->faces[i][j];
		if (!face)
			continue;

		const struct DPG_Solver_Face* dpg_s_face = (struct DPG_Solver_Face*) face;
		if (!face->boundary)
			increment_rlhs_internal_face(dpg_s_vol,dpg_s_face,lhs,&rhs_M,&ind_dof,sim);
		else
			increment_rlhs_boundary_face(dpg_s_vol,dpg_s_face,lhs,&rhs_M,sim);
	}}

//print_Matrix_d(lhs);
	*lhs_ptr = lhs;
}

static void increment_rhs_source
	(struct Vector_d* rhs, const struct Solver_Volume* s_vol, const struct Simulation* sim)
{
	const int n_eq = sim->test_case->n_eq;
	ptrdiff_t extents[2] = { (rhs->ext_0)/n_eq, n_eq, };
	struct Multiarray_d rhs_Ma =
		{ .layout = 'C', .order = 2, .extents = extents, .owns_data = rhs->owns_data, .data = rhs->data };

	sim->test_case->compute_source_rhs(sim,s_vol,&rhs_Ma);
}

static void add_to_petsc_Mat_Vec_dpg
	(const struct Solver_Volume* s_vol, const struct const_Vector_d* rhs_neg, const struct const_Matrix_d* lhs,
	 struct Solver_Storage_Implicit* ssi, const struct Simulation* sim)
{
	assert(sizeof(int) == sizeof(PetscInt));
	assert(sizeof(double) == sizeof(PetscScalar));

	const ptrdiff_t ext_0 = rhs_neg->ext_0;

	const struct const_Vector_i* idxm = constructor_petsc_idxm_dpg(ext_0,s_vol); // destructed.

	if (sim->test_case->use_schur_complement) {
		const ptrdiff_t dof_s = compute_size(s_vol->sol_coef->order,s_vol->sol_coef->extents),
		                dof_g = compute_size(s_vol->grad_coef->order,s_vol->grad_coef->extents);
		invert_sub_block_Matrix_d((struct Matrix_d*)lhs,0,0,dof_s);
		assert(dof_g == 0); // Add support.
	}

	MatSetValues(ssi->A,ext_0,idxm->data,ext_0,idxm->data,lhs->data,ADD_VALUES);
	VecSetValues(ssi->b,ext_0,idxm->data,rhs_neg->data,ADD_VALUES);

	destructor_const_Vector_i(idxm);
}

// Level 3 ********************************************************************************************************** //

/// \brief Scale required \ref Numerical_Flux terms by the face Jacobian.
static void scale_by_Jacobian
	(const struct Numerical_Flux* num_flux, ///< \ref Numerical_Flux.
	 const struct Solver_Face* s_face       ///< The current \ref Solver_Face.
	);

/// \brief Increment the rhs terms with the contribution of the boundary face.
static void increment_rhs_boundary_face
	(struct Matrix_d* rhs,                  ///< Holds the rhs terms.
	 const struct Numerical_Flux* num_flux, ///< \ref Numerical_Flux.
	 const struct Solver_Face* s_face,      ///< The current \ref Solver_Face.
	 const struct Simulation* sim           ///< \ref Simulation.
	);

/// \brief Increment the lhs terms with the contribution of the boundary face.
static void increment_lhs_boundary_face
	(struct Matrix_d* lhs,                  ///< Holds the lhs terms.
	 const struct Numerical_Flux* num_flux, ///< \ref Numerical_Flux.
	 const struct Solver_Face* s_face,      ///< The current \ref Solver_Face.
	 const struct Simulation* sim           ///< \ref Simulation.
	);

static void increment_rlhs_internal_face
	(const struct DPG_Solver_Volume* dpg_s_vol, const struct DPG_Solver_Face* dpg_s_face, struct Matrix_d* lhs,
	 struct Matrix_d* rhs, int* ind_dof, const struct Simulation* sim)
{
	/// As the rhs is **always** linear wrt the trace unknowns, the rhs and lhs are computed together.
	const struct const_Matrix_d* lhs_l = constructor_lhs_l_internal_face_dpg(dpg_s_vol,dpg_s_face); // destructed

//print_const_Matrix_d(tw0_vt_fc_op->op_std);
	const struct Solver_Face* s_face  = (struct Solver_Face*) dpg_s_face;
	struct Matrix_d nf_coef = interpret_Multiarray_as_Matrix_d(s_face->nf_coef);

	mm_d('N','N',1.0,1.0,lhs_l,(struct const_Matrix_d*)&nf_coef,rhs);

	const int n_eq = sim->test_case->n_eq,
	          n_vr = sim->test_case->n_var;

	const struct Solver_Volume* s_vol = (struct Solver_Volume*) dpg_s_vol;
	const ptrdiff_t n_dof_test = (lhs->ext_0)/n_eq,
	                n_dof_nf   = compute_n_dof_nf(s_vol);

//printf("lhs_l\n");
//print_const_Matrix_d(lhs_l);
	for (int vr = 0; vr < n_vr; ++vr)
		set_block_Matrix_d(lhs,lhs_l,vr*n_dof_test,*ind_dof+vr*n_dof_nf,'i');
	destructor_const_Matrix_d(lhs_l);

	*ind_dof += nf_coef.ext_0;
}

static void increment_rlhs_boundary_face
	(const struct DPG_Solver_Volume* dpg_s_vol, const struct DPG_Solver_Face* dpg_s_face, struct Matrix_d* lhs,
	 struct Matrix_d* rhs, const struct Simulation* sim)
{
	struct Numerical_Flux_Input* num_flux_i = constructor_Numerical_Flux_Input(sim); // destructed

	const struct Solver_Face* s_face = (struct Solver_Face*) dpg_s_face;
	constructor_Numerical_Flux_Input_data(num_flux_i,s_face,sim); // destructed

	struct Numerical_Flux* num_flux = constructor_Numerical_Flux(num_flux_i); // destructed
	destructor_Numerical_Flux_Input_data(num_flux_i);
	destructor_Numerical_Flux_Input(num_flux_i);

	scale_by_Jacobian(num_flux,s_face);

	increment_rhs_boundary_face(rhs,num_flux,s_face,sim);
	increment_lhs_boundary_face(lhs,num_flux,s_face,sim);

	destructor_Numerical_Flux(num_flux);

	UNUSED(dpg_s_vol);
}

// Level 4 ********************************************************************************************************** //

static void scale_by_Jacobian (const struct Numerical_Flux* num_flux, const struct Solver_Face* s_face)
{
	const struct Face* face = (struct Face*) s_face;

	assert(face->boundary);
	assert(num_flux->neigh_info[0].dnnf_ds != NULL || num_flux->neigh_info[0].dnnf_dg != NULL);

	const struct const_Vector_d jacobian_det_fc = interpret_const_Multiarray_as_Vector_d(s_face->jacobian_det_fc);
	scale_Multiarray_by_Vector_d('L',1.0,(struct Multiarray_d*)num_flux->nnf,&jacobian_det_fc,false);

	if (num_flux->neigh_info[0].dnnf_ds)
		scale_Multiarray_by_Vector_d('L',1.0,(struct Multiarray_d*)num_flux->neigh_info[0].dnnf_ds,&jacobian_det_fc,false);
	if (num_flux->neigh_info[0].dnnf_dg)
		EXIT_ADD_SUPPORT;
}

static void increment_rhs_boundary_face
	(struct Matrix_d* rhs, const struct Numerical_Flux* num_flux, const struct Solver_Face* s_face,
	 const struct Simulation* sim)
{
	ptrdiff_t extents[2] = { rhs->ext_0, rhs->ext_1, };
	struct Multiarray_d rhs_Ma =
		{ .layout = 'C', .order = 2, .extents = extents, .owns_data = false, .data = rhs->data, };

	const struct Operator* tw0_vt_fc = get_operator__tw0_vt_fc(0,s_face);

UNUSED(sim);
	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	const char op_format = 'd';
	mm_NNC_Operator_Multiarray_d(-1.0,1.0,tw0_vt_fc,num_flux->nnf,&rhs_Ma,op_format,2,NULL,NULL);
}

static void increment_lhs_boundary_face
	(struct Matrix_d* lhs, const struct Numerical_Flux* num_flux, const struct Solver_Face* s_face,
	 const struct Simulation* sim)
{
	UNUSED(sim);
	assert(((struct Face*)s_face)->boundary);

	struct Matrix_d* lhs_ll = constructor_lhs_f_1((int[]){0,0},num_flux,s_face); // destructed

	set_block_Matrix_d(lhs,(struct const_Matrix_d*)lhs_ll,0,0,'a');
//printf("lhs\n");
//print_Matrix_d(lhs_ll);
//print_Matrix_d(lhs);
	destructor_Matrix_d(lhs_ll);
}
