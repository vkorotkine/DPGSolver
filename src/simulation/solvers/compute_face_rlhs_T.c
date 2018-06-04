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

#include "macros.h"


#include "def_templates_compute_face_rlhs.h"

#include "def_templates_face_solver.h"
#include "def_templates_volume_solver.h"

#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"

#include "def_templates_boundary.h"
#include "def_templates_numerical_flux.h"
#include "def_templates_operators.h"
#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

#define PRINT_OPERATORS 0 ///< Enable to print operators to terminal.

/// \brief Finalize the rhs term contribution from the \ref Face.
static void finalize_face_rhs_dg_like_T
	(const int side_index,                         ///< The index of the side of the face under consideration.
	 const struct Numerical_Flux_T*const num_flux, ///< Defined for \ref compute_rlhs_f_fptr_T.
	 struct Solver_Face_T*const s_face             ///< Defined for \ref compute_rlhs_f_fptr_T.
	);

// Interface functions ********************************************************************************************** //

const struct Operator* get_operator__tw0_vt_fc_T (const int side_index, const struct Solver_Face_T* s_face)
{
	const struct Face* face             = (struct Face*) s_face;
	const struct Volume* vol            = face->neigh_info[side_index].volume;
	const struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) vol;
	const struct Solver_Element* e      = (struct Solver_Element*) vol->element;

	const int ind_lf   = face->neigh_info[side_index].ind_lf,
	          ind_href = face->neigh_info[side_index].ind_href;
	const int p_v = s_vol->p_ref,
	          p_f = s_face->p_ref;

	const int curved = ( (s_face->cub_type == 's') ? 0 : 1 );
	return get_Multiarray_Operator(e->tw0_vt_fc[curved],(ptrdiff_t[]){ind_lf,ind_href,0,p_f,p_v});
}

const struct Operator* get_operator__cv0_vs_fc_T (const int side_index, const struct Solver_Face_T* s_face)
{
	const struct Face* face             = (struct Face*) s_face;
	const struct Volume* vol            = face->neigh_info[side_index].volume;
	const struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) vol;
	const struct Solver_Element* e      = (struct Solver_Element*) vol->element;

	const int ind_lf   = face->neigh_info[side_index].ind_lf,
	          ind_href = face->neigh_info[side_index].ind_href;
	const int p_v = s_vol->p_ref,
	          p_f = s_face->p_ref;

	const int curved = ( (s_face->cub_type == 's') ? 0 : 1 );
	return get_Multiarray_Operator(e->cv0_vs_fc[curved],(ptrdiff_t[]){ind_lf,ind_href,0,p_f,p_v});
}

const struct Operator* get_operator__cv0_vr_fc_T (const int side_index, const struct Solver_Face_T* s_face)
{
	const struct Face* face             = (struct Face*) s_face;
	const struct Volume* vol            = face->neigh_info[side_index].volume;
	const struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) vol;
	const struct Solver_Element* s_e    = (struct Solver_Element*) vol->element;

	const int ind_lf   = face->neigh_info[side_index].ind_lf,
	          ind_href = face->neigh_info[side_index].ind_href;
	const int p_v = s_vol->p_ref,
	          p_f = s_face->p_ref;

	const int curved = ( (s_face->cub_type == 's') ? 0 : 1 );
	return get_Multiarray_Operator(s_e->cv0_vr_fc[curved],(ptrdiff_t[]){ind_lf,ind_href,0,p_f,p_v});
}

void permute_Matrix_T_fc
	(struct Matrix_T* data, const char perm_layout, const int side_index_dest, const struct Solver_Face_T* s_face)
{
	assert(perm_layout == 'R');
	assert(data->layout == 'R');

	const struct const_Vector_i* nc_fc = get_operator__nc_fc_T(side_index_dest,s_face);
	permute_Matrix_T_V(data,nc_fc);
}
#if TYPE_RC == TYPE_COMPLEX
void permute_Matrix_R_fc
	(struct Matrix_R* data, const char perm_layout, const int side_index_dest, const struct Solver_Face_T* s_face)
{
	assert(perm_layout == 'R');
	assert(data->layout == 'R');

	const struct const_Vector_i* nc_fc = get_operator__nc_fc_T(side_index_dest,s_face);
	permute_Matrix_R_V(data,nc_fc);
}
#endif

const struct const_Vector_i* get_operator__nc_fc_T (const int side_index_dest, const struct Solver_Face_T* s_face)
{
	const struct Neigh_Info* neigh_info = &((struct Face*)s_face)->neigh_info[side_index_dest];

	struct Volume* vol = neigh_info->volume;
	const struct Solver_Element* e = (struct Solver_Element*) vol->element;

	const int ind_ord = neigh_info->ind_ord,
	          ind_e   = get_face_element_index((struct Face*)s_face),
	          p_f     = s_face->p_ref;
	const int curved = ( (s_face->cub_type == 's') ? 0 : 1 );

	return get_const_Multiarray_Vector_i(e->nc_fc[curved],(ptrdiff_t[]){ind_ord,ind_e,ind_e,0,0,p_f,p_f});
}

void constructor_Numerical_Flux_Input_data_T
	(struct Numerical_Flux_Input_T* num_flux_i, const struct Solver_Face_T* s_face, const struct Simulation* sim)
{
	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;

	test_case->constructor_Boundary_Value_Input_face_fcl(&num_flux_i->bv_l,s_face,sim);     // destructed
	s_face->constructor_Boundary_Value_fcl(&num_flux_i->bv_r,&num_flux_i->bv_l,s_face,sim); // destructed
}

void destructor_Numerical_Flux_Input_data_T (struct Numerical_Flux_Input_T* num_flux_i)
{
	destructor_Boundary_Value_Input_T(&num_flux_i->bv_l);
	destructor_Boundary_Value_T(&num_flux_i->bv_r);
}

struct Matrix_T* constructor_lhs_f_1_T
	(const int side_index[2], const struct Numerical_Flux_T* num_flux, const struct Solver_Face_T* s_face)
{
	const struct Operator* tw0_vt_fc    = get_operator__tw0_vt_fc_T(side_index[0],s_face);
	const struct Operator* cv0_vs_fc_op = get_operator__cv0_vs_fc_T(side_index[1],s_face);
#if PRINT_OPERATORS == 1
printf("face\n");
print_const_Matrix_d(tw0_vt_fc->op_std);
print_const_Matrix_d(cv0_vs_fc_op->op_std);
#endif

	const struct const_Matrix_R* cv0_vs_fc = cv0_vs_fc_op->op_std;
	bool need_free_cv0 = false;
	if (side_index[0] != side_index[1]) {
		need_free_cv0 = true;
		cv0_vs_fc = constructor_copy_const_Matrix_R(cv0_vs_fc); // destructed
		permute_Matrix_R_fc((struct Matrix_R*)cv0_vs_fc,'R',side_index[0],s_face);
	}

	const struct const_Multiarray_T* dnnf_ds_Ma = num_flux->neigh_info[side_index[1]].dnnf_ds;

	const ptrdiff_t ext_0 = tw0_vt_fc->op_std->ext_0,
	                ext_1 = tw0_vt_fc->op_std->ext_1;
	const int n_eq = (int)dnnf_ds_Ma->extents[1],
	          n_vr = (int)dnnf_ds_Ma->extents[2];

	struct Matrix_T* tw0_nf = constructor_empty_Matrix_T('R',ext_0,ext_1);                         // destructed
	struct Matrix_T* lhs_l  = constructor_empty_Matrix_T('R',ext_0,cv0_vs_fc->ext_1);              // destructed
	struct Matrix_T* lhs    = constructor_empty_Matrix_T('R',n_eq*lhs_l->ext_0,n_vr*lhs_l->ext_1); // returned
	set_to_value_Matrix_T(tw0_nf,0.0);

	struct Vector_T dnnf_ds = { .ext_0 = dnnf_ds_Ma->extents[0], .owns_data = false, .data = NULL, };

	for (int vr = 0; vr < n_vr; ++vr) {
	for (int eq = 0; eq < n_eq; ++eq) {
		const ptrdiff_t ind =
			compute_index_sub_container(dnnf_ds_Ma->order,1,dnnf_ds_Ma->extents,(ptrdiff_t[]){eq,vr});
		dnnf_ds.data = (Type*)&dnnf_ds_Ma->data[ind];
		mm_diag_T('R',1.0,0.0,tw0_vt_fc->op_std,(struct const_Vector_T*)&dnnf_ds,tw0_nf,false);

		mm_TRT('N','N',-1.0,0.0,(struct const_Matrix_T*)tw0_nf,cv0_vs_fc,lhs_l);

		set_block_Matrix_T(lhs,eq*lhs_l->ext_0,vr*lhs_l->ext_1,
		                   (struct const_Matrix_T*)lhs_l,0,0,lhs_l->ext_0,lhs_l->ext_1,'i');
#if 0
printf("eq,var: %d %d\n",eq,vr);
if (eq == 0 && vr == 0)
print_Matrix_T(lhs_l);
#endif
	}}
	destructor_Matrix_T(tw0_nf);
	destructor_Matrix_T(lhs_l);

	if (need_free_cv0)
		destructor_const_Matrix_R(cv0_vs_fc);

	return lhs;
}

struct Matrix_T* constructor_lhs_p_f_2_T
	(const int side_index[2], const struct Numerical_Flux_T* num_flux, const struct Solver_Face_T* s_face)
{
	const struct Operator* tw0_vt_fc    = get_operator__tw0_vt_fc_T(side_index[0],s_face);
	const struct Operator* cv0_vr_fc_op = get_operator__cv0_vr_fc_T(side_index[1],s_face);

	const struct const_Matrix_R* cv0_vr_fc = cv0_vr_fc_op->op_std;
	bool need_free_cv0 = false;
	if (side_index[0] != side_index[1]) {
		need_free_cv0 = true;
		cv0_vr_fc = constructor_copy_const_Matrix_R(cv0_vr_fc); // destructed
		permute_Matrix_R_fc((struct Matrix_R*)cv0_vr_fc,'R',side_index[0],s_face);
	}

	const struct const_Multiarray_T* dnnf_dg_Ma = num_flux->neigh_info[side_index[1]].dnnf_dg;

	const ptrdiff_t ext_0 = tw0_vt_fc->op_std->ext_0,
	                ext_1 = tw0_vt_fc->op_std->ext_1;
	const int n_eq = (int)dnnf_dg_Ma->extents[1],
	          n_vr = (int)dnnf_dg_Ma->extents[2];

	struct Matrix_T* tw0_nf = constructor_empty_Matrix_T('R',ext_0,ext_1);                             // destructed
	struct Matrix_T* lhs_l  = constructor_empty_Matrix_T('R',ext_0,cv0_vr_fc->ext_1);                  // destructed
	struct Matrix_T* lhs_p  = constructor_empty_Matrix_T('R',n_eq*lhs_l->ext_0,n_vr*DIM*lhs_l->ext_1); // returned
	set_to_value_Matrix_T(tw0_nf,0.0);

	struct Vector_T dnnf_dg = { .ext_0 = dnnf_dg_Ma->extents[0], .owns_data = false, .data = NULL, };

	for (int d_g = 0; d_g < DIM; ++d_g) {
	for (int vr = 0; vr < n_vr; ++vr) {
	for (int eq = 0; eq < n_eq; ++eq) {
		const ptrdiff_t ind =
			compute_index_sub_container(dnnf_dg_Ma->order,1,dnnf_dg_Ma->extents,(ptrdiff_t[]){eq,vr,d_g});
		dnnf_dg.data = (Type*)&dnnf_dg_Ma->data[ind];
		mm_diag_T('R',1.0,0.0,tw0_vt_fc->op_std,(struct const_Vector_T*)&dnnf_dg,tw0_nf,false);

		mm_TRT('N','N',-1.0,0.0,(struct const_Matrix_T*)tw0_nf,cv0_vr_fc,lhs_l);
//print_Matrix_T(lhs_l);

		set_block_Matrix_T(lhs_p,eq*lhs_l->ext_0,(vr+n_vr*(d_g))*lhs_l->ext_1,
		                   (struct const_Matrix_T*)lhs_l,0,0,lhs_l->ext_0,lhs_l->ext_1,'i');
	}}}
//EXIT_UNSUPPORTED;
	destructor_Matrix_T(tw0_nf);
	destructor_Matrix_T(lhs_l);

	if (need_free_cv0)
		destructor_const_Matrix_R(cv0_vr_fc);

	return lhs_p;
}

void add_to_flux_imbalance_face_nf_w_T
	(const struct const_Matrix_T*const nf_fc, const struct const_Vector_R*const w_fc,
	 const struct Solver_Face_T*const s_face)
{
	const struct const_Matrix_T* nnf_integral =
		constructor_mm_diag_const_Matrix_T_R(1.0,nf_fc,w_fc,'L',false); // destructed

	const struct const_Vector_T* nnf_integral_sum =
		constructor_sum_const_Vector_T_const_Matrix_T('R',nnf_integral); // destructed
	destructor_const_Matrix_T(nnf_integral);

	const struct Face* face = (struct Face*) s_face;
	const int n_neigh    = ( face->boundary ? 1 : 2 );
	for (int n = 0; n < n_neigh; ++n) {
		struct Solver_Volume_T*const s_vol = (struct Solver_Volume_T*) face->neigh_info[n].volume;
		const Real normal_scale = ( n == 0 ? -1.0 : 1.0 );
		const ptrdiff_t n_vr = s_vol->flux_imbalance->ext_0;
		for (int vr = 0; vr < n_vr; ++vr)
			s_vol->flux_imbalance->data[vr] += normal_scale*nnf_integral_sum->data[vr];
	}
	destructor_const_Vector_T(nnf_integral_sum);
}

void compute_rhs_f_dg_like_T
	(const struct Numerical_Flux_T*const num_flux, struct Solver_Face_T*const s_face,
	 struct Solver_Storage_Implicit*const ssi)
{
	UNUSED(ssi);

	const struct Face* face = (struct Face*) s_face;
	finalize_face_rhs_dg_like_T(0,num_flux,s_face);
	if (!face->boundary) {
		permute_Multiarray_T_fc((struct Multiarray_T*)num_flux->nnf,'R',1,(struct Solver_Face_T*)face);
		scale_Multiarray_T((struct Multiarray_T*)num_flux->nnf,-1.0); // Use "-ve" normal.
		finalize_face_rhs_dg_like_T(1,num_flux,s_face);
	}
}

void permute_Multiarray_T_fc
	(struct Multiarray_T* data, const char perm_layout, const int side_index_dest,
	 const struct Solver_Face_T*const s_face)
{
	const struct const_Vector_i* nc_fc = get_operator__nc_fc_T(side_index_dest,s_face);
	permute_Multiarray_T_V(data,nc_fc,perm_layout);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void finalize_face_rhs_dg_like_T
	(const int side_index, const struct Numerical_Flux_T*const num_flux, struct Solver_Face_T*const s_face)
{
	const struct Face*const face = (struct Face*) s_face;

	const struct Operator* tw0_vt_fc = get_operator__tw0_vt_fc_T(side_index,s_face);

	const char op_format = get_set_op_format(0);

	struct Solver_Volume_T*const s_vol = (struct Solver_Volume_T*) face->neigh_info[side_index].volume;

	mm_NNC_Operator_Multiarray_T(-1.0,1.0,tw0_vt_fc,num_flux->nnf,s_vol->rhs,op_format,2,NULL,NULL);
}
