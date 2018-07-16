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
#include "definitions_numerical_flux.h"

#include "def_templates_compute_face_rlhs_fr_split_form.h"

#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"


#include "def_templates_face_solver.h"
#include "def_templates_face_solver_fr_split_form.h"
#include "def_templates_volume_solver.h"
#include "def_templates_volume_solver_fr_split_form.h"

#include "def_templates_boundary.h"
#include "def_templates_compute_face_rlhs.h"
#include "def_templates_numerical_flux.h"
#include "def_templates_operators.h"
#include "def_templates_solve_fr_split_form.h"
#include "def_templates_test_case.h"

#include "def_templates_flux.h"

// Static function declarations ************************************************************************************* //

/// \brief Container for solver-related parameters.
struct S_Params_T {
	scale_by_Jacobian_fptr_T scale_by_Jacobian; ///< Pointer to the appropriate function.

	compute_rlhs_f_fptr_T compute_rlhs; ///< Pointer to the appropriate function.
};

/// \brief Container for numerical flux related parameters.
struct Num_Flux_T {
	const struct const_Multiarray_T* n_dot_nf; ///< Unit normal dotted with the numerical flux.
};

/** \brief Set the parameters of \ref S_Params_T.
 *  \return A statically allocated \ref S_Params_T container. */
static struct S_Params_T set_s_params_T
	(const struct Simulation* sim ///< \ref Simulation.
	);
//im adding on now to attempt to compile
static void constructor_Boundary_Value_Input_g_face_fcl
	(struct Boundary_Value_Input_T*const bv_i,     ///< \ref Boundary_Value_Input_T.
	 const struct FRSF_Solver_Face_T*const frsf_s_face ///< \ref FRSF_Solver_Face_T.
	 );
static void constructor_Boundary_Value_g_face_fcl
	(struct Boundary_Value_T*const bv,             ///< \ref Boundary_Value_Input_T.
	 const struct FRSF_Solver_Face_T*const frsf_s_face ///< \ref FRSF_Solver_Face_T.
	 );
//my attempt pt 2
static void compute_rhs_f_frsf_like_T
 	(const struct Numerical_Flux_T*const num_flux, struct Solver_Face_T*const s_face,
 	 struct Solver_Storage_Implicit*const ssi);
static void finalize_face_rhs_frsf_like_T
	(const int side_index,                         ///< The index of the side of the face under consideration.
	 const struct Numerical_Flux_T*const num_flux, ///< Defined for \ref compute_rlhs_f_fptr_T.
	 struct Solver_Face_T*const s_face             ///< Defined for \ref compute_rlhs_f_fptr_T.
	 );
//end of the added attempt
/* // Interface functions ********************************************************************************************** // */

void compute_face_rlhs_fr_split_form_T
	(const struct Simulation* sim, struct Solver_Storage_Implicit* ssi, struct Intrusive_List* faces)
{
	assert(sim->elements->name == IL_ELEMENT_SOLVER_FRSF);
	assert(sim->faces->name    == IL_FACE_SOLVER_FRSF);
	assert(sim->volumes->name  == IL_VOLUME_SOLVER_FRSF);


	//Trying an idea
	/* UNUSED(ssi); */
	/* UNUSED(faces); */
	/* set_s_params_T(sim); */

	const bool has_2nd_order = get_set_has_1st_2nd_order(NULL)[1];

	struct S_Params_T s_params = set_s_params_T(sim);
	struct Numerical_Flux_Input_T* num_flux_i = constructor_Numerical_Flux_Input_T(sim); // destructed

	for (struct Intrusive_Link* curr = faces->first; curr; curr = curr->next) {
		struct Solver_Face_T*const s_face       = (struct Solver_Face_T*) curr;
		struct FRSF_Solver_Face_T*const frsf_s_face = (struct FRSF_Solver_Face_T*) curr;

		/* if (curr == faces->first) */
		/* 	continue; */
		/* else */
		/* 	if (curr == faces->last) */
		/* 		break; */
		/* 	else{ */
		constructor_Numerical_Flux_Input_data_frsf_T(num_flux_i,frsf_s_face,sim,has_2nd_order); // destructed

		struct Numerical_Flux_T* num_flux = constructor_Numerical_Flux_T(num_flux_i); // destructed
		destructor_Numerical_Flux_Input_data_T(num_flux_i);

		/* s_params.scale_by_Jacobian(num_flux,s_face); */
		s_params.compute_rlhs(num_flux,s_face,ssi);
		destructor_Numerical_Flux_T(num_flux);
			/* } */
	}
	//add on send
	/* struct Intrusive_Link* curr = faces->first; */
	/* struct Solver_Face_T*const s_face       = (struct Solver_Face_T*) curr; */
	/* struct FRSF_Solver_Face_T*const frsf_s_face = (struct FRSF_Solver_Face_T*) curr; */

	/* constructor_Numerical_Flux_Input_data_frsf_T(num_flux_i,frsf_s_face,sim,has_2nd_order); // destructed */

	/* struct Numerical_Flux_T* num_flux = constructor_Numerical_Flux_T(num_flux_i); // destructed */
	/* destructor_Numerical_Flux_Input_data_T(num_flux_i); */

	/* /\* s_params.scale_by_Jacobian(num_flux,s_face); *\/ */
	/* const struct Face* face = (struct Face*) s_face; */
	/* permute_Multiarray_T_fc((struct Multiarray_T*)num_flux_i->nnf,'R',1,(struct Solver_Face_T*)face); */
	/* scale_Multiarray_T((struct Multiarray_T*)num_flux_i->nnf,-1.0); // Use "-ve" normal. */
	/* finalize_face_rhs_frsf_like_T(1,num_flux_i,s_face); */
	/* destructor_Numerical_Flux_T(num_flux); */



	//end of full send




	destructor_Numerical_Flux_Input_T(num_flux_i);
}


void constructor_Numerical_Flux_Input_data_frsf_T
	(struct Numerical_Flux_Input_T*const num_flux_i, const struct FRSF_Solver_Face_T*const frsf_s_face,
	 const struct Simulation*const sim, const bool compute_gradient)
{
	/* UNUSED(compute_gradient); */
	const struct Solver_Face_T*const s_face = (struct Solver_Face_T*) frsf_s_face;

	if (compute_gradient) {
		constructor_Boundary_Value_Input_g_face_fcl(&num_flux_i->bv_l,frsf_s_face); // destructed
		constructor_Boundary_Value_g_face_fcl(&num_flux_i->bv_r,frsf_s_face);       // destructed
	}
	constructor_Numerical_Flux_Input_data_T(num_flux_i,s_face,sim); // destructed
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Constructor for the partially corrected weak gradient interpolated to the face cubature nodes as seen from
 *         the volume of input "side_index".
 *  \return See brief. */
struct Multiarray_T* constructor_copy_Multiarray_T_frsf (struct Multiarray_T* src, Type dif)
{
	const ptrdiff_t size = compute_size(src->order,src->extents);
	Type* data = malloc((size_t)size * sizeof *data); // moved
	for (int i = 0; i < size; ++i)
		data[i]=dif;

	return constructor_move_Multiarray_T_T(src->layout,src->order,src->extents,true,data);
}

const struct const_Multiarray_T* constructor_copy_const_Multiarray_T_frsf (const struct const_Multiarray_T*const src, Type dif)
{
	return (const struct const_Multiarray_T*) constructor_copy_Multiarray_T_frsf((struct Multiarray_T*)src, dif);
}
static void finalize_face_rhs_frsf_like_T
	(const int side_index, const struct Numerical_Flux_T*const num_flux, struct Solver_Face_T*const s_face)
{
	const struct Face*const face = (struct Face*) s_face;

	const struct Operator* tw0_vt_fc = get_operator__tw0_vt_fc_T(side_index,s_face);

	const char op_format = get_set_op_format(0);
	/* printf("%lf\n",num_flux->nnf->data[0]); */
	/* print_const_Multiarray_T(num_flux->nnf); */
	struct Solver_Volume_T*const s_vol = (struct Solver_Volume_T*) face->neigh_info[side_index].volume;

	struct Multiarray_d* sol_coef = s_vol->sol_coef;

	ptrdiff_t size = compute_size(sol_coef->order,sol_coef->extents);
	struct Vector_d* u2 = constructor_empty_Vector_d(size);
	for (int i = 0; i < size; ++i) {
		u2->data[i] = 0.5*sol_coef->data[i]*sol_coef->data[i];
	}
	/* print_Vector_d(u2); */
	Type dif;//this calculates the difference for DG strong form fnum-R*f at the faces' boundaries
	/* printf("side %d\n",side_index); */
	if (face->index == 0) {
		if (side_index==0){
			dif=num_flux->nnf->data[0]+u2->data[0];}
		else{
			dif=num_flux->nnf->data[0]-u2->data[size-1];}
		dif*=-1;
	} else {
		if (side_index==0){
			dif=num_flux->nnf->data[0]-u2->data[size-1];}
		else{
			dif=num_flux->nnf->data[0]+u2->data[0];}
	}
	/* printf("The dif %f\n",dif); */
	/* print_Multiarray_d(s_vol->rhs); */


	const struct const_Multiarray_T* num_flux_new =constructor_copy_const_Multiarray_T_frsf(num_flux->nnf,dif);//copies but changes that data value temporarily to my calculated dif
	// print_const_Multiarray_d(num_flux_new);
	/* printf("new flux %lf\n",num_flux_new->data[0]); */

	/* mm_NNC_Operator_Multiarray_T(-1.0,1.0,tw0_vt_fc,num_flux->nnf,s_vol->rhs,op_format,2,NULL,NULL); */
	mm_NNC_Operator_Multiarray_T(-1.0,1.0,tw0_vt_fc,num_flux_new,s_vol->rhs,op_format,2,NULL,NULL);
	/* mm_NNC_Operator_Multiarray_T(-1.0,0.0,tw0_vt_fc,num_flux_new,s_vol->rhs,op_format,2,NULL,NULL); */
	/* print_const_Matrix_d(tw0_vt_fc->op_std); */
	/* EXIT; */

	/* printf("%d\n",face->index); */
	/* print_Multiarray_d(s_vol->rhs); */
	/* EXIT; */
	/* struct FRSF_Solver_Volume_T* frsf_s_vol = (struct FRSF_Solver_Volume_T*) s_vol; */
	/* /\* print_const_Matrix_T(frsf_s_vol->m_inv);//confirmed mass matrix perf *\/ */
	/* for (int i = 0; i < size; ++i){ */
	/* 	double k=(size)*i; */
	/* 	int j; */
	/* 	j = (int)(k); */
	/* 	s_vol->rhs->data[i] = frsf_s_vol->m_inv->data[j]*s_vol->rhs->data[i];//multiply rhs by M_inv, works b/c M diag */
	/* } */
	/* print_Multiarray_d(s_vol->sol_coef); */
	/* print_Multiarray_d(s_vol->rhs); */
	//currently printing to see my values and that everything is working well
	destructor_const_Multiarray_T(num_flux_new);//destructed
}

static void compute_rhs_f_frsf_like_T
	(const struct Numerical_Flux_T*const num_flux, struct Solver_Face_T*const s_face,
	 struct Solver_Storage_Implicit*const ssi)
{
	UNUSED(ssi);

	const struct Face* face = (struct Face*) s_face;
	finalize_face_rhs_frsf_like_T(0,num_flux,s_face);
	if (!face->boundary) {
		permute_Multiarray_T_fc((struct Multiarray_T*)num_flux->nnf,'R',1,(struct Solver_Face_T*)face);
		scale_Multiarray_T((struct Multiarray_T*)num_flux->nnf,-1.0); // Use "-ve" normal.
		finalize_face_rhs_frsf_like_T(1,num_flux,s_face);
	}
}

static struct S_Params_T set_s_params_T (const struct Simulation* sim)
{
	struct S_Params_T s_params;

	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	switch (test_case->solver_method_curr) {
	case 'e':
		/* s_params.scale_by_Jacobian = scale_by_Jacobian_e_T; */
		s_params.compute_rlhs      = compute_rhs_f_frsf_like_T;

		//me added
		/* const struct const_Matrix_d*m_inv = sim->m_inv; */
		/* s_params.compute_rlhs      = constructor_mm_const_Matrix_d('N','N',1.0,m_inv,s_params.compute_rlhs,'R'); */
		/* scale_rhs_by_m_inv(s_params); */
		//end me added
		break;
	case 'i':
	default:
		EXIT_ERROR("Unsupported: %c (type_rc: %d)\n",test_case->solver_method_curr,TYPE_RC);
		break;
	}

	return s_params;
}
//add on attempt
static struct Multiarray_T* constructor_partial_grad_fc_interp
	(const int side_index,                         ///< The index of the side under consideration.
	 const struct FRSF_Solver_Face_T*const frsf_s_face ///< \ref FRSF_Solver_Face_T.
	 );
static void constructor_Boundary_Value_Input_g_face_fcl
	(struct Boundary_Value_Input_T*const bv_i, const struct FRSF_Solver_Face_T*const frsf_s_face)
{
	bv_i->g = (struct const_Multiarray_T*) constructor_partial_grad_fc_interp(0,frsf_s_face); // destructed
}

static void constructor_Boundary_Value_g_face_fcl
	(struct Boundary_Value_T*const bv, const struct FRSF_Solver_Face_T*const frsf_s_face)
{
	const struct Face*const face = (struct Face*) frsf_s_face;
	if (face->boundary)
		return;

	const int side_index = 1;
	const struct Solver_Face_T*const s_face = (struct Solver_Face_T*) frsf_s_face;

	struct Multiarray_T* grad_r_fcr = constructor_partial_grad_fc_interp(side_index,frsf_s_face); // moved
	permute_Multiarray_T_fc(grad_r_fcr,'R',side_index,s_face);

	bv->g = (struct const_Multiarray_T*) grad_r_fcr; // destructed
}
// Level 1 ********************************************************************************************************** //

static Real compute_scaling_weak_gradient
	(const struct FRSF_Solver_Face_T*const frsf_s_face ///< \ref FRSF_Solver_Face_T.
	 );
static struct Multiarray_T* constructor_partial_grad_fc_interp
	(const int side_index, const struct FRSF_Solver_Face_T*const frsf_s_face)
{

	const struct Face*const face              = (struct Face*) frsf_s_face;
	const struct Solver_Face_T*const s_face   = (struct Solver_Face_T*) frsf_s_face;
	const struct Volume*const vol             = (struct Volume*) face->neigh_info[side_index].volume;
	const struct FRSF_Solver_Volume_T* frsf_s_vol = (struct FRSF_Solver_Volume_T*) vol;

	const struct Neigh_Info_FRSF*const ni = &frsf_s_face->neigh_info[side_index];

	const Real scale = compute_scaling_weak_gradient(frsf_s_face);
	const struct Multiarray_T*const g_coef_p =
		constructor_sum_Multiarrays_Multiarray_T(1.0,frsf_s_vol->grad_coef_v,scale,ni->grad_coef_f); // destructed


	const struct Operator*const cv0_vr_fc = get_operator__cv0_vr_fc_T(side_index,s_face);

	const char op_format = get_set_op_format(0);

	struct Multiarray_T* ret =
		constructor_mm_NN1_Operator_Multiarray_T(cv0_vr_fc,g_coef_p,'C',op_format,g_coef_p->order,NULL); // returned
	destructor_Multiarray_T((struct Multiarray_T*)g_coef_p);
	return ret;
}

// Level 2 ********************************************************************************************************** //

#define PENALTY_SCALING 1.01

static Real compute_scaling_weak_gradient (const struct FRSF_Solver_Face_T*const frsf_s_face)
{
	const int ind_num_flux_2nd = get_set_ind_num_flux(NULL)[1];

	const struct Face*const face = (struct Face*) frsf_s_face;
	switch (ind_num_flux_2nd) {
	case NUM_FLUX_BR2_STABLE:
		return PENALTY_SCALING*NFMAX;
		break;
	case NUM_FLUX_CDG2:
		if (face->boundary)
			return PENALTY_SCALING*NFMAX;
		else
			EXIT_ADD_SUPPORT; // 0.0 for one side, 2.0*PENALTY_SCALING*NFMAX for larger area side.
		break;
	default:
		EXIT_ERROR("Unsupported: %d",ind_num_flux_2nd);
		break;
	}
}

#include "undef_templates_compute_face_rlhs_fr_split_form.h"

#include "undef_templates_matrix.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"

#include "undef_templates_face_solver.h"
#include "undef_templates_face_solver_fr_split_form.h"
#include "undef_templates_volume_solver.h"
#include "undef_templates_volume_solver_fr_split_form.h"

#include "undef_templates_boundary.h"
#include "undef_templates_compute_face_rlhs.h"
#include "undef_templates_numerical_flux.h"
#include "undef_templates_operators.h"
#include "undef_templates_solve_fr_split_form.h"
#include "undef_templates_test_case.h"
