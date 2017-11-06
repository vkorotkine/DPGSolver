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

#include "compute_face_rlhs_dg.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"
#include "definitions_intrusive.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "face_solver_dg.h"
#include "element_solver_dg.h"
#include "volume.h"
#include "volume_solver_dg.h"

#include "multiarray_operator.h"
#include "numerical_flux.h"
#include "operator.h"
#include "simulation.h"
#include "solve_dg.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

struct Num_Flux;

/** \brief Function pointer to the function used to scale by the face Jacobian.
 *
 *  \param num_flux \ref Numerical_Flux.
 *  \param face     \ref Face.
 *  \param sim      \ref Simulation.
 */
typedef void (*scale_by_Jacobian_fptr)
	(const struct Numerical_Flux* num_flux,
	 struct Face* face,
	 const struct Simulation* sim
	);

/** \brief Function pointer to the function used to evaluate the rhs (and optionally lhs) terms.
 *
 *  \param num_flux  \ref Numerical_Flux.
 *  \param face      \ref Face.
 *  \param s_store_i \ref Solver_Storage_Implicit.
 *  \param sim       \ref Simulation.
 */
typedef void (*compute_rlhs_fptr)
	(const struct Numerical_Flux* num_flux,
	 struct Face* face,
	 struct Solver_Storage_Implicit* s_store_i,
	 const struct Simulation* sim
	);

/// \brief Container for solver-related parameters.
struct S_Params {
	scale_by_Jacobian_fptr scale_by_Jacobian; ///< Pointer to the appropriate function.

	compute_rlhs_fptr compute_rlhs; ///< Pointer to the appropriate function.
};

/// \brief Container for numerical flux related parameters.
struct Num_Flux {
	const struct const_Multiarray_d* n_dot_nf; ///< Unit normal dotted with the numerical flux.
};

/** \brief Constructor for the solution evaluated at the face cubature nodes.
 *  \return Standard. */
const struct const_Multiarray_d* constructor_sol_fc
	(struct Face* face,           ///< \ref Face.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Set the parameters of \ref S_Params.
 *  \return A statically allocated \ref S_Params container. */
static struct S_Params set_s_params
	(const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Construct the data members of the \ref Numerical_Flux_Input container which are specific to the face under
 *         consideration. */
void constructor_Numerical_Flux_Input_data
	(struct Numerical_Flux_Input* num_flux_i, ///< \ref Numerical_Flux_Input.
	 const struct Face* face,                 ///< \ref Face.
	 const struct Simulation* sim             ///< \ref Simulation.
	);

/// \brief Destructor for the data members of the \ref Numerical_Flux_Input container.
void destructor_Numerical_Flux_Input_data
	(struct Numerical_Flux_Input* num_flux_i ///< \ref Numerical_Flux_Input.
	);

// Interface functions ********************************************************************************************** //

void compute_face_rlhs_dg (const struct Simulation* sim, struct Solver_Storage_Implicit* s_store_i)
{
	assert(sim->elements->name == IL_ELEMENT_SOLVER_DG);
	assert(sim->faces->name    == IL_FACE_SOLVER_DG);
	assert(sim->volumes->name  == IL_VOLUME_SOLVER_DG);

	struct S_Params s_params = set_s_params(sim);
	struct Numerical_Flux_Input* num_flux_i = constructor_Numerical_Flux_Input(sim); // destructed

	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Face* face = (struct Face*) curr;
//printf("face: %d\n",face->index);

		constructor_Numerical_Flux_Input_data(num_flux_i,face,sim); // destructed
//print_const_Multiarray_d(num_flux_i->bv_l.s);
//print_const_Multiarray_d(num_flux_i->bv_r.s);

		struct Numerical_Flux* num_flux = constructor_Numerical_Flux(num_flux_i); // destructed
		destructor_Numerical_Flux_Input_data(num_flux_i);
//print_const_Multiarray_d(num_flux->nnf);
//print_const_Multiarray_d(num_flux->neigh_info[0].dnnf_ds);
//if (!face->boundary)
//	print_const_Multiarray_d(num_flux->neigh_info[1].dnnf_ds);

		s_params.scale_by_Jacobian(num_flux,face,sim);
//print_const_Multiarray_d(num_flux->nnf);
//print_const_Multiarray_d(num_flux->neigh_info[0].dnnf_ds);
//if (!face->boundary)
//	print_const_Multiarray_d(num_flux->neigh_info[1].dnnf_ds);
		s_params.compute_rlhs(num_flux,face,s_store_i,sim);
		destructor_Numerical_Flux(num_flux);
//if (face->index == 2)
//break;
	}
//EXIT_UNSUPPORTED;
	destructor_Numerical_Flux_Input(num_flux_i);
}

const struct Operator* get_operator__tw0_vs_fc__rlhs_dg (const int side_index, struct Face* face)
{
	struct Solver_Face* s_face = (struct Solver_Face*) face;
	struct Volume* vol         = (struct Volume*) face->neigh_info[side_index].volume;

	const struct DG_Solver_Element* e = (const struct DG_Solver_Element*) vol->element;

	const int ind_lf   = face->neigh_info[side_index].ind_lf,
	          ind_href = face->neigh_info[side_index].ind_href;
	const int p_v = ((struct Solver_Volume*)vol)->p_ref,
	          p_f = s_face->p_ref;

	const int curved = ( (s_face->cub_type == 's') ? 0 : 1 );
	return get_Multiarray_Operator(e->tw0_vs_fc[curved],(ptrdiff_t[]){ind_lf,ind_href,0,p_f,p_v});
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Scale \ref Numerical_Flux::nnf by the face Jacobian (i.e. only the explicit term).
static void scale_by_Jacobian_e
	(const struct Numerical_Flux* num_flux, ///< Defined for \ref scale_by_Jacobian_fptr.
	 struct Face* face,                     ///< Defined for \ref scale_by_Jacobian_fptr.
	 const struct Simulation* sim           ///< Defined for \ref scale_by_Jacobian_fptr.
	);

/// \brief Scale \ref Numerical_Flux::nnf and \ref Numerical_Flux::Neigh_Info_NF::dnnf_ds by the face Jacobian.
static void scale_by_Jacobian_i1
	(const struct Numerical_Flux* num_flux, ///< Defined for \ref scale_by_Jacobian_fptr.
	 struct Face* face,                     ///< Defined for \ref scale_by_Jacobian_fptr.
	 const struct Simulation* sim           ///< Defined for \ref scale_by_Jacobian_fptr.
	);

/// \brief Compute only the rhs term.
static void compute_rhs_f_dg
	(const struct Numerical_Flux* num_flux,     ///< Defined for \ref compute_rlhs_fptr.
	 struct Face* face,                         ///< Defined for \ref compute_rlhs_fptr.
	 struct Solver_Storage_Implicit* s_store_i, ///< Defined for \ref compute_rlhs_fptr.
	 const struct Simulation* sim               ///< Defined for \ref compute_rlhs_fptr.
	);

/// \brief Compute the rhs and first order lhs terms.
static void compute_rlhs_1
	(const struct Numerical_Flux* num_flux,     ///< Defined for \ref compute_rlhs_fptr.
	 struct Face* face,                         ///< Defined for \ref compute_rlhs_fptr.
	 struct Solver_Storage_Implicit* s_store_i, ///< Defined for \ref compute_rlhs_fptr.
	 const struct Simulation* sim               ///< Defined for \ref compute_rlhs_fptr.
	);

static struct S_Params set_s_params (const struct Simulation* sim)
{
	struct S_Params s_params;

	struct Test_Case* test_case = sim->test_case;
	switch (test_case->solver_method_curr) {
	case 'e':
		s_params.scale_by_Jacobian = scale_by_Jacobian_e;
		s_params.compute_rlhs      = compute_rhs_f_dg;
		break;
	case 'i':
		if (test_case->has_1st_order && !test_case->has_2nd_order) {
			s_params.scale_by_Jacobian = scale_by_Jacobian_i1;
			s_params.compute_rlhs      = compute_rlhs_1;
		} else if (!test_case->has_1st_order && test_case->has_2nd_order) {
			EXIT_ADD_SUPPORT; //s_params.compute_rlhs = compute_rlhs_2;
		} else if (test_case->has_1st_order && test_case->has_2nd_order) {
			EXIT_ADD_SUPPORT; //s_params.compute_rlhs = compute_rlhs_12;
		} else {
			EXIT_ERROR("Unsupported: %d %d\n",test_case->has_1st_order,test_case->has_2nd_order);
		}
		break;
	default:
		EXIT_ERROR("Unsupported: %c\n",test_case->solver_method_curr);
		break;
	}

	return s_params;
}

void constructor_Numerical_Flux_Input_data
	(struct Numerical_Flux_Input* num_flux_i, const struct Face* face, const struct Simulation* sim)
{
	struct Test_Case* test_case = sim->test_case;
	struct Solver_Face* s_face       = (struct Solver_Face*) face;
	struct DG_Solver_Face* dg_s_face = (struct DG_Solver_Face*) face;

	test_case->constructor_Boundary_Value_Input_face_fcl(&num_flux_i->bv_l,s_face,sim);        // destructed
	dg_s_face->constructor_Boundary_Value_fcl(&num_flux_i->bv_r,&num_flux_i->bv_l,s_face,sim); // destructed
}

void destructor_Numerical_Flux_Input_data (struct Numerical_Flux_Input* num_flux_i)
{
	destructor_Boundary_Value_Input(&num_flux_i->bv_l);
	destructor_Boundary_Value(&num_flux_i->bv_r);
}

// Level 1 ********************************************************************************************************** //

/// \brief Finalize the rhs term contribution from the \ref Face.
static void finalize_face_rhs_dg
	(const int side_index,                  ///< The index of the side of the face under consideration.
	 const struct Numerical_Flux* num_flux, ///< Defined for \ref compute_rlhs_fptr.
	 struct Face* face,                     ///< Defined for \ref compute_rlhs_fptr.
	 const struct Simulation* sim           ///< Defined for \ref compute_rlhs_fptr.
	);

/// \brief Finalize the lhs term contribution from the \ref Face.
static void finalize_face_lhs
	(const int side_index[2],                  /**< The indices of the affectee, affector, respectively. See the
	                                            *   comments in \ref compute_face_rlhs_dg.h for the convention. */
	 const struct Numerical_Flux* num_flux,    ///< Defined for \ref compute_rlhs_fptr.
	 struct Face* face,                        ///< Defined for \ref compute_rlhs_fptr.
	 struct Solver_Storage_Implicit* s_store_i ///< Defined for \ref compute_rlhs_fptr.
	);

static void scale_by_Jacobian_e (const struct Numerical_Flux* num_flux, struct Face* face, const struct Simulation* sim)
{
UNUSED(sim);
	struct Solver_Face* s_face = (struct Solver_Face*)face;

	const struct const_Vector_d jacobian_det_fc = interpret_const_Multiarray_as_Vector_d(s_face->jacobian_det_fc);
	scale_Multiarray_by_Vector_d('L',1.0,(struct Multiarray_d*)num_flux->nnf,&jacobian_det_fc,false);
}

static void scale_by_Jacobian_i1
	(const struct Numerical_Flux* num_flux, struct Face* face, const struct Simulation* sim)
{
UNUSED(sim);
	assert((!face->boundary && num_flux->neigh_info[1].dnnf_ds != NULL) ||
	       ( face->boundary && num_flux->neigh_info[1].dnnf_ds == NULL));

	struct Solver_Face* s_face = (struct Solver_Face*)face;

	const struct const_Vector_d jacobian_det_fc = interpret_const_Multiarray_as_Vector_d(s_face->jacobian_det_fc);
	scale_Multiarray_by_Vector_d('L',1.0,(struct Multiarray_d*)num_flux->nnf,&jacobian_det_fc,false);
	scale_Multiarray_by_Vector_d('L',1.0,(struct Multiarray_d*)num_flux->neigh_info[0].dnnf_ds,&jacobian_det_fc,false);
	if (!face->boundary)
		scale_Multiarray_by_Vector_d('L',1.0,(struct Multiarray_d*)num_flux->neigh_info[1].dnnf_ds,&jacobian_det_fc,false);
}

static void compute_rhs_f_dg
	(const struct Numerical_Flux* num_flux, struct Face* face, struct Solver_Storage_Implicit* s_store_i,
	 const struct Simulation* sim)
{
	UNUSED(s_store_i);
	assert(sim->elements->name == IL_ELEMENT_SOLVER_DG);

	finalize_face_rhs_dg(0,num_flux,face,sim);
	if (!face->boundary) {
		permute_Multiarray_d_fc((struct Multiarray_d*)num_flux->nnf,'R',1,(struct Solver_Face*)face);
// Can remove both of the scalings (here and and finalize_face_rhs_dg).
		scale_Multiarray_d((struct Multiarray_d*)num_flux->nnf,-1.0); // Use "-ve" normal.
		finalize_face_rhs_dg(1,num_flux,face,sim);
	}
}

static void compute_rlhs_1
	(const struct Numerical_Flux* num_flux, struct Face* face, struct Solver_Storage_Implicit* s_store_i,
	 const struct Simulation* sim)
{
	/// See \ref compute_face_rlhs_dg.h for the `lhs_**` notation.
	compute_rhs_f_dg(num_flux,face,s_store_i,sim);

	finalize_face_lhs((int[]){0,0},num_flux,face,s_store_i); // lhs_ll
	if (!face->boundary) {
		finalize_face_lhs((int[]){0,1},num_flux,face,s_store_i); // lhs_lr

		for (int i = 0; i < 2; ++i) {
			const struct Neigh_Info_NF* n_i = &num_flux->neigh_info[i];
			permute_Multiarray_d_fc((struct Multiarray_d*)n_i->dnnf_ds,'R',1,(struct Solver_Face*)face);
			scale_Multiarray_d((struct Multiarray_d*)n_i->dnnf_ds,-1.0); // Use "-ve" normal.
		}

		finalize_face_lhs((int[]){1,0},num_flux,face,s_store_i); // lhs_rl
		finalize_face_lhs((int[]){1,1},num_flux,face,s_store_i); // lhs_rr
	}
}

// Level 2 ********************************************************************************************************** //

static void finalize_face_rhs_dg
	(const int side_index, const struct Numerical_Flux* num_flux, struct Face* face, const struct Simulation* sim)
{
	const struct Operator* tw0_vs_fc = get_operator__tw0_vs_fc__rlhs_dg(side_index,face);

UNUSED(sim);
	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	const char op_format = 'd';

	struct DG_Solver_Volume* dg_s_vol = (struct DG_Solver_Volume*) face->neigh_info[side_index].volume;

//printf("%d\n",vol->index);
	mm_NNC_Operator_Multiarray_d(-1.0,1.0,tw0_vs_fc,num_flux->nnf,dg_s_vol->rhs,op_format,2,NULL,NULL);
//print_Multiarray_d(dg_s_vol->rhs);
}

static void finalize_face_lhs
	(const int side_index[2], const struct Numerical_Flux* num_flux, struct Face* face,
	 struct Solver_Storage_Implicit* s_store_i)
{
	struct Solver_Face* s_face     = (struct Solver_Face*) face;
	struct Volume* vol[2]          = { face->neigh_info[side_index[0]].volume,
	                                   face->neigh_info[side_index[1]].volume, };
	struct Solver_Volume* s_vol[2] = { (struct Solver_Volume*) vol[0],
	                                   (struct Solver_Volume*) vol[1], };

	const struct const_DG_Solver_Element* e[2] = { (struct const_DG_Solver_Element*) vol[0]->element,
	                                               (struct const_DG_Solver_Element*) vol[1]->element, };

	const int ind_lf[2]   = { face->neigh_info[side_index[0]].ind_lf,   face->neigh_info[side_index[1]].ind_lf, },
	          ind_href[2] = { face->neigh_info[side_index[0]].ind_href, face->neigh_info[side_index[1]].ind_href, };
	const int p_v[2] = { s_vol[0]->p_ref, s_vol[1]->p_ref, },
	          p_f = s_face->p_ref;

	const int curved = ( (s_face->cub_type == 's') ? 0 : 1 );
	const struct Operator* tw0_vs_fc =
		get_Multiarray_Operator(e[0]->tw0_vs_fc[curved],(ptrdiff_t[]){ind_lf[0],ind_href[0],0,p_f,p_v[0]});

	const struct Operator* cv0_vs_fc_op =
		get_Multiarray_Operator(e[1]->cv0_vs_fc[curved],(ptrdiff_t[]){ind_lf[1],ind_href[1],0,p_f,p_v[1]});

	const struct const_Matrix_d* cv0_vs_fc = cv0_vs_fc_op->op_std;
	bool need_free_cv0 = false;
	if (side_index[0] != side_index[1]) {
		need_free_cv0 = true;
		cv0_vs_fc = constructor_copy_const_Matrix_d(cv0_vs_fc); // destructed
		permute_Matrix_d_fc((struct Matrix_d*)cv0_vs_fc,'R',side_index[0],s_face);
	}


	const ptrdiff_t ext_0 = tw0_vs_fc->op_std->ext_0,
	                ext_1 = tw0_vs_fc->op_std->ext_1;

	struct Matrix_d* tw0_nf = constructor_empty_Matrix_d('R',ext_0,ext_1);            // destructed
	struct Matrix_d* lhs    = constructor_empty_Matrix_d('R',ext_0,cv0_vs_fc->ext_1); // destructed
	set_to_value_Matrix_d(tw0_nf,0.0);

	const struct const_Multiarray_d* dnnf_ds_Ma = num_flux->neigh_info[side_index[1]].dnnf_ds;
	struct Vector_d dnnf_ds = { .ext_0 = dnnf_ds_Ma->extents[0], .owns_data = false, .data = NULL, };

	const int n_eq = dnnf_ds_Ma->extents[1],
	          n_vr = dnnf_ds_Ma->extents[2];
	for (int vr = 0; vr < n_vr; ++vr) {
	for (int eq = 0; eq < n_eq; ++eq) {
		const ptrdiff_t ind =
			compute_index_sub_container(dnnf_ds_Ma->order,1,dnnf_ds_Ma->extents,(ptrdiff_t[]){eq,vr});
		dnnf_ds.data = (double*)&dnnf_ds_Ma->data[ind];
		mm_diag_d('R',1.0,0.0,tw0_vs_fc->op_std,(struct const_Vector_d*)&dnnf_ds,tw0_nf,false);

		mm_d('N','N',-1.0,0.0,(struct const_Matrix_d*)tw0_nf,cv0_vs_fc,lhs);

//printf("var, eq: %d %d\n",vr,eq);
//print_Matrix_d(lhs);
		set_petsc_Mat_row_col(s_store_i,s_vol[side_index[1]],eq,s_vol[side_index[0]],vr);
		add_to_petsc_Mat(s_store_i,(struct const_Matrix_d*)lhs);
	}}
	destructor_Matrix_d(tw0_nf);
	destructor_Matrix_d(lhs);

	if (need_free_cv0)
		destructor_const_Matrix_d(cv0_vs_fc);
}
