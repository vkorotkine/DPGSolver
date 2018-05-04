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


#include "def_templates_compute_face_rlhs_dg.h"

#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"

#include "def_templates_face_solver.h"
#include "def_templates_face_solver_dg.h"
#include "def_templates_volume_solver.h"
#include "def_templates_volume_solver_dg.h"

#include "def_templates_boundary.h"
#include "def_templates_compute_face_rlhs.h"
#include "def_templates_numerical_flux.h"
#include "def_templates_operators.h"
#include "def_templates_solve_dg.h"
#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

struct Num_Flux_T;

/** \brief Function pointer to the function used to scale by the face Jacobian.
 *
 *  \param num_flux \ref Numerical_Flux_T.
 *  \param s_face   \ref Solver_Face_T.
 */
typedef void (*scale_by_Jacobian_fptr_T)
	(struct Numerical_Flux_T*const num_flux,
	 const struct Solver_Face_T*const s_face
	);

/** \brief Function pointer to the function used to evaluate the rhs (and optionally lhs) terms.
 *
 *  \param num_flux  \ref Numerical_Flux_T.
 *  \param dg_s_face \ref DG_Solver_Face_T.
 *  \param s_store_i \ref Solver_Storage_Implicit.
 *  \param sim       \ref Simulation.
 */
typedef void (*compute_rlhs_fptr_T)
	(const struct Numerical_Flux_T*const num_flux,
	 struct DG_Solver_Face_T*const dg_s_face,
	 struct Solver_Storage_Implicit*const s_store_i,
	 const struct Simulation*const sim
	);

/// \brief Container for solver-related parameters.
struct S_Params_T {
	scale_by_Jacobian_fptr_T scale_by_Jacobian; ///< Pointer to the appropriate function.

	compute_rlhs_fptr_T compute_rlhs; ///< Pointer to the appropriate function.
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

/// \brief Add the current face contribution to \ref Solver_Volume_T::flux_imbalance from both sides.
static void add_to_flux_imbalance
	(const struct Numerical_Flux_T*const num_flux_w_J, ///< \ref Numerical_Flux_T scaled by face jacobian term.
	 const struct Solver_Face_T*const s_face,          ///< The current \ref Solver_Face_T.
	 const struct Simulation*const sim                 ///< \ref Simulation.
	);

/** \brief Constructor for \ref Boundary_Value_Input_T::g using interpolation of the partially corrected weak gradient
 *         to the 'f'ace 'c'ubature nodes as seen from the 'l'eft. */
static void constructor_Boundary_Value_Input_g_face_fcl
	(struct Boundary_Value_Input_T*const bv_i,       ///< \ref Boundary_Value_Input_T.
	 const struct DG_Solver_Face_T*const dg_s_face,  ///< \ref DG_Solver_Face_T.
	 const struct Simulation*const sim               ///< \ref Simulation.
	);

/** \brief Constructor for \ref Boundary_Value_T::g using interpolation of the partially corrected weak gradient to the
 *         to the 'f'ace 'c'ubature nodes as seen from the 'l'eft, **only if the face is not on a domain boundary**. */
static void constructor_Boundary_Value_g_face_fcl
	(struct Boundary_Value_T*const bv,              ///< \ref Boundary_Value_Input_T.
	 const struct DG_Solver_Face_T*const dg_s_face, ///< \ref DG_Solver_Face_T.
	 const struct Simulation*const sim              ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void compute_face_rlhs_dg_T
	(const struct Simulation* sim, struct Solver_Storage_Implicit* ssi, struct Intrusive_List* faces)
{
	assert(sim->elements->name == IL_ELEMENT_SOLVER_DG);
	assert(sim->faces->name    == IL_FACE_SOLVER_DG);
	assert(sim->volumes->name  == IL_VOLUME_SOLVER_DG);

	const struct Test_Case_T*const test_case = (struct Test_Case_T*) sim->test_case_rc->tc;
	const int has_2nd_order = test_case->has_2nd_order;

	struct S_Params_T s_params = set_s_params_T(sim);
	struct Numerical_Flux_Input_T* num_flux_i = constructor_Numerical_Flux_Input_T(sim); // destructed

	for (struct Intrusive_Link* curr = faces->first; curr; curr = curr->next) {
		struct Face* face                  = (struct Face*) curr;
		struct Solver_Face_T* s_face       = (struct Solver_Face_T*) curr;
		struct DG_Solver_Face_T* dg_s_face = (struct DG_Solver_Face_T*) curr;
UNUSED(face);
//printf("face: %d\n",face->index);

		constructor_Numerical_Flux_Input_data_dg_T(num_flux_i,dg_s_face,sim,has_2nd_order); // destructed
#if 0
print_const_Multiarray_T(num_flux_i->bv_l.s);
print_const_Multiarray_T(num_flux_i->bv_r.s);
#endif

		struct Numerical_Flux_T* num_flux = constructor_Numerical_Flux_T(num_flux_i); // destructed
		destructor_Numerical_Flux_Input_data_T(num_flux_i);
#if 0
print_const_Multiarray_T(num_flux->nnf);
print_const_Multiarray_T(num_flux->neigh_info[0].dnnf_ds);
if (!face->boundary)
	print_const_Multiarray_T(num_flux->neigh_info[1].dnnf_ds);
#endif

		s_params.scale_by_Jacobian(num_flux,s_face);
#if 0
print_const_Multiarray_T(num_flux->nnf);
print_const_Multiarray_T(num_flux->neigh_info[0].dnnf_ds);
if (!face->boundary)
	print_const_Multiarray_T(num_flux->neigh_info[1].dnnf_ds);
#endif
		s_params.compute_rlhs(num_flux,dg_s_face,ssi,sim);
		destructor_Numerical_Flux_T(num_flux);
//if (face->index == 2)
//break;
//EXIT_UNSUPPORTED;
	}
//EXIT_UNSUPPORTED;
	destructor_Numerical_Flux_Input_T(num_flux_i);
}

void compute_flux_imbalances_faces_dg_T (const struct Simulation*const sim)
{
return;
	assert(list_is_derived_from("solver",'v',sim));
	assert(list_is_derived_from("solver",'f',sim));
	assert(list_is_derived_from("solver",'e',sim));

	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	test_case->solver_method_curr = 'e';

	struct S_Params_T s_params = set_s_params_T(sim);
	struct Numerical_Flux_Input_T* num_flux_i = constructor_Numerical_Flux_Input_T(sim); // destructed

	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Solver_Face_T* s_face = (struct Solver_Face_T*) curr;

		constructor_Numerical_Flux_Input_data_T(num_flux_i,s_face,sim); // destructed

		struct Numerical_Flux_T* num_flux = constructor_Numerical_Flux_T(num_flux_i); // destructed
		destructor_Numerical_Flux_Input_data_T(num_flux_i);

		s_params.scale_by_Jacobian(num_flux,s_face);

		add_to_flux_imbalance(num_flux,s_face,sim);
		destructor_Numerical_Flux_T(num_flux);
	}
	destructor_Numerical_Flux_Input_T(num_flux_i);

	test_case->solver_method_curr = 0;
}

void constructor_Numerical_Flux_Input_data_dg_T
	(struct Numerical_Flux_Input_T*const num_flux_i, const struct DG_Solver_Face_T*const dg_s_face,
	 const struct Simulation*const sim, const bool compute_gradient)
{
	const struct Solver_Face_T*const s_face = (struct Solver_Face_T*) dg_s_face;

	if (compute_gradient) {
		constructor_Boundary_Value_Input_g_face_fcl(&num_flux_i->bv_l,dg_s_face,sim); // destructed
		constructor_Boundary_Value_g_face_fcl(&num_flux_i->bv_r,dg_s_face,sim);       // destructed
	}
	constructor_Numerical_Flux_Input_data_T(num_flux_i,s_face,sim); // destructed
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Scale \ref Numerical_Flux_T::nnf by the face Jacobian (i.e. only the explicit term).
static void scale_by_Jacobian_e_T
	(struct Numerical_Flux_T*const num_flux, ///< See brief.
	 const struct Solver_Face_T*const s_face ///< See brief.
	);

/// \brief Version of \ref compute_rlhs_fptr_T computing only the rhs term.
static void compute_rhs_f_dg_T
	(const struct Numerical_Flux_T*const num_flux,   ///< See brief.
	 struct DG_Solver_Face_T*const dg_s_face,        ///< See brief.
	 struct Solver_Storage_Implicit*const s_store_i, ///< See brief.
	 const struct Simulation*const sim               ///< See brief.
	);

/** \brief Constructor for the partially corrected weak gradient interpolated to the face cubature nodes as seen from
 *         the volume of input "side_index".
 *  \return See brief. */
static struct Multiarray_T* constructor_partial_grad_fc_interp
	(const int side_index,                          ///< The index of the side under consideration.
	 const struct DG_Solver_Face_T*const dg_s_face, ///< \ref DG_Solver_Face_T.
	 const struct Simulation*const sim              ///< \ref Simulation.
	);

static struct S_Params_T set_s_params_T (const struct Simulation* sim)
{
	struct S_Params_T s_params;

	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	switch (test_case->solver_method_curr) {
	case 'e':
		s_params.scale_by_Jacobian = scale_by_Jacobian_e_T;
		s_params.compute_rlhs      = compute_rhs_f_dg_T;
		break;
#if TYPE_RC == TYPE_REAL
	case 'i':
		if (test_case->has_1st_order && !test_case->has_2nd_order) {
			s_params.scale_by_Jacobian = scale_by_Jacobian_i1;
			s_params.compute_rlhs      = compute_rlhs_1;
		} else if (!test_case->has_1st_order && test_case->has_2nd_order) {
			s_params.scale_by_Jacobian = scale_by_Jacobian_i2;
			s_params.compute_rlhs      = compute_rlhs_2;
		} else if (test_case->has_1st_order && test_case->has_2nd_order) {
			s_params.scale_by_Jacobian = scale_by_Jacobian_i12;
			s_params.compute_rlhs      = compute_rlhs_12;
		} else {
			EXIT_ERROR("Unsupported: %d %d\n",test_case->has_1st_order,test_case->has_2nd_order);
		}
		break;
#endif
	default:
		EXIT_ERROR("Unsupported: %c (type_rc: %d)\n",test_case->solver_method_curr,TYPE_RC);
		break;
	}

	return s_params;
}

static void add_to_flux_imbalance
	(const struct Numerical_Flux_T*const num_flux_w_J, const struct Solver_Face_T*const s_face,
	 const struct Simulation*const sim)
{
	UNUSED(sim);
	const struct const_Matrix_T nnf_M = interpret_const_Multiarray_as_Matrix_T(num_flux_w_J->nnf);
	const struct const_Vector_R* w_fc = get_operator__w_fc__s_e_T(s_face);

	add_to_flux_imbalance_face_nf_w_T(&nnf_M,w_fc,s_face);
}

static void constructor_Boundary_Value_Input_g_face_fcl
	(struct Boundary_Value_Input_T*const bv_i, const struct DG_Solver_Face_T*const dg_s_face,
	 const struct Simulation*const sim)
{
	bv_i->g = (struct const_Multiarray_T*) constructor_partial_grad_fc_interp(0,dg_s_face,sim); // destructed
}

static void constructor_Boundary_Value_g_face_fcl
	(struct Boundary_Value_T*const bv, const struct DG_Solver_Face_T*const dg_s_face,
	 const struct Simulation*const sim)
{
	const struct Face*const face = (struct Face*) dg_s_face;
	if (face->boundary)
		return;

	const int side_index = 1;
	const struct Solver_Face_T*const s_face = (struct Solver_Face_T*) dg_s_face;

	struct Multiarray_T* grad_r_fcr = constructor_partial_grad_fc_interp(side_index,dg_s_face,sim); // moved
	permute_Multiarray_T_fc(grad_r_fcr,'R',side_index,s_face);

	bv->g = (struct const_Multiarray_T*) grad_r_fcr; // destructed
}

// Level 1 ********************************************************************************************************** //

/// \brief Finalize the rhs term contribution from the \ref Face.
static void finalize_face_rhs_dg_T
	(const int side_index,                         ///< The index of the side of the face under consideration.
	 const struct Numerical_Flux_T*const num_flux, ///< Defined for \ref compute_rlhs_fptr_T.
	 struct DG_Solver_Face_T*const dg_s_face,      ///< Defined for \ref compute_rlhs_fptr_T.
	 const struct Simulation*const sim             ///< Defined for \ref compute_rlhs_fptr_T.
	);

/** \brief Return the scaling for the face contribution to the weak gradient used to compute the numerical flux.
 *  \return See brief. */
static Real compute_scaling_weak_gradient
	(const struct DG_Solver_Face_T*const dg_s_face, ///< \ref DG_Solver_Face_T.
	 const struct Test_Case_T*const test_case       ///< \ref Test_Case_T.
	);

static void scale_by_Jacobian_e_T
	(struct Numerical_Flux_T*const num_flux, const struct Solver_Face_T*const s_face)
{
	const struct const_Vector_R jacobian_det_fc = interpret_const_Multiarray_as_Vector_d(s_face->jacobian_det_fc);
	scale_Multiarray_T_by_Vector_R('L',1.0,(struct Multiarray_T*)num_flux->nnf,&jacobian_det_fc,false);
}

static void compute_rhs_f_dg_T
	(const struct Numerical_Flux_T*const num_flux, struct DG_Solver_Face_T*const dg_s_face,
	 struct Solver_Storage_Implicit*const s_store_i, const struct Simulation*const sim)
{
	UNUSED(s_store_i);
	assert(sim->elements->name == IL_ELEMENT_SOLVER_DG);

	const struct Face* face = (struct Face*) dg_s_face;
	finalize_face_rhs_dg_T(0,num_flux,dg_s_face,sim);
	if (!face->boundary) {
		permute_Multiarray_T_fc((struct Multiarray_T*)num_flux->nnf,'R',1,(struct Solver_Face_T*)face);
		scale_Multiarray_T((struct Multiarray_T*)num_flux->nnf,-1.0); // Use "-ve" normal.
		finalize_face_rhs_dg_T(1,num_flux,dg_s_face,sim);
	}
}

static struct Multiarray_T* constructor_partial_grad_fc_interp
	(const int side_index, const struct DG_Solver_Face_T*const dg_s_face, const struct Simulation*const sim)
{

	const struct Face*const face              = (struct Face*) dg_s_face;
	const struct Solver_Face_T*const s_face   = (struct Solver_Face_T*) dg_s_face;
	const struct Volume*const vol             = (struct Volume*) face->neigh_info[side_index].volume;
	const struct DG_Solver_Volume_T* dg_s_vol = (struct DG_Solver_Volume_T*) vol;

	const struct Neigh_Info_DG*const ni = &dg_s_face->neigh_info[side_index];
	const struct Test_Case_T*const test_case = (struct Test_Case_T*) sim->test_case_rc->tc;

	const Real scale = compute_scaling_weak_gradient(dg_s_face,test_case);
	const struct Multiarray_T*const g_coef_p =
		constructor_sum_Multiarrays_Multiarray_T(1.0,dg_s_vol->grad_coef_v,scale,ni->grad_coef_f); // destructed


	const struct Operator*const cv0_vr_fc = get_operator__cv0_vr_fc_T(side_index,s_face);

	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	UNUSED(sim);
	const char op_format = 'd';

	struct Multiarray_T* ret =
		constructor_mm_NN1_Operator_Multiarray_T(cv0_vr_fc,g_coef_p,'C',op_format,g_coef_p->order,NULL); // returned
	destructor_Multiarray_T((struct Multiarray_T*)g_coef_p);
	return ret;
}

// Level 2 ********************************************************************************************************** //

/** Scaling factor for the penalty term used to compute the partially corrected weak gradient. Must be greater than 1
 *  for the scheme to be stable (Theorem 2, \cite Brdar2012). */
#define PENALTY_SCALING 1.01

static void finalize_face_rhs_dg_T
	(const int side_index, const struct Numerical_Flux_T*const num_flux, struct DG_Solver_Face_T*const dg_s_face,
	 const struct Simulation*const sim)
{
	const struct Face* face            = (struct Face*) dg_s_face;
	const struct Solver_Face_T* s_face = (struct Solver_Face_T*) face;

	const struct Operator* tw0_vt_fc = get_operator__tw0_vt_fc_T(side_index,s_face);

UNUSED(sim);
	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	const char op_format = 'd';

	struct DG_Solver_Volume_T* dg_s_vol = (struct DG_Solver_Volume_T*) face->neigh_info[side_index].volume;

//printf("%d\n",vol->index);
	mm_NNC_Operator_Multiarray_T(-1.0,1.0,tw0_vt_fc,num_flux->nnf,dg_s_vol->rhs,op_format,2,NULL,NULL);
#if 0
#if TYPE_RC == TYPE_REAL
const int bc = face->bc % BC_STEP_SC;
if (bc == BC_SLIPWALL) {
//print_const_Multiarray_T(num_flux->nnf);
//print_Multiarray_T(dg_s_vol->rhs);
//print_const_Matrix_R(tw0_vt_fc->op_std);
}
#endif
#endif
//print_Multiarray_T(dg_s_vol->rhs);
}

static Real compute_scaling_weak_gradient
	(const struct DG_Solver_Face_T*const dg_s_face, const struct Test_Case_T*const test_case)
{
	const int ind_num_flux_2nd = test_case->ind_num_flux[1];

	const struct Face*const face = (struct Face*) dg_s_face;
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

