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


#include "def_templates_compute_face_rlhs_dg.h"

#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"

#include "def_templates_face_solver.h"
#include "def_templates_face_solver_dg.h"
#include "def_templates_volume_solver.h"
#include "def_templates_volume_solver_dg.h"

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
 *  \param face     \ref Face.
 *  \param sim      \ref Simulation.
 */
typedef void (*scale_by_Jacobian_fptr_T)
	(const struct Numerical_Flux_T*const num_flux,
	 struct Face*const face,
	 const struct Simulation*const sim
	);

/** \brief Function pointer to the function used to evaluate the rhs (and optionally lhs) terms.
 *
 *  \param num_flux  \ref Numerical_Flux_T.
 *  \param dg_s_face \ref DG_Solver_Face_T.
 *  \param s_store_i \ref Solver_Storage_Implicit.
 *  \param sim       \ref Simulation.
 */
typedef void (*compute_rlhs_fptr_T)
	(const struct Numerical_Flux_T* num_flux,
	 struct DG_Solver_Face_T* dg_s_face,
	 struct Solver_Storage_Implicit* s_store_i,
	 const struct Simulation* sim
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

// Interface functions ********************************************************************************************** //

void compute_face_rlhs_dg_T
	(const struct Simulation* sim, struct Solver_Storage_Implicit* ssi, struct Intrusive_List* faces)
{
	assert(sim->elements->name == IL_ELEMENT_SOLVER_DG);
	assert(sim->faces->name    == IL_FACE_SOLVER_DG);
	assert(sim->volumes->name  == IL_VOLUME_SOLVER_DG);

	struct S_Params_T s_params = set_s_params_T(sim);
	struct Numerical_Flux_Input_T* num_flux_i = constructor_Numerical_Flux_Input_T(sim); // destructed

	for (struct Intrusive_Link* curr = faces->first; curr; curr = curr->next) {
		struct Face* face                  = (struct Face*) curr;
		struct Solver_Face_T* s_face       = (struct Solver_Face_T*) curr;
		struct DG_Solver_Face_T* dg_s_face = (struct DG_Solver_Face_T*) curr;
//printf("face: %d\n",face->index);

		constructor_Numerical_Flux_Input_data_T(num_flux_i,s_face,sim); // destructed
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

		s_params.scale_by_Jacobian(num_flux,face,sim);
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
	assert(list_is_derived_from("solver",'v',sim));
	assert(list_is_derived_from("solver",'f',sim));
	assert(list_is_derived_from("solver",'e',sim));

	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	test_case->solver_method_curr = 'e';

	struct S_Params_T s_params = set_s_params_T(sim);
	struct Numerical_Flux_Input_T* num_flux_i = constructor_Numerical_Flux_Input_T(sim); // destructed

	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Face* face            = (struct Face*) curr;
		struct Solver_Face_T* s_face = (struct Solver_Face_T*) curr;

		constructor_Numerical_Flux_Input_data_T(num_flux_i,s_face,sim); // destructed

		struct Numerical_Flux_T* num_flux = constructor_Numerical_Flux_T(num_flux_i); // destructed
		destructor_Numerical_Flux_Input_data_T(num_flux_i);

		s_params.scale_by_Jacobian(num_flux,face,sim);

		add_to_flux_imbalance(num_flux,s_face,sim);
		destructor_Numerical_Flux_T(num_flux);
	}
	destructor_Numerical_Flux_Input_T(num_flux_i);

	test_case->solver_method_curr = 0;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Scale \ref Numerical_Flux_T::nnf by the face Jacobian (i.e. only the explicit term).
static void scale_by_Jacobian_e_T
	(const struct Numerical_Flux_T* num_flux, ///< Defined for \ref scale_by_Jacobian_fptr_T.
	 struct Face* face,                       ///< Defined for \ref scale_by_Jacobian_fptr_T.
	 const struct Simulation* sim             ///< Defined for \ref scale_by_Jacobian_fptr_T.
	);

/// \brief Compute only the rhs term.
static void compute_rhs_f_dg_T
	(const struct Numerical_Flux_T* num_flux,   ///< Defined for \ref compute_rlhs_fptr_T.
	 struct DG_Solver_Face_T* dg_s_face,        ///< Defined for \ref compute_rlhs_fptr_T.
	 struct Solver_Storage_Implicit* s_store_i, ///< Defined for \ref compute_rlhs_fptr_T.
	 const struct Simulation* sim               ///< Defined for \ref compute_rlhs_fptr_T.
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
			EXIT_ADD_SUPPORT; //s_params.compute_rlhs = compute_rlhs_12;
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

// Level 1 ********************************************************************************************************** //

/// \brief Finalize the rhs term contribution from the \ref Face.
static void finalize_face_rhs_dg_T
	(const int side_index,                    ///< The index of the side of the face under consideration.
	 const struct Numerical_Flux_T* num_flux, ///< Defined for \ref compute_rlhs_fptr_T.
	 struct DG_Solver_Face_T* dg_s_face,      ///< Defined for \ref compute_rlhs_fptr_T.
	 const struct Simulation* sim             ///< Defined for \ref compute_rlhs_fptr_T.
	);

static void scale_by_Jacobian_e_T (const struct Numerical_Flux_T* num_flux, struct Face* face, const struct Simulation* sim)
{
UNUSED(sim);
	struct Solver_Face_T* s_face = (struct Solver_Face_T*)face;

	const struct const_Vector_R jacobian_det_fc = interpret_const_Multiarray_as_Vector_d(s_face->jacobian_det_fc);
	scale_Multiarray_T_by_Vector_R('L',1.0,(struct Multiarray_T*)num_flux->nnf,&jacobian_det_fc,false);
}

static void compute_rhs_f_dg_T
	(const struct Numerical_Flux_T* num_flux, struct DG_Solver_Face_T* dg_s_face,
	 struct Solver_Storage_Implicit* s_store_i, const struct Simulation* sim)
{
	UNUSED(s_store_i);
	assert(sim->elements->name == IL_ELEMENT_SOLVER_DG);

	const struct Face* face = (struct Face*) dg_s_face;
	finalize_face_rhs_dg_T(0,num_flux,dg_s_face,sim);
	if (!face->boundary) {
		permute_Multiarray_T_fc((struct Multiarray_T*)num_flux->nnf,'R',1,(struct Solver_Face_T*)face);
// Can remove both of the scalings (here and and finalize_face_rhs_dg_T).
		scale_Multiarray_T((struct Multiarray_T*)num_flux->nnf,-1.0); // Use "-ve" normal.
		finalize_face_rhs_dg_T(1,num_flux,dg_s_face,sim);
	}
}

// Level 2 ********************************************************************************************************** //

static void finalize_face_rhs_dg_T
	(const int side_index, const struct Numerical_Flux_T* num_flux, struct DG_Solver_Face_T* dg_s_face,
	 const struct Simulation* sim)
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
//print_Multiarray_T(dg_s_vol->rhs);
}
