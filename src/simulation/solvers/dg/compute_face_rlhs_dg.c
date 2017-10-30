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
 *  \param num_flux \ref Numerical_Flux.
 *  \param face     \ref Face.
 *  \param sim      \ref Simulation.
 */
typedef void (*compute_rlhs_fptr)
	(const struct Numerical_Flux* num_flux,
	 struct Face* face,
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
UNUSED(s_store_i);
	assert(sim->faces->name == IL_FACE_SOLVER_DG);

	struct S_Params s_params = set_s_params(sim);
	struct Numerical_Flux_Input* num_flux_i = constructor_Numerical_Flux_Input(sim); // destructed

	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Face* face = (struct Face*) curr;

		constructor_Numerical_Flux_Input_data(num_flux_i,face,sim); // destructed
/// \todo clean-up
		// Compute the normal numerical fluxes (and optionally their Jacobians) at the face cubature nodes.
		struct Numerical_Flux* num_flux = constructor_Numerical_Flux(num_flux_i); // destructed
		destructor_Numerical_Flux_Input_data(num_flux_i);
//print_const_Multiarray_d(num_flux->nnf);

		s_params.scale_by_Jacobian(num_flux,face,sim);
//print_const_Multiarray_d(num_flux->nnf);

		s_params.compute_rlhs(num_flux,face,sim);
		destructor_Numerical_Flux(num_flux);
	}
	destructor_Numerical_Flux_Input(num_flux_i);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Scale \ref Numerical_Flux::nnf by the face Jacobian (i.e. only the explicit term).
static void scale_by_Jacobian_e
	(const struct Numerical_Flux* num_flux, ///< Defined for \ref scale_by_Jacobian_fptr.
	 struct Face* face,                     ///< Defined for \ref scale_by_Jacobian_fptr.
	 const struct Simulation* sim           ///< Defined for \ref scale_by_Jacobian_fptr.
	);

/// \brief Compute only the rhs term.
static void compute_rhs
	(const struct Numerical_Flux* num_flux, ///< Defined for \ref compute_rlhs_fptr.
	 struct Face* face,                     ///< Defined for \ref compute_rlhs_fptr.
	 const struct Simulation* sim           ///< Defined for \ref compute_rlhs_fptr.
	);

static struct S_Params set_s_params (const struct Simulation* sim)
{
	struct S_Params s_params;

	struct Test_Case* test_case = sim->test_case;
	switch (test_case->solver_method_curr) {
	case 'e':
		s_params.scale_by_Jacobian = scale_by_Jacobian_e;
		s_params.compute_rlhs = compute_rhs;
		break;
	case 'i':
		if (test_case->has_1st_order && !test_case->has_2nd_order)
			EXIT_ADD_SUPPORT; //s_params.compute_rlhs = compute_rlhs_1;
		else if (!test_case->has_1st_order && test_case->has_2nd_order)
			EXIT_ADD_SUPPORT; //s_params.compute_rlhs = compute_rlhs_2;
		else if (test_case->has_1st_order && test_case->has_2nd_order)
			EXIT_ADD_SUPPORT; //s_params.compute_rlhs = compute_rlhs_12;
		else
			EXIT_ERROR("Unsupported: %d %d\n",test_case->has_1st_order,test_case->has_2nd_order);
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

	test_case->constructor_Boundary_Value_Input_face_fcl(&num_flux_i->bv_l,s_face,sim); // destructed

	dg_s_face->constructor_Boundary_Value_fcl(
		&num_flux_i->bv_r,(const struct Boundary_Value_Input*)&num_flux_i->bv_l,s_face,sim); // destructed
}

void destructor_Numerical_Flux_Input_data (struct Numerical_Flux_Input* num_flux_i)
{
	destructor_Numerical_Flux_Input_mem(num_flux_i);
}

// Level 1 ********************************************************************************************************** //

/// \brief Finalize the rhs term contribution from the \ref Face.
static void finalize_face_rhs
	(const int side_index,                  ///< The index of the side of the face under consideration.
	 const struct Numerical_Flux* num_flux, ///< Defined for \ref compute_rlhs_fptr.
	 struct Face* face,                     ///< Defined for \ref compute_rlhs_fptr.
	 const struct Simulation* sim           ///< Defined for \ref compute_rlhs_fptr.
	);

static void scale_by_Jacobian_e (const struct Numerical_Flux* num_flux, struct Face* face, const struct Simulation* sim)
{
UNUSED(sim);
	struct Solver_Face* s_face = (struct Solver_Face*)face;

	const struct const_Vector_d jacobian_det_fc = interpret_const_Multiarray_as_Vector_d(s_face->jacobian_det_fc);
	scale_Multiarray_by_Vector_d('L',1.0,(struct Multiarray_d*)num_flux->nnf,&jacobian_det_fc,false);
}

static void compute_rhs (const struct Numerical_Flux* num_flux, struct Face* face, const struct Simulation* sim)
{
	assert(sim->elements->name == IL_ELEMENT_SOLVER_DG);

	finalize_face_rhs(0,num_flux,face,sim);
	if (!face->boundary) {
		permute_Multiarray_d_fc((struct Multiarray_d*)num_flux->nnf,'R',1,(struct Solver_Face*)face);
// Can remove both of the scalings (here and and finalize_face_rhs).
		scale_Multiarray_d((struct Multiarray_d*)num_flux->nnf,-1.0); // Use "-ve" normal.
		finalize_face_rhs(1,num_flux,face,sim);
	}
}

// Level 2 ********************************************************************************************************** //

static void finalize_face_rhs
	(const int side_index, const struct Numerical_Flux* num_flux, struct Face* face, const struct Simulation* sim)
{
UNUSED(sim);
	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	const char op_format = 'd';

	struct Solver_Face* s_face        = (struct Solver_Face*) face;
	struct Volume* vol                = (struct Volume*) face->neigh_info[side_index].volume;
	struct DG_Solver_Volume* dg_s_vol = (struct DG_Solver_Volume*) vol;

	const struct DG_Solver_Element* e = (const struct DG_Solver_Element*) vol->element;

	const int ind_lf   = face->neigh_info[side_index].ind_lf,
	          ind_href = face->neigh_info[side_index].ind_href;
	const int p_v = ((struct Solver_Volume*)vol)->p_ref,
	          p_f = s_face->p_ref;

	const struct Operator* tw0_vs_fc = ( (s_face->cub_type == 's')
		? get_Multiarray_Operator(e->tw0_vs_fcs,(ptrdiff_t[]){ind_lf,ind_href,0,p_f,p_v})
		: get_Multiarray_Operator(e->tw0_vs_fcc,(ptrdiff_t[]){ind_lf,ind_href,0,p_f,p_v}) );

	mm_NNC_Operator_Multiarray_d(-1.0,1.0,tw0_vs_fc,num_flux->nnf,dg_s_vol->rhs,op_format,2,NULL,NULL);
}
