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

#include "multiarray.h"
#include "vector.h"

#include "face_solver_dg.h"
#include "element_solver_dg.h"

#include "simulation.h"
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
 *  \param num_flux \ref Num_Flux.
 *  \param face     \ref Face.
 *  \param sim      \ref Simulation.
 */
typedef void (*compute_rlhs_fptr)
	(const struct Num_Flux* num_flux,
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

void compute_face_rlhs_dg (const struct Simulation* sim)
{
	assert(sim->faces->name == IL_FACE_SOLVER_DG);

	struct S_Params s_params = set_s_params(sim);
UNUSED(s_params);
	struct Numerical_Flux_Input* num_flux_i = constructor_Numerical_Flux_Input(sim); // destructed

	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Face*        face         = (struct Face*) curr;
//		struct Solver_Face* s_face       = (struct Solver_Face*) curr;
//		struct DG_Solver_Face* dg_s_face = (struct DG_Solver_Face*) curr;

		constructor_Numerical_Flux_Input_data(num_flux_i,face,sim);

//print_const_Multiarray_d(num_flux_i->neigh_info[0].s);
//print_const_Multiarray_d(num_flux_i->neigh_info[1].s);

		// Compute the normal numerical fluxes (and optionally their Jacobians) at the face cubature nodes.
		struct Numerical_Flux* num_flux = constructor_Numerical_Flux(num_flux_i);
		destructor_Numerical_Flux_Input_data(num_flux_i);
print_const_Multiarray_d(num_flux->nnf);

		s_params.scale_by_Jacobian(num_flux,face,sim);
print_const_Multiarray_d(num_flux->nnf);

UNUSED(num_flux);



EXIT_UNSUPPORTED;
	}

	destructor_Numerical_Flux_Input(num_flux_i);

	EXIT_ADD_SUPPORT;
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
	(const struct Num_Flux* num_flux, ///< Defined for \ref compute_rlhs_fptr.
	 struct Face* face,                        ///< Defined for \ref compute_rlhs_fptr.
	 const struct Simulation* sim              ///< Defined for \ref compute_rlhs_fptr.
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

	num_flux_i->neigh_info[0].s = test_case->constructor_s_l_fcl(face,sim);
	num_flux_i->neigh_info[0].g = test_case->constructor_g_l_fcl(face,sim);
	num_flux_i->neigh_info[1].s = dg_s_face->constructor_s_r_fcl(face,sim);
	num_flux_i->neigh_info[1].g = dg_s_face->constructor_g_r_fcl(face,sim);

	num_flux_i->normals_l = s_face->normals_fc;
	num_flux_i->xyz_l     = s_face->xyz_fc;
}

void destructor_Numerical_Flux_Input_data (struct Numerical_Flux_Input* num_flux_i)
{
	for (int i = 0; i < 2; ++i) {
		if (num_flux_i->neigh_info[i].s)
			destructor_const_Multiarray_d(num_flux_i->neigh_info[i].s);
		if (num_flux_i->neigh_info[i].g)
			destructor_const_Multiarray_d(num_flux_i->neigh_info[i].g);
	}
}

// Level 1 ********************************************************************************************************** //

static void scale_by_Jacobian_e (const struct Numerical_Flux* num_flux, struct Face* face, const struct Simulation* sim)
{
UNUSED(sim);
	struct Solver_Face* s_face = (struct Solver_Face*)face;

	const struct const_Multiarray_d*const jacobian_det_fc_Ma = s_face->jacobian_det_fc;
	assert(jacobian_det_fc_Ma->order == 1);
	struct const_Vector_d jacobian_det_fc =
		{ .ext_0     = jacobian_det_fc_Ma->extents[0],
		  .owns_data = false,
		  .data      = jacobian_det_fc_Ma->data, };

	scale_Multiarray_by_Vector_d('L',1.0,(struct Multiarray_d*)num_flux->nnf,&jacobian_det_fc,false);
}

static void compute_rhs (const struct Num_Flux* num_flux, struct Face* face, const struct Simulation* sim)
{
UNUSED(sim);
UNUSED(face);
UNUSED(num_flux);
EXIT_ADD_SUPPORT;
}
