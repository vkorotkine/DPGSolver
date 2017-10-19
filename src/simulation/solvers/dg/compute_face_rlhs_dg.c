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

#include "face_solver_dg.h"
#include "element_solver_dg.h"

#include "simulation.h"

// Static function declarations ************************************************************************************* //

struct Num_Flux;

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
struct S_Params set_s_params
	(const struct Simulation* sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void compute_face_rlhs_dg (const struct Simulation* sim)
{
	assert(sim->faces->name == IL_FACE_SOLVER_DG);

	struct S_Params s_params = set_s_params(sim);

	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Face*        face   = (struct Face*) curr;
		struct Solver_Face* s_face = (struct Solver_Face*) curr;

		// Compute the solution and gradients at the face cubature nodes on each side of the face.
// place these immediately in the n_num_flux input container. Make sure that the input container has data for both a
// left and right volume.
		const struct const_Multiarray_d* left_s_fc = constructor_sol_fc(face,sim);
		// always interp

		const struct const_Multiarray_d* right_s_fc = constructor_sol_fc
		// either interp or boundary.

	}

	EXIT_ADD_SUPPORT;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

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

// Level 1 ********************************************************************************************************** //

static void compute_rhs (const struct Num_Flux* num_flux, struct Face* face, const struct Simulation* sim)
{
UNUSED(sim);
UNUSED(face);
UNUSED(num_flux);
EXIT_ADD_SUPPORT;
}
