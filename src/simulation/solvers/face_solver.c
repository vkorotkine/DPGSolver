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

#include "face_solver.h"

#include <assert.h>
#include <string.h>

#include "macros.h"
#include "definitions_test_case.h"
#include "definitions_bc.h"

#include "volume.h"

#include "multiarray.h"

#include "const_cast.h"
#include "geometry.h"
#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Checks if one of the neighbouring volumes to the current face is curved.
 *  \return `true` if curved volume is found; `false` otherwise. */
bool check_for_curved_neigh
	(struct Face* face ///< \ref Face.
	);

/** \brief Set the function pointers to the appropriate functions to compute values needed for the numerical flux
 *         computation. */
static void set_function_pointers_num_flux
	(struct Solver_Face* s_face,  ///< Pointer to the \ref Solver_Face.
	 const struct Simulation* sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void constructor_derived_Solver_Face (struct Face* face_ptr, const struct Simulation* sim)
{
	struct Solver_Face* s_face = (struct Solver_Face*) face_ptr;

	s_face->ind_dof = -1;
	const_cast_i(&s_face->p_ref,sim->p_ref[0]);
	const_cast_c(&s_face->cub_type,(check_for_curved_neigh((struct Face*)s_face) ? 'c' : 's'));

	s_face->nf_coef = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){0,0});   // destructed

	const_constructor_move_Multiarray_d(
		&s_face->xyz_fc,constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){0,0}));        // destructed
	const_constructor_move_Multiarray_d(
		&s_face->normals_fc,constructor_empty_Multiarray_d('R',2,(ptrdiff_t[]){0,0}));    // destructed
	const_constructor_move_Multiarray_d(
		&s_face->jacobian_det_fc,constructor_empty_Multiarray_d('C',1,(ptrdiff_t[]){0})); // destructed

	set_function_pointers_num_flux(s_face,sim);
}

void destructor_derived_Solver_Face (struct Face* face_ptr)
{
	struct Solver_Face* face = (struct Solver_Face*) face_ptr;

	destructor_Multiarray_d(face->nf_coef);

	destructor_const_Multiarray_d(face->xyz_fc);
	destructor_const_Multiarray_d(face->normals_fc);
	destructor_const_Multiarray_d(face->jacobian_det_fc);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Set the function pointers to the appropriate functions to compute boundary values needed for the numerical
 *         flux computation. */
static void set_function_pointers_num_flux_bc
	(struct Solver_Face* s_face,  ///< Defined for \ref set_function_pointers_num_flux.
	 const struct Simulation* sim ///< Defined for \ref set_function_pointers_num_flux.
	);

bool check_for_curved_neigh (struct Face* face)
{
	if (face->neigh_info[0].volume->curved || (face->neigh_info[1].volume && face->neigh_info[1].volume->curved))
		return true;
	return false;
}

static void set_function_pointers_num_flux (struct Solver_Face* s_face, const struct Simulation* sim)
{
	const struct Face* face = (struct Face*) s_face;
	if (!face->boundary) {
		struct Test_Case* test_case = sim->test_case;
		switch (test_case->pde_index) {
		case PDE_ADVECTION:
		case PDE_EULER:
			s_face->constructor_Boundary_Value_fcl = constructor_Boundary_Value_s_fcl_interp;
			break;
		case PDE_POISSON:
			EXIT_UNSUPPORTED;
			break;
		case PDE_NAVIER_STOKES:
			EXIT_UNSUPPORTED;
			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",test_case->pde_index);
			break;
		}
	} else {
		set_function_pointers_num_flux_bc(s_face,sim);
	}
}

// Level 1 ********************************************************************************************************** //

/** \brief Set the function pointers to the appropriate functions to compute boundary values needed for the numerical
 *         flux computation. */
static void set_function_pointers_num_flux_bc_advection
	(struct Solver_Face* s_face,  ///< Defined for \ref set_function_pointers_num_flux.
	 const struct Simulation* sim ///< Defined for \ref set_function_pointers_num_flux.
	);

static void set_function_pointers_num_flux_bc (struct Solver_Face* s_face, const struct Simulation* sim)
{
	switch (sim->test_case->pde_index) {
	case PDE_ADVECTION: set_function_pointers_num_flux_bc_advection(s_face,sim); break;
	default:
		EXIT_ERROR("Unsupported: %d\n",sim->test_case->pde_index);
		break;
	}
}

// Level 2 ********************************************************************************************************** //

#include "boundary_advection.h"

static void set_function_pointers_num_flux_bc_advection (struct Solver_Face* s_face, const struct Simulation* sim)
{
UNUSED(sim);
	const struct Face* face = (struct Face*) s_face;

	const int bc = face->bc % BC_STEP_SC;
	switch (bc) {
	case BC_INFLOW:
		s_face->constructor_Boundary_Value_fcl = constructor_Boundary_Value_advection_inflow;
		break;
	case BC_OUTFLOW:
		s_face->constructor_Boundary_Value_fcl = constructor_Boundary_Value_advection_outflow;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",face->bc);
		break;
	}
}
