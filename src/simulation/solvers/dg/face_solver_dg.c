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

#include "face_solver_dg.h"

#include <assert.h>
#include <string.h>

#include "macros.h"
#include "definitions_test_case.h"
#include "definitions_bc.h"

#include "face.h"
#include "volume_solver_dg.h"

#include "multiarray.h"

#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Set the function pointers to the appropriate functions to compute values needed for the numerical flux
 *         computation. */
static void set_function_pointers_num_flux
	(struct Face* face_ptr,       ///< Pointer to the \ref Face.
	 const struct Simulation* sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void constructor_derived_DG_Solver_Face (struct Face* face_ptr, const struct Simulation* sim)
{
//	struct DG_Solver_Face* face = (struct DG_Solver_Face*) face_ptr;

	set_function_pointers_num_flux(face_ptr,sim);
}

void destructor_derived_DG_Solver_Face (struct Face* face_ptr)
{
	UNUSED(face_ptr);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Set the function pointers to the appropriate functions to compute boundary values needed for the numerical
 *         flux computation. */
static void set_function_pointers_num_flux_bc
	(struct Face* face_ptr,       ///< Defined for \ref set_function_pointers_num_flux.
	 const struct Simulation* sim ///< Defined for \ref set_function_pointers_num_flux.
	);

static void set_function_pointers_num_flux (struct Face* face_ptr, const struct Simulation* sim)
{
	assert(sim->method == METHOD_DG); // can be made flexible in future.

	struct DG_Solver_Face* dg_s_face = (struct DG_Solver_Face*) face_ptr;
	if (!face_ptr->boundary) {
		struct Test_Case* test_case = sim->test_case;
		switch (test_case->pde_index) {
		case PDE_ADVECTION:
		case PDE_EULER:
			dg_s_face->constructor_Boundary_Value_fcl = constructor_Boundary_Value_s_fcl_interp;
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
		set_function_pointers_num_flux_bc(face_ptr,sim);
	}
}

// Level 1 ********************************************************************************************************** //

/** \brief Set the function pointers to the appropriate functions to compute boundary values needed for the numerical
 *         flux computation. */
static void set_function_pointers_num_flux_bc_advection
	(struct Face* face_ptr,       ///< Defined for \ref set_function_pointers_num_flux.
	 const struct Simulation* sim ///< Defined for \ref set_function_pointers_num_flux.
	);

static void set_function_pointers_num_flux_bc (struct Face* face_ptr, const struct Simulation* sim)
{
	switch (sim->test_case->pde_index) {
	case PDE_ADVECTION: set_function_pointers_num_flux_bc_advection(face_ptr,sim); break;
	default:
		EXIT_ERROR("Unsupported: %d\n",sim->test_case->pde_index);
		break;
	}
}

// Level 2 ********************************************************************************************************** //

#include "boundary_advection.h"

static void set_function_pointers_num_flux_bc_advection (struct Face* face_ptr, const struct Simulation* sim)
{
UNUSED(sim);
	struct DG_Solver_Face* dg_s_face = (struct DG_Solver_Face*) face_ptr;

	const int bc = face_ptr->bc % BC_STEP_SC;
	switch (bc) {
	case BC_INFLOW:
		dg_s_face->constructor_Boundary_Value_fcl = constructor_Boundary_Value_advection_inflow;
		break;
	case BC_OUTFLOW:
		dg_s_face->constructor_Boundary_Value_fcl = constructor_Boundary_Value_advection_outflow;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",face_ptr->bc);
		break;
	}
}
