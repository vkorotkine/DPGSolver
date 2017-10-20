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

#include <string.h>

#include "macros.h"
#include "definitions_test_case.h"

#include "face.h"
#include "volume_solver_dg.h"

#include "multiarray.h"

#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Set the function pointers to the appropriate functions to compute values needed for the numerical flux
 *         computation. */
static void set_function_pointers_num_flux
	(struct Face* face_ptr,       ///< Defined for \ref constructor_derived_DG_Solver_Face.
	 const struct Simulation* sim ///< Defined for \ref constructor_derived_DG_Solver_Face.
	);

// Interface functions ********************************************************************************************** //

void constructor_derived_DG_Solver_Face (struct Face* face_ptr, const struct Simulation* sim)
{
	struct DG_Solver_Face* face = (struct DG_Solver_Face*) face_ptr;

	for (int i = 0; i < 2; ++i) {
		struct DG_Solver_Volume* volume = (struct DG_Solver_Volume*) face_ptr->neigh_info[i].volume;
		face->rhs[i] = ( volume ? volume->rhs : NULL );
	}

	set_function_pointers_num_flux(face_ptr,sim);
}

void destructor_derived_DG_Solver_Face (struct Face* face_ptr)
{
	UNUSED(face_ptr);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void set_function_pointers_num_flux (struct Face* face_ptr, const struct Simulation* sim)
{
	struct DG_Solver_Face* dg_s_face = (struct DG_Solver_Face*) face_ptr;
	if (!face_ptr->boundary) {
// make external and move to another file after usage is set.
		struct Test_Case* test_case = sim->test_case;
		switch (test_case->pde_index) {
		case PDE_ADVECTION:
		case PDE_EULER:
			dg_s_face->constructor_s_r_fcl = constructor_s_r_fcl_interp;
			dg_s_face->constructor_g_r_fcl = constructor_sg_fc_null;
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
		EXIT_ADD_SUPPORT; // pointer to boundary condition function (depends only on bc index).
	}
}
