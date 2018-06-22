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
#include <string.h>

#include "macros.h"
#include "definitions_core.h"
#include "definitions_numerical_flux.h"
#include "definitions_test_case.h"


#include "def_templates_solution_burgers_inviscid.h"

#include "def_templates_boundary.h"
#include "def_templates_flux.h"
#include "def_templates_geometry.h"
#include "def_templates_numerical_flux.h"
#include "def_templates_solution.h"
#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void set_function_pointers_solution_burgers_inviscid_T (struct Test_Case_T* test_case, const struct Simulation*const sim)
{
	test_case->set_grad = set_sg_do_nothing_T;
	if (strstr(sim->pde_spec,"periodic/trigonometric")) {
		test_case->constructor_xyz              = constructor_xyz_fixed_cube_parametric_T;
		test_case->constructor_sol              = constructor_const_sol_trigonometric_T;
		test_case->set_sol                      = set_sol_trigonometric_T;
		test_case->compute_source_rhs           = compute_source_rhs_do_nothing_T;
		test_case->add_to_flux_imbalance_source = add_to_flux_imbalance_source_do_nothing_T;
		test_case->constructor_Error_CE         = constructor_Error_CE_NULL;
		const_cast_b(&test_case->has_analytical,false);
		const_cast_b(&test_case->copy_initial_rhs,false);
	} else {
		EXIT_ERROR("Unsupported: %s\n",sim->pde_spec);
	}

	test_case->compute_Flux = compute_Flux_1_T;
	test_case->compute_Flux_iv[0] = compute_Flux_T_burgers_inviscid;
	test_case->compute_Flux_iv[1] = NULL;

	test_case->compute_Numerical_Flux = compute_Numerical_Flux_1_T;
	switch (test_case->ind_num_flux[0]) {
	case NUM_FLUX_LAX_FRIEDRICHS:
		test_case->compute_Numerical_Flux_e[0] = compute_Numerical_Flux_T_burgers_inviscid_lax_friedrichs;
		test_case->compute_Numerical_Flux_i[0] = compute_Numerical_Flux_T_burgers_inviscid_lax_friedrichs_jacobian;
		break;
	default:
		EXIT_ERROR("Unsupported: %d.\n",test_case->ind_num_flux[0]);
		break;
	}

	test_case->constructor_Boundary_Value_Input_face_fcl = constructor_Boundary_Value_Input_face_s_fcl_interp_T;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

#include "undef_templates_solution_burgers_inviscid.h"

#include "undef_templates_boundary.h"
#include "undef_templates_flux.h"
#include "undef_templates_geometry.h"
#include "undef_templates_numerical_flux.h"
#include "undef_templates_solution.h"
#include "undef_templates_test_case.h"
