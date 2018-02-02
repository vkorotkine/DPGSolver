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


#include "def_templates_solution_advection.h"

#include "def_templates_boundary.h"
#include "def_templates_flux.h"
#include "def_templates_numerical_flux.h"
#include "def_templates_solution.h"
#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Set the function pointers to the numerical flux computing functions.
static void set_function_pointers_num_flux_T
	(struct Test_Case_T* test_case,    ///< \ref Test_Case_T.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void set_function_pointers_solution_advection_T (struct Test_Case_T* test_case, const struct Simulation*const sim)
{
	test_case->set_grad = set_sg_do_nothing_T;
	if (strstr(sim->pde_spec,"peterson")) {
		test_case->constructor_sol              = constructor_const_sol_peterson_T;
		test_case->set_sol                      = set_sol_peterson_T;
		test_case->compute_source_rhs           = compute_source_rhs_do_nothing_T;
		test_case->add_to_flux_imbalance_source = add_to_flux_imbalance_source_do_nothing_T;
		test_case->constructor_Error_CE         = constructor_Error_CE_advection_all;
	} else if ((strcmp(sim->pde_spec,"demkowicz_dpg_ii") == 0) ||
	           (strcmp(sim->pde_spec,"steady/default") == 0)) {
		test_case->constructor_sol              = constructor_const_sol_advection_default_T;
		test_case->set_sol                      = set_sol_advection_default_T;
		test_case->compute_source_rhs           = compute_source_rhs_advection_default_T;
		test_case->add_to_flux_imbalance_source = add_to_flux_imbalance_source_advection_default_T;
		test_case->constructor_Error_CE         = constructor_Error_CE_advection_all;
	} else {
		EXIT_ERROR("Unsupported: %s\n",sim->pde_spec);
	}

	test_case->compute_Flux = compute_Flux_1_T;
	test_case->compute_Flux_iv[0] = compute_Flux_T_advection;

	set_function_pointers_num_flux_T(test_case,sim);

	test_case->constructor_Boundary_Value_Input_face_fcl = constructor_Boundary_Value_Input_face_s_fcl_interp_T;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void set_function_pointers_num_flux_T (struct Test_Case_T* test_case, const struct Simulation*const sim)
{
	switch (sim->method) {
	case METHOD_DG: // fallthrough
	case METHOD_DPG:
		test_case->compute_Numerical_Flux = compute_Numerical_Flux_1_T;
		switch (test_case->ind_num_flux[0]) {
		case NUM_FLUX_UPWIND:
			test_case->compute_Numerical_Flux_e[0] = compute_Numerical_Flux_T_advection_upwind;
			test_case->compute_Numerical_Flux_i[0] = compute_Numerical_Flux_T_advection_upwind_jacobian;
			break;
		default:
			EXIT_ERROR("Unsupported: %d.\n",test_case->ind_num_flux[0]);
			break;
		}
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",sim->method);
		break;
	}
}
