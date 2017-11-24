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

#include "solution_advection.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "macros.h"
#include "definitions_core.h"
#include "definitions_test_case.h"

#include "multiarray.h"

#include "boundary.h"
#include "compute_error.h"
#include "compute_error_advection.h"
#include "const_cast.h"
#include "file_processing.h"
#include "flux_advection.h"
#include "numerical_flux_advection.h"
#include "simulation.h"
#include "solution.h"
#include "solution_advection_default.h"
#include "peterson/solution_peterson.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Set the function pointers to the numerical flux computing functions.
void set_function_pointers_num_flux
	(struct Test_Case* test_case,      ///< \ref Test_Case.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void set_function_pointers_solution_advection (struct Test_Case* test_case, const struct Simulation*const sim)
{
	test_case->set_grad = set_sg_do_nothing;
	if (strstr(sim->pde_spec,"peterson")) {
		test_case->constructor_sol = constructor_const_sol_peterson;
		test_case->set_sol         = set_sol_peterson;
		test_case->compute_source_rhs = compute_source_rhs_do_nothing;
		test_case->constructor_Error_CE = constructor_Error_CE_advection_all;
	} else if ((strcmp(sim->pde_spec,"demkowicz_dpg_ii") == 0) ||
	           (strcmp(sim->pde_spec,"steady/default") == 0)) {
		test_case->constructor_sol = constructor_const_sol_advection_default;
		test_case->set_sol         = set_sol_advection_default;
		test_case->compute_source_rhs = compute_source_rhs_advection_default;
		test_case->constructor_Error_CE = constructor_Error_CE_advection_all;
	} else {
		EXIT_ERROR("Unsupported: %s\n",sim->pde_spec);
	}

	const bool* flux_comp_mem_e = (bool[]){1,0,0},
	          * flux_comp_mem_i = (bool[]){1,1,0};
	for (int i = 0; i < MAX_NUM_FLUX_OUT; ++i) {
		const_cast_b(&test_case->flux_comp_mem_e[i],flux_comp_mem_e[i]);
		const_cast_b(&test_case->flux_comp_mem_i[i],flux_comp_mem_i[i]);
	}

	test_case->compute_Flux = compute_Flux_1;
	test_case->compute_Flux_e[0] = compute_Flux_advection;
	test_case->compute_Flux_i[0] = compute_Flux_advection_jacobian;

	set_function_pointers_num_flux(test_case,sim);

	test_case->constructor_Boundary_Value_Input_face_fcl = constructor_Boundary_Value_Input_face_s_fcl_interp;
}

struct Sol_Data__Advection get_sol_data_advection (const struct Simulation* sim)
{
	static bool need_input = true;

	static struct Sol_Data__Advection sol_data;
	if (need_input) {
		need_input = false;
		read_data_advection(sim->input_path,&sol_data);
	}

	return sol_data;
}

void read_data_advection (const char*const input_path, struct Sol_Data__Advection*const sol_data)
{
	const int count_to_find = 1;

	FILE* input_file = fopen_input(input_path,'s'); // closed

	int count_found = 0;
	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),input_file)) {
		if (strstr(line,"b_adv")) {
			read_skip_d_1(line,1,sol_data->b_adv,DMAX);
			++count_found;
		}
	}
	fclose(input_file);

	if (count_found != count_to_find)
		EXIT_ERROR("Did not find the required number of variables");
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

void set_function_pointers_num_flux (struct Test_Case* test_case, const struct Simulation*const sim)
{
	switch (sim->method) {
	case METHOD_DG: // fallthrough
	case METHOD_DPG:
		test_case->compute_Numerical_Flux = compute_Numerical_Flux_1;
		switch (test_case->ind_num_flux[0]) {
		case NUM_FLUX_UPWIND:
			test_case->compute_Numerical_Flux_e[0] = compute_Numerical_Flux_advection_upwind;
			test_case->compute_Numerical_Flux_i[0] = compute_Numerical_Flux_advection_upwind_jacobian;
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
