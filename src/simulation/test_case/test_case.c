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

#include "test_case.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "macros.h"
#include "definitions_test_case.h"

#include "const_cast.h"
#include "file_processing.h"
#include "simulation.h"
#include "solution_euler.h"

// Static function declarations ************************************************************************************* //

/// \brief Set associations between `char*` and `int` variables.
static void set_string_associations
	(struct Test_Case* test_case,      ///< \ref Test_Case.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

/// \brief Set pde related parameters.
static void set_pde_related
	(struct Test_Case* test_case,      ///< \ref Test_Case.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

/// \brief Set the function pointer members of \ref Test_Case.
static void set_function_pointers
	(struct Test_Case* test_case,      ///< \ref Test_Case.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

/// \brief Read members of \ref Test_Case from input file.
static void read_test_case_parameters
	(struct Test_Case* test_case,      ///< \ref Test_Case.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

struct Test_Case* constructor_Test_Case (const struct Simulation* sim)
{
	struct Test_Case* test_case = calloc(1,sizeof *test_case); // returned

	set_string_associations(test_case,sim);
	set_pde_related(test_case,sim);

	read_test_case_parameters(test_case,sim);
	set_function_pointers(test_case,sim);

	test_case->solver_method_curr = 0;

	return test_case;
}

void destructor_Test_Case (const struct Test_Case* test_case)
{
	free((void*)test_case);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void set_string_associations (struct Test_Case* test_case, const struct Simulation*const sim)
{
	// pde_index
	if (strstr(sim->pde_name,"advection"))
		const_cast_i(&test_case->pde_index,PDE_ADVECTION);
	else if (strstr(sim->pde_name,"poisson"))
		const_cast_i(&test_case->pde_index,PDE_POISSON);
	else if (strstr(sim->pde_name,"euler"))
		const_cast_i(&test_case->pde_index,PDE_EULER);
	else if (strstr(sim->pde_name,"navier_stokes"))
		const_cast_i(&test_case->pde_index,PDE_NAVIER_STOKES);
	else
		EXIT_ERROR("Unsupported: %s\n",sim->pde_name);
}

static void set_pde_related (struct Test_Case* test_case, const struct Simulation* sim)
{
	switch (test_case->pde_index) {
	case PDE_ADVECTION:
		const_cast_i(&test_case->n_var,1);
		const_cast_i(&test_case->n_eq,1);
		const_cast_bool(&test_case->has_1st_order,true);
		const_cast_bool(&test_case->has_2nd_order,false);
		break;
	case PDE_POISSON:
		const_cast_i(&test_case->n_var,1);
		const_cast_i(&test_case->n_eq,1);
		const_cast_bool(&test_case->has_1st_order,false);
		const_cast_bool(&test_case->has_2nd_order,true);
		break;
	case PDE_EULER:
		const_cast_i(&test_case->n_var,sim->d+2);
		const_cast_i(&test_case->n_eq,sim->d+2);
		const_cast_bool(&test_case->has_1st_order,true);
		const_cast_bool(&test_case->has_2nd_order,false);
		break;
	case PDE_NAVIER_STOKES:
		const_cast_i(&test_case->n_var,sim->d+2);
		const_cast_i(&test_case->n_eq,sim->d+2);
		const_cast_bool(&test_case->has_1st_order,true);
		const_cast_bool(&test_case->has_2nd_order,true);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",test_case->pde_index);
		break;
	}
}

static void set_function_pointers (struct Test_Case* test_case, const struct Simulation*const sim)
{
	switch (test_case->pde_index) {
//		case PDE_ADVECTION:     set_function_pointers_solution_advection(test_case,sim);     break;
//		case PDE_POISSON:       set_function_pointers_solution_poisson(test_case,sim);       break;
		case PDE_EULER:
			set_function_pointers_solution_euler(test_case,sim);
//			test_case->flux_fptr =
			break;
//		case PDE_NAVIER_STOKES: set_function_pointers_solution_navier_stokes(test_case,sim); break;
		default: EXIT_ERROR("Unsupported: %d\n",test_case->pde_index); break;
	}
}

static void read_test_case_parameters (struct Test_Case* test_case, const struct Simulation*const sim)
{
	const int count_to_find = 1;

	FILE* input_file = fopen_input(sim->input_path,'t'); // closed

	int count_found = 0;
	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),input_file)) {
		if (strstr(line,"solver_proc")) {
			++count_found;
			read_skip_const_i(line,&test_case->solver_proc);
		}
		if (strstr(line,"solver_type_e")) read_skip_const_i(line,&test_case->solver_type_e);

		if (strstr(line,"time_final")) read_skip_const_d(line,&test_case->time_final,1,false);
		if (strstr(line,"time_step"))  read_skip_const_d(line,&test_case->dt,1,false);

		if (strstr(line,"display_progress")) read_skip_const_b(line,&test_case->display_progress);
	}
	fclose(input_file);

	if (count_found != count_to_find)
		EXIT_ERROR("Did not find the required number of variables");
}
