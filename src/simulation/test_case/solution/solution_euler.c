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

#include "solution_euler.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "macros.h"

#include "simulation.h"
#include "solution.h"
#include "solution_periodic_vortex.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void set_function_pointers_solution_euler (struct Test_Case* test_case, const struct Simulation*const sim)
{
	test_case->compute_init_grad_coef_v = compute_grad_coef_v_do_nothing;
	test_case->compute_init_grad_coef_f = compute_grad_coef_f_do_nothing;
	if (strstr(sim->pde_spec,"periodic_vortex")) {
		test_case->compute_init_sol_coef_v = compute_sol_coef_v_periodic_vortex;
		switch (sim->method) {
		case METHOD_DG:
			test_case->compute_init_grad_coef_f = compute_grad_coef_f_do_nothing;
			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",sim->method);
			break;
		}
	} else {
		EXIT_ERROR("Unsupported: %s\n",sim->pde_spec);
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
