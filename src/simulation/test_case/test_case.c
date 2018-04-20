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

#include "definitions_dpg.h"
#include "definitions_test_case.h"

#include "const_cast.h"
#include "file_processing.h"
#include "restart.h"
#include "simulation.h"
#include "solution_advection.h"
#include "solution_diffusion.h"
#include "solution_euler.h"
#include "solution_navier_stokes.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "test_case_T.c"

struct Test_Case_rc* constructor_Test_Case_rc_real (const struct Simulation* sim)
{
	struct Test_Case_rc* test_case_rc = malloc(sizeof *test_case_rc); // free

	const_cast_b(&test_case_rc->is_real,true);
	test_case_rc->tc = (void*)constructor_Test_Case(sim); // destructed

	return test_case_rc;
}

void destructor_Test_Case_rc_real (struct Test_Case_rc* test_case_rc)
{
	assert(test_case_rc->is_real);

	destructor_Test_Case(test_case_rc->tc);
	free(test_case_rc);
}

bool test_case_requires_positivity (const struct Test_Case*const test_case)
{
	switch (test_case->pde_index) {
	case PDE_ADVECTION: // fallthrough
	case PDE_DIFFUSION:
		return false;
		break;
	case PDE_EULER:         // fallthrough
	case PDE_NAVIER_STOKES:
		return true;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",test_case->pde_index);
		break;
	}
}

bool test_case_explicitly_enforces_conservation (const struct Simulation*const sim)
{
	switch (sim->method) {
	case METHOD_DG:
		break; // Do nothing.
	case METHOD_DPG: {
		const struct Test_Case*const test_case = (struct Test_Case*) sim->test_case_rc->tc;
		if (test_case->ind_conservation > CONSERVATION_NOT_ENFORCED)
			return true;
		break;
	} default:
		EXIT_ERROR("Unsupported: %d",sim->method);
		break;
	}
	return false;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

