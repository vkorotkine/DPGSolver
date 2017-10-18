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

#include "solve.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"
#include "definitions_test_case.h"
#include "definitions_intrusive.h"

#include "computational_elements.h"

#include "intrusive.h"
#include "simulation.h"
#include "solve_explicit.h"
#include "solve_implicit.h"
#include "solution.h"
#include "test_case.h"

#include "compute_grad_coef_dg.h"
#include "compute_volume_rlhs_dg.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void solve_for_solution (struct Simulation* sim)
{
	assert(sim->volumes->name == IL_SOLVER_VOLUME);
	assert(sim->faces->name   == IL_SOLVER_FACE);

	constructor_derived_Elements(sim,IL_SOLUTION_ELEMENT);
	set_initial_solution(sim);
	destructor_derived_Elements(sim,IL_ELEMENT);

	const struct Test_Case* test_case = sim->test_case;
	switch (test_case->solver_proc) {
	case SOLVER_E:
		solve_explicit(sim);
		break;
	case SOLVER_I:
		solve_implicit(sim);
		break;
	case SOLVER_EI:
		solve_explicit(sim);
		solve_implicit(sim);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",sim->test_case->solver_proc);
		break;
	}
EXIT_ADD_SUPPORT;
// derive DG Solver_Volumes/Faces and start solving. Include RES in the solver volume for simplicity.

}

double compute_rhs (const struct Simulation* sim)
{
	double max_rhs = 0.0;

	switch (sim->method) {
	case METHOD_DG:
		compute_grad_coef_dg(sim);
// name change here: rhs -> rlhs
		compute_volume_rhs_dg(sim);
// make face->rhs point to the same memory as volume->rhs.
EXIT_UNSUPPORTED;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",sim->method);
		break;
	}

	return max_rhs;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

