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

#include "volume_solver.h"

#include "multiarray.h"

#include "geometry.h"
#include "intrusive.h"
#include "simulation.h"
#include "solve_explicit.h"
#include "solve_implicit.h"
#include "solution.h"
#include "test_case.h"

#include "solve_dg.h"

// Static function declarations ************************************************************************************* //

/** \brief Compute the number of 'd'egrees 'o'f 'f'reedom in the volume computational elements.
 *  \return See brief. */
static ptrdiff_t compute_dof_volumes
	(const struct Simulation* sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void solve_for_solution (struct Simulation* sim)
{
	assert(sim->volumes->name == IL_SOLVER_VOLUME);
	assert(sim->faces->name   == IL_SOLVER_FACE);

	set_up_solver_geometry(sim);
	set_initial_solution(sim);

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
}

double compute_rhs (const struct Simulation* sim)
{
/// \todo Add assertions relevant to rhs.
	double max_rhs = 0.0;

	switch (sim->method) {
	case METHOD_DG:
		max_rhs = compute_rhs_dg(sim);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",sim->method);
		break;
	}

	return max_rhs;
}

double compute_rlhs (const struct Simulation* sim, struct Solver_Storage_Implicit* s_store_i)
{
	double max_rhs = 0.0;

	switch (sim->method) {
		case METHOD_DG: max_rhs = compute_rlhs_dg(sim,s_store_i);    break;
		default:        EXIT_ERROR("Unsupported: %d\n",sim->method); break;
	}

	return max_rhs;
}

ptrdiff_t compute_dof (const struct Simulation* sim)
{
	ptrdiff_t dof = 0;
	switch (sim->method) {
	case METHOD_DG:
		dof += compute_dof_volumes(sim);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",sim->method);
		break;
	}
	return dof;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static ptrdiff_t compute_dof_volumes (const struct Simulation* sim)
{
	ptrdiff_t dof = 0;
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next)
		dof += ((struct Solver_Volume*)curr)->sol_coef->extents[0];
	return dof;
}
