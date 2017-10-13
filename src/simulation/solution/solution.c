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

#include "solution.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"
#include "definitions_intrusive.h"

#include "multiarray.h"

#include "computational_elements.h"
#include "solver_volume.h"

#include "simulation.h"

// Static function declarations ************************************************************************************* //

/** \brief Set up the initial solution for the simulation. Computes:
 *	- \ref Solver_Volume::sol_coef;
 *	- \ref Solver_Volume::grad_coef (if applicable).
 */
static void set_up_solution
	(struct Simulation* sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void compute_solution (struct Simulation* sim)
{
	assert(sim->volumes->name == IL_SOLVER_VOLUME);
	assert(sim->faces->name   == IL_SOLVER_FACE);

	constructor_derived_Elements(sim,IL_SOLVER_ELEMENT);

	set_up_solution(sim);

	destructor_derived_Elements(sim,IL_ELEMENT);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Compute the initial \ref Solver_Volume::sol_coef.
static void compute_sol_coef
	(const struct Simulation*const sim, ///< \ref Simulation.
	 struct Solver_Volume*const volume  ///< \ref Volume.
	);

static void set_up_solution (struct Simulation* sim)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume* volume = (struct Solver_Volume*) curr;

		compute_sol_coef(sim,volume);

printf("sol: %d\n",((struct Volume*)volume)->index);
print_Multiarray_d(volume->sol_coef);
	}

EXIT_UNSUPPORTED;
}

// Level 1 ********************************************************************************************************** //

static void compute_sol_coef (const struct Simulation*const sim, struct Solver_Volume*const volume)
{
UNUSED(sim);
UNUSED(volume);
EXIT_UNSUPPORTED;
	// Find geometry coefficients at cubature nodes
	// Compute the mass matrix and find the L2 projection as inv(mass)*phi_vc'*w_v*j_v*sol_vc
	// Convert to coefficients

	// Compute solution at volume cubature nodes
}
