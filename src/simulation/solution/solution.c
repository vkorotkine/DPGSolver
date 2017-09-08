// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "solution.h"

#include <stdlib.h>
#include <stdio.h>

#include "macros.h"

#include "multiarray.h"

#include "solver_volume.h"

// Static function declarations ************************************************************************************* //

/// \brief Compute the initial \ref Volume::sol_coef.
static void compute_sol_coef
	(const struct Simulation*const sim, ///< \ref Simulation.
	 struct Solver_Volume*const volume  ///< \ref Volume.
	);

// Interface functions ********************************************************************************************** //

void set_up_solution (struct Simulation* sim, struct Intrusive_List* solver_volumes)
{
	for (struct Intrusive_Link* curr = solver_volumes->first; curr; curr = curr->next) {
		struct Solver_Volume* volume = (struct Solver_Volume*) curr;

		compute_sol_coef(sim,volume);

printf("sol: %d\n",((struct Volume*)volume)->index);
print_Multiarray_d(volume->sol_coef,1e-10);
	}
EXIT_UNSUPPORTED;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

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
