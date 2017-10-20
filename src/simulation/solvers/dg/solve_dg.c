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

#include "solve_dg.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"

#include "volume_solver_dg.h"

#include "multiarray.h"

#include "compute_grad_coef_dg.h"
#include "compute_volume_rlhs_dg.h"
#include "compute_face_rlhs_dg.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

/// \brief Set the memory of the rhs and lhs (if applicable) terms to zero for the volumes.
static void zero_memory_volumes
	(const struct Simulation* sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

double compute_rhs_dg (const struct Simulation* sim)
{
/// \todo Add assertions relevant to rhs?
	double max_rhs = 0.0;

	zero_memory_volumes(sim);
	compute_grad_coef_dg(sim);
	compute_volume_rlhs_dg(sim);
	compute_face_rlhs_dg(sim);
EXIT_UNSUPPORTED;

	return max_rhs;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void zero_memory_volumes (const struct Simulation* sim)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		set_to_value_Multiarray_d(((struct DG_Solver_Volume*)curr)->rhs,0.0);
	}
}
