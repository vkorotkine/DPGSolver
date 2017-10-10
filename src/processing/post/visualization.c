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

#include "visualization.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"
#include "definitions_intrusive.h"
#include "definitions_visualization.h"

#include "computational_elements.h"
#include "element_plotting.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

/// \brief Output the visualization of the specified output in a format suitable for Paraview.
static void output_visualization_paraview
	(const struct Simulation* sim, ///< \ref Simulation.
	 const int vis_type            ///< The type of visualization. Options: see \ref definitions_visualization.h.
	);

// Interface functions ********************************************************************************************** //

void output_visualization (struct Simulation* sim, const int vis_type)
{
	assert(sim->volumes->name == IL_SOLVER_VOLUME);
	assert(sim->faces->name   == IL_SOLVER_FACE);

	constructor_derived_Elements(sim,IL_PLOTTING_ELEMENT);

	output_visualization_paraview(sim,vis_type);

/// \todo change input here to desired base list.
	destructor_derived_Elements(sim,IL_PLOTTING_ELEMENT);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void output_visualization_paraview (const struct Simulation* sim, const int vis_type)
{
UNUSED(sim);
	switch (vis_type) {
	default:
		EXIT_ERROR("Unsupported: %d\n",vis_type);
		break;
	}
}
