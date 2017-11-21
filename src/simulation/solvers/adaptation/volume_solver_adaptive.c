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

#include "volume_solver_adaptive.h"

#include "macros.h"

#include "simulation/simulation.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void constructor_derived_Adaptive_Solver_Volume (struct Volume* volume_ptr, const struct Simulation* sim)
{
	struct Adaptive_Solver_Volume* a_s_volume = (struct Adaptive_Solver_Volume*) volume_ptr;
	UNUSED(a_s_volume);
	UNUSED(sim);
}

void destructor_derived_Adaptive_Solver_Volume (struct Volume* volume_ptr)
{
	struct Adaptive_Solver_Volume* a_s_volume = (struct Adaptive_Solver_Volume*) volume_ptr;
	UNUSED(a_s_volume);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
