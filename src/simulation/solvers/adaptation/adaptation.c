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

#include "adaptation.h"

#include <assert.h>
#include <stdlib.h>

#include "macros.h"
#include "definitions_intrusive.h"

#include "simulation/simulation.h"
#include "simulation/computational_elements/computational_elements.h"
#include "simulation/solvers/adaptation/definitions_adaptation.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void adapt_hp (struct Simulation* sim, const int adapt_strategy)
{
	constructor_derived_computational_elements(sim,IL_SOLVER_ADAPTIVE); // destructed

	destructor_derived_computational_elements(sim,IL_SOLVER);
UNUSED(adapt_strategy);
EXIT_UNSUPPORTED;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
