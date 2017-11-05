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
 *  \todo Attempt to template these functions.
 */

#include "test_complex_compute_face_rhs_dg.h"

#include <assert.h>

#include "macros.h"
#include "definitions_intrusive.h"

#include "test_complex_numerical_flux.h"

#include "intrusive.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void compute_face_rhs_dg_c (const struct Simulation* sim, struct Intrusive_List* faces)
{
	assert(sim->elements->name == IL_ELEMENT_SOLVER_DG);
	assert(sim->faces->name    == IL_FACE_SOLVER_DG_COMPLEX);
	assert(sim->volumes->name  == IL_VOLUME_SOLVER_DG_COMPLEX);

	struct Numerical_Flux_Input_c* num_flux_i = constructor_Numerical_Flux_Input_c(sim); // destructed

	for (struct Intrusive_Link* curr = faces->first; curr; curr = curr->next) {
UNUSED(curr);
EXIT_ADD_SUPPORT;
	}
	destructor_Numerical_Flux_Input_c(num_flux_i);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
