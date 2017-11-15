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

#include "test_complex_compute_all_rhs_dpg.h"

#include <assert.h>

#include "macros.h"
#include "definitions_intrusive.h"

#include "test_complex_flux.h"
#include "test_complex_operators.h"
#include "test_support_math_functions.h"


#include "complex_multiarray.h"
#include "multiarray.h"

#include "intrusive.h"
#include "multiarray_operator.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void compute_all_rhs_dpg_c
	(const struct Complex_DPG_Solver_Volume* c_dpg_s_vol, struct Solver_Storage_Implicit* ssi,
	 const struct Simulation* sim)
{
	assert(sim->volumes->name == IL_VOLUME_SOLVER_DPG_COMPLEX);
	assert(sim->elements->name == IL_ELEMENT_SOLVER_DPG);

	struct Flux_Input_c* flux_i = constructor_Flux_Input_c(sim); // destructed

	destructor_Flux_Input_c(flux_i);
UNUSED(c_dpg_s_vol);
UNUSED(ssi);
EXIT_ADD_SUPPORT;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
