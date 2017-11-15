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

#include "volume_solver_dpg_complex.h"

#include "macros.h"

#include "volume.h"
#include "volume_solver.h"

#include "complex_multiarray_minimal.h"
#include "multiarray.h"

#include "simulation.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void constructor_derived_Complex_DPG_Solver_Volume (struct Volume* volume_ptr, const struct Simulation* sim)
{
	UNUSED(sim);
	const struct Solver_Volume* s_vol             = (struct Solver_Volume*) volume_ptr;
	struct Complex_DPG_Solver_Volume* c_dpg_s_vol = (struct Complex_DPG_Solver_Volume*) volume_ptr;

	const int order = s_vol->sol_coef->order;
	ptrdiff_t* extents = s_vol->sol_coef->extents;

	c_dpg_s_vol->sol_coef = constructor_empty_Multiarray_c('C',order,extents); // destructed
}

void destructor_derived_Complex_DPG_Solver_Volume (struct Volume* volume_ptr)
{
	struct Complex_DPG_Solver_Volume* c_dpg_s_vol = (struct Complex_DPG_Solver_Volume*) volume_ptr;

	destructor_Multiarray_c(c_dpg_s_vol->sol_coef);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
