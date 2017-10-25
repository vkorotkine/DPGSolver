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

#include "compute_error.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "definitions_intrusive.h"

#include "multiarray.h"
#include "vector.h"

#include "computational_elements.h"
#include "element_error.h"
#include "volume_solver.h"

#include "intrusive.h"
#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Destructor for a \ref Error_CE container.
static void destructor_Error_CE
	(struct Error_CE* error_ce ///< Standard.
	);

// Interface functions ********************************************************************************************** //

void output_error (const struct Simulation* sim)
{
	assert(sim->volumes->name == IL_SOLVER_VOLUME);
	assert(sim->faces->name   == IL_SOLVER_FACE);

	constructor_derived_Elements((struct Simulation*)sim,IL_SOLUTION_ELEMENT);
	constructor_derived_Elements((struct Simulation*)sim,IL_ELEMENT_ERROR);

	struct Error_CE* error_ce = sim->test_case->constructor_Error_CE(sim);

	destructor_derived_Elements((struct Simulation*)sim,IL_SOLUTION_ELEMENT);
	destructor_derived_Elements((struct Simulation*)sim,IL_ELEMENT);

	destructor_Error_CE(error_ce);
}

double compute_volume (const struct Solver_Volume* s_vol)
{
	struct Volume* b_vol = (struct Volume*)s_vol;
	struct Error_Element* e = (struct Error_Element*) b_vol->element;

	const int curved = b_vol->curved,
	          p      = s_vol->p_ref;
	const struct const_Vector_d* w_vc = get_const_Multiarray_Vector_d(e->w_vc[curved],(ptrdiff_t[]){0,0,p,p});
	const struct const_Vector_d jacobian_det_vc = interpret_const_Multiarray_as_Vector_d(s_vol->jacobian_det_vc);
	assert(w_vc->ext_0 == jacobian_det_vc.ext_0);

	double volume = 0.0;
	const ptrdiff_t ext_0 = w_vc->ext_0;
	for (int i = 0; i < ext_0; ++i)
		volume += w_vc->data[i]*jacobian_det_vc.data[i];

	return volume;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void destructor_Error_CE (struct Error_CE* error_ce)
{
	destructor_const_Vector_d(error_ce->sol_L2);
}
