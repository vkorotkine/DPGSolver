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

#include "numerical_flux.h"

#include <assert.h>

#include "definitions_numerical_flux.h"

#include "multiarray.h"
#include "vector.h"

#include "element_solver_dg.h"
#include "face.h"
#include "face_solver.h"
#include "volume.h"
#include "volume_solver.h"

#include "const_cast.h"
#include "geometry_normals.h"
#include "math_functions.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "solve_dg.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "numerical_flux_T.c"

const int* get_set_ind_num_flux (const int*const new_vals)
{
	static int ind_num_flux[2] = { -1, -1, };
	if (new_vals) {
		assert((new_vals[0] == NUM_FLUX_INVALID)    ||
		       (new_vals[0] == NUM_FLUX_UPWIND)     ||
		       (new_vals[0] == NUM_FLUX_ROE_PIKE));
		assert((new_vals[1] == NUM_FLUX_INVALID)    ||
		       (new_vals[1] == NUM_FLUX_BR2_STABLE) ||
		       (new_vals[1] == NUM_FLUX_CDG2));

		for (int i = 0; i < 2; ++i)
			ind_num_flux[i] = new_vals[i];
	}
	return ind_num_flux;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
