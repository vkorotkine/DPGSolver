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

#include "compute_error_advection.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "macros.h"
#include "definitions_error.h"

#include "volume.h"
#include "volume_solver.h"

#include "multiarray.h"
#include "vector.h"

#include "compute_error.h"
#include "const_cast.h"
#include "element.h"
#include "intrusive.h"
#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Return a statically allocated `char*` holding the specific header for all of the Euler variables.
 *  \return See brief. */
static const char* compute_header_spec_advection_all
	(const struct Simulation* sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

struct Error_CE* constructor_Error_CE_advection_all (const struct Simulation* sim)
{
	const int n_out = 1;

	struct Error_CE_Helper* e_ce_h = constructor_Error_CE_Helper(sim,n_out);
	e_ce_h->header_spec = compute_header_spec_advection_all(sim);
	const_cast_i(&e_ce_h->error_type,ERROR_STANDARD);

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		e_ce_h->s_vol[0] = (struct Solver_Volume*) curr;
		struct Error_CE_Data* e_ce_d = constructor_Error_CE_Data(e_ce_h,sim); // destructed

		increment_sol_L2(e_ce_h,e_ce_d);
		destructor_Error_CE_Data(e_ce_d);

		update_domain_order(e_ce_h);
	}

	struct Error_CE* error_ce = constructor_Error_CE(e_ce_h,sim); // returned
	destructor_Error_CE_Helper(e_ce_h);

	return error_ce;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static const char* compute_header_spec_advection_all (const struct Simulation* sim)
{
	UNUSED(sim);
	static char header_spec[STRLEN_MAX];

	sprintf(header_spec,"%-14s","$u$");

	return header_spec;
}
