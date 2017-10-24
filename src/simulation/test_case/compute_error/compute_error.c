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

// Static function declarations ************************************************************************************* //

/// \brief Destructor for a \ref Error_CE container.
static void destructor_Error_CE
	(struct Error_CE* error_ce ///< Standard.
	);

// Interface functions ********************************************************************************************** //

void output_error (const struct Simulation* sim)
{
	struct Error_CE* error_ce = sim->test_case->constructor_Error_CE(sim);

	destructor_Error_CE(error_ce);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void destructor_Error_CE (struct Error_CE* error_ce)
{
	destructor_const_Vector_d(sol_L2);
}
