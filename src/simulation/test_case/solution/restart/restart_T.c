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

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"

#include "def_templates_restart.h"

#include "def_templates_multiarray.h"
#include "def_templates_solution.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void set_sol_restart_T (const struct Simulation*const sim, struct Solution_Container_T sol_cont)
{
	EXIT_ADD_SUPPORT; UNUSED(sim); UNUSED(sol_cont);
}

const struct const_Multiarray_T* constructor_const_sol_restart_T
	(const struct const_Multiarray_R*const xyz, const struct Simulation*const sim)
{
	EXIT_ADD_SUPPORT; UNUSED(sim); UNUSED(xyz); return NULL;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
