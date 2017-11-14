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

#include "face_solver_dg.h"

#include "macros.h"

#include "simulation.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void constructor_derived_DG_Solver_Face (struct Face* face_ptr, const struct Simulation* sim)
{
	UNUSED(face_ptr);
	UNUSED(sim);
}

void destructor_derived_DG_Solver_Face (struct Face* face_ptr)
{
	UNUSED(face_ptr);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
