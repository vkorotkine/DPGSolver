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

#include "face_solver_dpg.h"

#include "simulation.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "face_solver_dpg_T.c"

void copy_members_r_to_c_DPG_Solver_Face
	(struct DPG_Solver_Face_c*const dpg_s_face, const struct DPG_Solver_Face*const dpg_s_face_r,
	 const struct Simulation*const sim)
{
	UNUSED(sim);
	UNUSED(dpg_s_face);
	UNUSED(dpg_s_face_r);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
