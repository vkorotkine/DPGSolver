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

#include "face_solver.h"

#include "volume.h"

#include "multiarray.h"

#include "const_cast.h"
#include "geometry.h"
#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Checks if one of the neighbouring volumes to the current face is curved.
 *  \return `true` if curved volume is found; `false` otherwise. */
bool check_for_curved_neigh
	(struct Face* face ///< \ref Face.
	);

// Templated functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "def_templates_multiarray.h"
#include "def_templates_boundary_d.h"
#include "def_templates_face_solver.h"
#include "face_solver_T.c"
#include "undef_templates_type.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_boundary.h"
#include "undef_templates_face_solver.h"

// Interface functions ********************************************************************************************** //

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

bool check_for_curved_neigh (struct Face* face)
{
	if (face->neigh_info[0].volume->curved || (face->neigh_info[1].volume && face->neigh_info[1].volume->curved))
		return true;
	return false;
}
