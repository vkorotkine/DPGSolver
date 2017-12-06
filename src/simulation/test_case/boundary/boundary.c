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

#include "multiarray.h"

#include "element_solver.h"
#include "face.h"
#include "face_solver.h"
#include "volume.h"
#include "volume_solver.h"

#include "compute_face_rlhs.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "solve_dg.h"

// Templated functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "def_templates_multiarray.h"
#include "def_templates_boundary_d.h"
#include "def_templates_operators_d.h"
#include "def_templates_face_solver.h"
#include "boundary_T.c"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
