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

#include "computational_elements_c.h"

#include "face_solver_adaptive.h"
#include "face_solver_dg.h"
#include "face_solver_dpg.h"
#include "face_solver_opg.h"
#include "volume_solver_adaptive.h"
#include "volume_solver_dg.h"
#include "volume_solver_dpg.h"
#include "volume_solver_opg.h"

#include "element.h"
#include "element_geometry.h"
#include "element_plotting.h"
#include "element_solution.h"
#include "element_adaptation.h"
#include "element_solver.h"
#include "element_solver_dg.h"
#include "element_solver_dpg.h"
#include "element_solver_opg.h"

#include "computational_elements.h"
#include "intrusive.h"
#include "simulation.h"

// Templated functions ********************************************************************************************** //

#include "def_templates_type_dc.h"
#include "computational_elements_T.c"
#include "undef_templates_type.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
