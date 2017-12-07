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

#include "test_complex_solution_euler.h"

#include "macros.h"
#include "definitions_core.h"
#include "definitions_test_case.h"

#include "multiarray.h"

#include "boundary.h"
#include "compute_error_euler.h"
#include "const_cast.h"
#include "simulation.h"
#include "test_case.h"

#include "complex_multiarray.h"
#include "test_complex_flux_euler.h"
#include "test_complex_geometry.h"
#include "test_complex_geometry_parametric.h"
#include "test_complex_numerical_flux_euler.h"
#include "test_complex_solution.h"
#include "periodic_vortex/test_complex_solution_periodic_vortex.h"
#include "supersonic_vortex/test_complex_solution_supersonic_vortex.h"
#include "test_complex_test_case.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

#include "def_templates_type_dc.h"
#include "solution_euler_T.c"

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
