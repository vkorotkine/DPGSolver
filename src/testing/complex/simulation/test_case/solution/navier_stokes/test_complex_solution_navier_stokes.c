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

#include "test_complex_solution_navier_stokes.h"

#include "macros.h"
#include "definitions_test_case.h"

#include "compute_error_euler.h"
#include "compute_error_navier_stokes.h"
#include "file_processing.h"
#include "simulation.h"
#include "solution_navier_stokes.h"

#include "complex_multiarray.h"

#include "test_complex_flux_euler.h"
#include "test_complex_flux_navier_stokes.h"
#include "test_complex_geometry.h"
#include "test_complex_geometry_parametric.h"
#include "test_complex_numerical_flux_euler.h"
#include "test_complex_numerical_flux_navier_stokes.h"
#include "test_complex_solution.h"
#include "test_complex_test_case.h"

#include "free_stream/test_complex_solution_free_stream.h"
#include "taylor_couette/test_complex_solution_taylor_couette.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

#include "def_templates_type_dc.h"
#include "solution_navier_stokes_T.c"

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
