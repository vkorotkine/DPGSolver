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

#include "test_complex_solution_advection.h"

#include "multiarray.h"

#include "boundary.h"
#include "compute_error.h"
#include "compute_error_advection.h"
#include "const_cast.h"
#include "file_processing.h"
#include "simulation.h"
#include "solution.h"

#include "test_complex_flux_advection.h"
#include "test_complex_numerical_flux_advection.h"
#include "test_complex_test_case.h"

#include "test_complex_solution_advection_default.h"
#include "free_stream_advection/test_complex_solution_free_stream_advection.h"
#include "peterson/test_complex_solution_peterson.h"
#include "vortex_advection/test_complex_solution_vortex_advection.h"
#include "hyperbolic_tan/test_complex_solution_hyperbolic_tan.h"
#include "gaussian_bump/test_complex_solution_gaussian_bump.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

#include "def_templates_type_dc.h"
#include "solution_advection_T.c"

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
