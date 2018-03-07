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

#include "test_complex_solution_diffusion.h"

#include "multiarray.h"

#include "boundary.h"
#include "compute_error.h"
#include "compute_error_diffusion.h"
#include "const_cast.h"
#include "file_processing.h"
#include "geometry.h"
#include "simulation.h"
#include "solution.h"

#include "test_complex_flux_diffusion.h"
#include "test_complex_numerical_flux_diffusion.h"
#include "default_steady/test_complex_solution_diffusion_default_steady.h"
#include "test_complex_test_case.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

#include "def_templates_type_dc.h"
#include "solution_diffusion_T.c"

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
