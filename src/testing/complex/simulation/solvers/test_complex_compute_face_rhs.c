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

#include "test_complex_compute_face_rhs.h"

#include "test_complex_boundary.h"
#include "test_complex_numerical_flux.h"

#include <assert.h>

#include "macros.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void destructor_Numerical_Flux_Input_c_data (struct Numerical_Flux_Input_c* num_flux_i)
{
	destructor_Boundary_Value_Input_c(&num_flux_i->bv_l);
	destructor_Boundary_Value_c(&num_flux_i->bv_r);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
