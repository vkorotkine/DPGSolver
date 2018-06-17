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

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"

#include "def_templates_penalty_opg.h"
#include "def_templates_flux.h"
#include "def_templates_numerical_flux.h"
#include "def_templates_face_solver.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void constructor_rlhs_f_test_penalty_unsupported_T
	(const struct Flux_T*const flux, const struct Numerical_Flux_T*const num_flux, struct Solver_Face_T*const s_face,
	 struct Solver_Storage_Implicit*const ssi)
{
	UNUSED(flux); UNUSED(num_flux); UNUSED(s_face); UNUSED(ssi);
	EXIT_UNSUPPORTED;
}

void constructor_rhs_f_test_penalty_do_nothing_T
	(const struct Flux_T*const flux, const struct Numerical_Flux_T*const num_flux, struct Solver_Face_T*const s_face,
	 struct Solver_Storage_Implicit*const ssi)
{
	UNUSED(flux); UNUSED(num_flux); UNUSED(s_face); UNUSED(ssi);
	return;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

#include "undef_templates_penalty_opg.h"
#include "undef_templates_flux.h"
#include "undef_templates_numerical_flux.h"
#include "undef_templates_face_solver.h"
