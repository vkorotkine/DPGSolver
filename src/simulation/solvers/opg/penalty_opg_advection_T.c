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

#include "def_templates_penalty_opg_advection.h"
#include "def_templates_flux.h"
#include "def_templates_numerical_flux.h"
#include "def_templates_face_solver.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void constructor_rlhs_f_test_penalty_advection_upwind_T
	(const struct Flux_T*const flux, const struct Numerical_Flux_T*const num_flux, struct Solver_Face_T*const s_face,
	 struct Solver_Storage_Implicit*const ssi)
{
	UNUSED(flux);

	/** It is assumed that the initial solution for \ref Solver_Volume_T::test_s_coef is equal to zero on all
	 *  boundary faces which require the addition of the penalty term (outflow faces). It is also currently assumed
	 *  that the exact solution test functions are equal to zero (i.e. that the value of `g` described in
	 *  \ref constructor_rlhs_f_b_test_penalty_T is equal to zero).
	 *
	 *  While this is consistent with the test space used by Brunken et al. \cite Brunken2018, it has the
	 *  disadvantage of the exact test functions having jump discontinuities on faces which have a combination of
	 *  inflow and outflow boundary conditions. It would likely be desirable to set the value of the boundary
	 *  condition for the test functions based on the values of the solution test functions interpolated to the
	 *  inflow boundary nodes. \todo Investigate when applicable.
	 *
	 *  It is also possible that imposing the boundary conditions for the solution test functions at the face
	 *  cubature nodes does not set the physically correct number of boundary conditions when the number of cubature
	 *  nodes is not equal to the number of face basis functions. It may be required to set values for the numerical
	 *  solution in a face basis and subsequently interpolate to the face cubature nodes. Note that an identical
	 *  problem may be present for the treatment of the face terms for the DG scheme as well...
	 */

	assert(num_flux->nnf != NULL);
	; // do nothing (currently assuming that \f$ g = 0 \f$.

	if (num_flux->neigh_info[0].dnnf_ds != NULL) {

	}
	UNUSED(s_face); UNUSED(ssi);
	EXIT_ADD_SUPPORT;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

#include "undef_templates_penalty_opg_advection.h"
#include "undef_templates_flux.h"
#include "undef_templates_numerical_flux.h"
#include "undef_templates_face_solver.h"
