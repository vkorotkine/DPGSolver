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
#include "definitions_alloc.h"

#include "def_templates_face_solver_opg.h"
#include "def_templates_volume_solver_opg.h"

#include "def_templates_penalty_opg.h"
#include "def_templates_compute_rlhs.h"
#include "def_templates_numerical_flux.h"
#include "def_templates_face_solver.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void reset_penalty_indicators_opg_T (const struct Intrusive_List*const faces)
{
	EXIT_UNSUPPORTED; /// \todo Delete if unused.
	UNUSED(faces);
	/* for (struct Intrusive_Link* curr = faces->first; curr; curr = curr->next) { */
	/* 	const struct Face*const face = (struct Face*) curr; */
	/* 	if (!face->boundary) */
	/* 		continue; */

	/* 	struct OPG_Solver_Face_T*const opg_s_face = (struct OPG_Solver_Face_T*) curr; */
	/* 	if (!boundary_0) { */
	/* 		/\** Possibly required to choose this boundary condition for the most downwind face; hopefully */
	/* 		 *  not the case as this would be difficult for more complicated cases than linear advection. */
	/* 		 *  \todo Delete this comment after checking. *\/ */
	/* 		boundary_0 = true; */
	/* 		opg_s_face->bc_test_s = BC_TEST_S_SPEC_CONST; */
	/* 	} else { */
	/* 		opg_s_face->bc_test_s = BC_TEST_S_FREE_CONST; */
	/* 	} */
	/* 	struct OPG_Solver_Volume_T*const opg_s_volume = (struct OPG_Solver_Volume_T*) face->neigh_info[0].volume; */
	/* 	opg_s_volume->bc_test_s = BC_TEST_S_NEEDED; */
	/* } */
}

void constructor_rlhs_f_test_penalty_unsupported_T
	(const struct Flux_Ref_T*const flux_r, const struct Numerical_Flux_T*const num_flux,
	 struct Solver_Face_T*const s_face, struct Solver_Storage_Implicit*const ssi)
{
	UNUSED(flux_r); UNUSED(num_flux); UNUSED(s_face); UNUSED(ssi);
	EXIT_UNSUPPORTED;
}

void constructor_rlhs_f_test_penalty_do_nothing_T
	(const struct Flux_Ref_T*const flux_r, const struct Numerical_Flux_T*const num_flux,
	 struct Solver_Face_T*const s_face, struct Solver_Storage_Implicit*const ssi)
{
	UNUSED(flux_r); UNUSED(num_flux); UNUSED(s_face); UNUSED(ssi);
	return;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

#include "undef_templates_face_solver_opg.h"
#include "undef_templates_volume_solver_opg.h"

#include "undef_templates_penalty_opg.h"
#include "undef_templates_compute_rlhs.h"
#include "undef_templates_numerical_flux.h"
#include "undef_templates_face_solver.h"
