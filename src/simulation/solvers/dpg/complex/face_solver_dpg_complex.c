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

#include "face_solver_dpg_complex.h"

#include "macros.h"

#include "complex_multiarray_minimal.h"
#include "multiarray.h"

#include "simulation.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void constructor_derived_Complex_DPG_Solver_Face (struct Face* face_ptr, const struct Simulation* sim)
{
	UNUSED(sim);
	const struct Solver_Face* s_face             = (struct Solver_Face*) face_ptr;
	struct Complex_DPG_Solver_Face* c_dpg_s_face = (struct Complex_DPG_Solver_Face*) face_ptr;

	const int order = s_face->nf_coef->order;
	ptrdiff_t* extents = s_face->nf_coef->extents;

	c_dpg_s_face->nf_coef = constructor_empty_Multiarray_c('C',order,extents); // destructed

	// Function pointer(s) are set as part of the complex solver functions such that test functions are not exposed
	// to the main code.
}

void destructor_derived_Complex_DPG_Solver_Face (struct Face* face_ptr)
{
	struct Complex_DPG_Solver_Face* c_dpg_s_face = (struct Complex_DPG_Solver_Face*) face_ptr;

	destructor_Multiarray_c(c_dpg_s_face->nf_coef);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
