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

#include "face_solver_dg.h"

#include <string.h>

#include "macros.h"
#include "definitions_test_case.h"

#include "face.h"
#include "volume_solver.h"

#include "multiarray.h"

#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void constructor_derived_DG_Solver_Face (struct Face* face_ptr, const struct Simulation* sim)
{
	UNUSED(sim);
	struct DG_Solver_Face* face = (struct DG_Solver_Face*) face_ptr;

	for (int i = 0; i < 2; ++i) {
		struct Solver_Volume* volume = (struct Solver_Volume*) face_ptr->neigh_info[i].volume;
		if (volume) {
			const int order    = volume->sol_coef->order;
			ptrdiff_t* extents = volume->sol_coef->extents;
			face->rhs[i] = constructor_empty_Multiarray_d('C',order,extents); // destructed
		} else {
			face->rhs[i] = NULL;
		}
	}
}

void destructor_derived_DG_Solver_Face (struct Face* face_ptr)
{
	struct DG_Solver_Face* face = (struct DG_Solver_Face*) face_ptr;

	destructor_Multiarray_d(face->rhs[0]);
	if (face->rhs[1])
		destructor_Multiarray_d(face->rhs[1]);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
