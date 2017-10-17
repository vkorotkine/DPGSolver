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

#include "face_solver.h"

#include <string.h>

#include "macros.h"
#include "definitions_intrusive.h"

#include "multiarray.h"

#include "const_cast.h"
#include "geometry.h"
#include "simulation.h"
#include "volume.h"

// Static function declarations ************************************************************************************* //

/** \brief Checks if one of the neighbouring volumes to the current face is curved.
 *  \return `true` if curved volume is found; `false` otherwise. */
bool check_for_curved_neigh
	(struct Face* face ///< \ref Face.
	);

// Interface functions ********************************************************************************************** //

void constructor_derived_Solver_Face (struct Face* face_ptr, const struct Simulation* sim)
{
	struct Solver_Face* face = (struct Solver_Face*) face_ptr;

	const_cast_i(&face->p_ref,sim->p_s_v[0]);
	const_cast_c(&face->cub_type,(check_for_curved_neigh((struct Face*)face) ? 'c' : 's'));

	face->sol_coef  = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){0,0});   // destructed
	face->grad_coef = constructor_empty_Multiarray_d('C',3,(ptrdiff_t[]){0,0,0}); // destructed

	const_constructor_move_Multiarray_d(
		&face->xyz_fc,constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){0,0}));        // destructed
	const_constructor_move_Multiarray_d(
		&face->normals_fc,constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){0,0}));    // destructed
	const_constructor_move_Multiarray_d(
		&face->jacobian_det_fc,constructor_empty_Multiarray_d('C',1,(ptrdiff_t[]){0})); // destructed
}

void destructor_derived_Solver_Face (struct Face* face_ptr)
{
	struct Solver_Face* face = (struct Solver_Face*) face_ptr;

	destructor_Multiarray_d(face->sol_coef);
	destructor_Multiarray_d(face->grad_coef);

	destructor_const_Multiarray_d(face->xyz_fc);
	destructor_const_Multiarray_d(face->normals_fc);
	destructor_const_Multiarray_d(face->jacobian_det_fc);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

bool check_for_curved_neigh (struct Face* face)
{
	if (face->neigh_info[0].volume->curved || (face->neigh_info[1].volume && face->neigh_info[1].volume->curved))
		return true;
	return false;
}
