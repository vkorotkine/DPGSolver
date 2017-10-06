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

#include "solver_face.h"

#include <string.h>

#include "macros.h"
#include "definitions_intrusive.h"

#include "multiarray.h"

#include "const_cast.h"
#include "geometry.h"
#include "simulation.h"
#include "solution.h"
#include "volume.h"

// Static function declarations ************************************************************************************* //

/** \brief Constructor for an individual \ref Solver_Face.
 *  \return Standard. */
static struct Solver_Face* constructor_Solver_Face
	(struct Face* face,           ///< \ref Face.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for an individual \ref Solver_Face.
static void destructor_Solver_Face
	(struct Solver_Face* face ///< Standard.
	);

// Interface functions ********************************************************************************************** //

struct Intrusive_List* constructor_Solver_Faces (struct Simulation*const sim)
{
	struct Intrusive_List* faces        = sim->faces;
	struct Intrusive_List* solver_faces = constructor_empty_IL(IL_SOLVER_FACE,faces);

	for (struct Intrusive_Link* curr = faces->first; curr; curr = curr->next)
		push_back_IL(solver_faces,(struct Intrusive_Link*) constructor_Solver_Face((struct Face*)curr,sim));

	return solver_faces;
}

void destructor_Solver_Faces (struct Intrusive_List* solver_faces)
{
	for (const struct Intrusive_Link* curr = solver_faces->first; curr; ) {
		struct Intrusive_Link* next = curr->next;
		destructor_Solver_Face((struct Solver_Face*) curr);
		curr = next;
	}
	destructor_IL(solver_faces);
}

static void destructor_Solver_Face (struct Solver_Face* face)
{
	destructor_const_Multiarray_d(face->xyz_fc);
	destructor_const_Multiarray_d(face->n_fc);
	destructor_const_Multiarray_d(face->jacobian_det_fc);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Checks if one of the neighbouring volumes to the current face is curved.
 *  \return `true` if curved volume is found; `false` otherwise. */
bool check_for_curved_neigh
	(struct Face* face ///< \ref Face.
	);

static struct Solver_Face* constructor_Solver_Face (struct Face* face, const struct Simulation* sim)
{
	struct Solver_Face* solver_face = calloc(1,sizeof *solver_face); // returned
	memcpy(&solver_face->face,face,sizeof *face); // shallow copy of the base.

	set_derived_link(face,solver_face);
	set_derived_link(solver_face,NULL);

	const_cast_i(&solver_face->p_ref,sim->p_s_v[0]);
	const_cast_c(&solver_face->cub_type,(check_for_curved_neigh(face) ? 'c' : 's'));

	const_constructor_move_Multiarray_d(
		&solver_face->xyz_fc,constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){0,0}));        // destructed
	const_constructor_move_Multiarray_d(
		&solver_face->n_fc,constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){0,0}));          // destructed
	const_constructor_move_Multiarray_d(
		&solver_face->jacobian_det_fc,constructor_empty_Multiarray_d('C',1,(ptrdiff_t[]){0})); // destructed

	return solver_face;
}

// Level 1 ********************************************************************************************************** //

bool check_for_curved_neigh (struct Face* face)
{
	if (face->neigh_info[0].volume->curved || face->neigh_info[1].volume->curved)
		return true;
	return false;
}
