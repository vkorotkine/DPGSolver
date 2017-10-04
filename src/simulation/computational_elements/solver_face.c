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

#include "macros.h"

#include "geometry.h"
#include "simulation.h"
#include "solution.h"

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
	struct Intrusive_List* solver_faces = constructor_empty_IL(IL_SOLVER_FACE);

	for (struct Intrusive_Link* curr = faces->first; curr; curr = curr->next)
		push_back_IL(solver_faces,(struct Intrusive_Link*) constructor_Solver_Face((struct Face*)curr,sim));

	destructor_IL(faces);
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
	EXIT_ADD_SUPPORT;
	destructor_Face((struct Face*) face);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct Solver_Face* constructor_Solver_Face (struct Face* face, const struct Simulation* sim)
{
	struct Solver_Face* solver_face = calloc(1,sizeof *solver_face); // returned

	memcpy(&solver_face->face,face,sizeof *face); // shallow copy of the base.
UNUSED(sim);
EXIT_ADD_SUPPORT;

	return solver_face;
}
