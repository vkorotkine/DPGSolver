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

#include "computational_elements.h"

#include <assert.h>

#include "macros.h"
#include "definitions_elements.h"

#include "volume.h"
#include "face.h"

#include "intrusive.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

/// \brief Update pointers in \ref Volume to the derived list links.
static void update_volume_pointers
	(struct Intrusive_Link* link ///< The current link.
	);

/// \brief Update pointers in \ref Face to the derived list links.
static void update_face_pointers
	(struct Intrusive_Link* link ///< The current link.
	);

// Interface functions ********************************************************************************************** //

void update_computational_element_list_pointers (const struct Simulation* sim)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next)
		update_volume_pointers(curr);

	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next)
		update_face_pointers(curr);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void update_volume_pointers (struct Intrusive_Link* link)
{
	struct Volume* volume = (struct Volume*) link;

	for (int i = 0; i < NFMAX; ++i) {
	for (int j = 0; j < NSUBFMAX; ++j) {
		struct Face* face = (struct Face*) volume->faces[i][j];
		if (face) {
			const struct Face*const derived = (const struct Face*const)((struct Intrusive_Link*)face)->derived;
			assert(derived != NULL);
			const_cast_Face(&volume->faces[i][j],derived);
		}
	}}
}

static void update_face_pointers (struct Intrusive_Link* link)
{
	struct Face* face = (struct Face*) link;

	for (int i = 0; i < 2; ++i) {
		struct Volume* volume = face->neigh_info[i].volume;
		if (volume) {
			struct Volume* derived = (struct Volume*)((struct Intrusive_Link*)volume)->derived;
			assert(derived != NULL);
			face->neigh_info[i].volume = derived;
		}
	}
}
