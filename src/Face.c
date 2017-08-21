// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "Face.h"
#include "Intrusive.h"
#include "Simulation.h"
#include "Mesh.h"

#include <stdlib.h>

#include "Macros.h"

// Static function declarations ************************************************************************************* //

/// \brief Destructor for an individual \ref Face.
static void destructor_Face
	(struct Face* face ///< Standard.
	);

// Interface functions ********************************************************************************************** //

struct Intrusive_List* constructor_Face_List (const struct Simulation*const sim, const struct Mesh*const mesh)
{
	struct Intrusive_List* Faces = constructor_empty_IL();

//	const struct const_Matrix_d*const nodes = mesh->mesh_data->nodes;
//	const struct const_Multiarray_Vector_ui*const node_nums = mesh->mesh_data->node_nums;

	const struct const_Multiarray_Vector_ui*const v_to_v  = mesh->mesh_conn->v_to_v,
	                                       *const v_to_lf = mesh->mesh_conn->v_to_lf;
UNUSED(v_to_v);
UNUSED(v_to_lf);
UNUSED(sim);

//	const unsigned int n_f =
/*	for (size_t f = 0; f < n_f; ++f) {
		push_back_IL(Faces,(struct Intrusive_Link*) constructor_Face(sim,&vol_mi));
	}
*/

	return Faces;
}

void destructor_Faces (struct Intrusive_List* faces)
{
	for (const struct Intrusive_Link* curr = faces->first; curr; ) {
		struct Intrusive_Link* next = curr->next;
		destructor_Face((struct Face*) curr);
		curr = next;
	}
	destructor_IL(faces);
}

static void destructor_Face (struct Face* face)
{
	UNUSED(face);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

