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

/// \brief Container for \ref Volume related mesh information.
struct Face_mesh_info {
	int elem_type; ///< The type of \ref Element associated with the volume.

	int bc; ///< The boundary condition of the face.
};

/// \brief Constructor for an individual \ref Face.
static struct Face* constructor_Face
	(const struct Simulation*const sim,        ///< The \ref Simulation.
	 const struct Face_mesh_info*const face_mi ///< The \ref Face_mesh_info.
	);

/// \brief Destructor for an individual \ref Face.
static void destructor_Face
	(struct Face* face ///< Standard.
	);

// Interface functions ********************************************************************************************** //

struct Intrusive_List* constructor_Face_List (const struct Simulation*const sim, const struct Mesh*const mesh)
{
	struct Intrusive_List* Faces = constructor_empty_IL();

//	const struct const_Matrix_d*const nodes = mesh->mesh_data->nodes;
//	const struct const_Multiarray_Vector_i*const node_nums = mesh->mesh_data->node_nums;

	const struct const_Multiarray_Vector_i*const v_to_v  = mesh->mesh_conn->v_to_v,
	                                      *const v_to_lf = mesh->mesh_conn->v_to_lf;

	struct Volume** volume_array = malloc(sim->n_v * sizeof *volume_array); // free

	ptrdiff_t ind_v = 0;
	for (const struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		volume_array[ind_v] = (struct Volume*) curr;
		++ind_v;
	}

print_const_Multiarray_Vector_i(v_to_v);
print_const_Multiarray_Vector_i(v_to_lf);
	ptrdiff_t f = 0;

	const ptrdiff_t v_max = sim->n_v;
	for (ptrdiff_t v = 0; v < v_max; ++v) {
		struct Volume*const volume_l = volume_array[v];

		const struct const_Vector_i*const v_to_v_V = v_to_v->data[v];

		const int f_max = v_to_v_V->extents[0];
		for (int f = 0; f < f_max; ++f) {
			struct Face_mesh_info vol_mi =
				{ .elem_type = get_element_by_face(volume_l->element,f);
				  .bc        = v_to_lf->data[v]->data[f],
				};

			const int v_neigh = v_to_v_V->data[f];

			if (v_neigh == -1) { // Boundary face
				push_back_IL(Faces,(struct Intrusive_Link*) constructor_Face(sim,&face_mi));
				volume_l->faces[f][0] =
				volume_array[v]
			}

//			if (volume_array[v_neigh]->faces[
		UNUSED(v_to_v_V);
		UNUSED(f);
	}
	free(volume_array);

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

static struct Volume* constructor_Face (const struct Simulation*const sim, const struct Face_mesh_info*const face_mi)
{
	struct Face* face = malloc(sizeof *face); // returned

	const_cast_bool(&volume->boundary,check_if_boundary(face_mi->to_lf));
	const_cast_bool(&volume->curved,  check_if_curved(face_mi->to_lf,sim->domain_type));

	const_cast_const_Element(&face->element,get_element_by_type(sim->elements,face_mi->elem_type));

	return face;
}
