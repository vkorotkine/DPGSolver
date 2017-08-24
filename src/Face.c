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

#include "Volume.h"

#include "const_cast.h"

#include "constants_mesh.h"
#include "constants_bc.h"

// Static function declarations ************************************************************************************* //

/// \brief Container for \ref Volume related mesh information.
struct Face_mesh_info {
	struct const_Element* element;        ///< The pointer to the \ref Element corresponding to the face.

	const struct const_Vector_i* ve_inds; ///< The indices of the vertices of the face.
	int to_lf;                            ///< The boundary condition of the face.

	struct Neigh_info_mi {
		struct Volume* volume; /**< The pointers to the two adjacent \ref Volume elements. The second pointer is NULL
	                            *   for a boundary face. */
	} neigh_info[2];
};

/// \brief Constructor for an individual \ref Face.
static struct Face* constructor_Face
	(const struct Simulation*const sim,        ///< \ref Simulation.
	 const struct Mesh*const mesh,             ///< \ref Mesh.
	 const struct Face_mesh_info*const face_mi ///< \ref Face_mesh_info.
	);

/// \brief Destructor for an individual \ref Face.
static void destructor_Face
	(struct Face* face ///< Standard.
	);

// Interface functions ********************************************************************************************** //

struct Intrusive_List* constructor_Face_List (struct Simulation*const sim, const struct Mesh*const mesh)
{
	struct Intrusive_List* faces = constructor_empty_IL();

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
	ptrdiff_t n_f = 0;

	const ptrdiff_t v_max = sim->n_v;
	for (ptrdiff_t v = 0; v < v_max; ++v) {
		struct Volume*const volume_l = volume_array[v];

		const struct const_Vector_i*const v_to_v_V = v_to_v->data[v];

		const int lf_max = v_to_v_V->extents[0];
		for (int lf = 0; lf < lf_max; ++lf) {
			if (volume_l->faces[lf][0] != NULL) // Already found this face.
				continue;

			const int v_neigh = v_to_v_V->data[lf];
// set ve_inds from volume->element->f_ve of node_nums->data[ind_v]

			struct Face_mesh_info face_mi =
				{ .element              = get_element_by_face(volume_l->element,lf),
				  .ve_inds              = ve_inds,
				  .to_lf                = v_to_lf->data[v]->data[lf],
				  .neigh_info[0].volume = volume_l,
				  .neigh_info[1].volume = (v_neigh == -1) ? NULL : volume_array[v_neigh],
				};

			push_back_IL(faces,(struct Intrusive_Link*) constructor_Face(sim,mesh,&face_mi));

			struct Face* face = (struct Face*) faces->last;
			const_cast_Face(&volume_l->faces[lf][0],face);
			if (v_neigh != -1)
				const_cast_Face(&face_mi.neigh_info[1].volume->faces[face_mi.to_lf][0],face);

			++n_f;
		}
	}
	free(volume_array);

	sim->n_f = n_f;

	return faces;
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

void const_cast_Face (const struct Face*const* dest, const struct Face*const src)
{
	*(struct Face**) dest = (struct Face*) src;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Check if the current face is on a boundary.
 *
 *	This function works similarly to \ref check_if_boundary_v.
 *
 *	\return True if on a boundary. */
static bool check_if_boundary_f
	(const int to_lf,                                   ///< Current face component of \ref Mesh_Connectivity::v_to_lf.
	 const struct const_Multiarray_Vector_i*const f_ve, ///< Defined in \ref Element.
	 const struct const_Vector_i*const ve_inds,         ///< The vertex indices for the volume.
	 const struct Mesh_Vertices*const mesh_vert         ///< \ref Mesh_Vertices.
	);

/** \brief Check if the current face is curved.
 *
 *	This function works similarly to \ref check_if_curved_v.
 *
 *	\return True if curved. */
static bool check_if_curved_f
	(const int to_lf,                                   ///< Current face component of \ref Mesh_Connectivity::v_to_lf.
	 const int domain_type,                             ///< \ref Simulation::domain_type.
	 const struct const_Multiarray_Vector_i*const f_ve, ///< Defined in \ref Element.
	 const struct const_Vector_i*const ve_inds,         ///< The vertex indices for the volume.
	 const struct Mesh_Vertices*const mesh_vert         ///< \ref Mesh_Vertices.
	);

static struct Face* constructor_Face
	(const struct Simulation*const sim, const struct Mesh*const mesh, const struct Face_mesh_info*const face_mi)
{
	const struct Mesh_Vertices*const mesh_vert = mesh->mesh_vert;

	struct Face* face = malloc(sizeof *face); // returned

	for (int i = 0; i < 2; ++i) {
		if (face_mi->neigh_info[i].volume) {
			face->neigh_info[i].ind_lf   = face_mi->to_lf;
			face->neigh_info[i].ind_href = 0;
			face->neigh_info[i].ind_sref = 0;
			face->neigh_info[i].ind_ord  = -1; // Set to invalid value, to be subsequently determined.
			face->neigh_info[i].volume   = face_mi->neigh_info[i].volume;
		} else {
			face->neigh_info[i].ind_lf   = -1;
			face->neigh_info[i].ind_href = -1;
			face->neigh_info[i].ind_sref = -1;
			face->neigh_info[i].ind_ord  = -1;
			face->neigh_info[i].volume   = NULL;
		}
	}

	const_cast_const_Element(&face->element,face_mi->element);

	const_cast_bool(&face->boundary,
	                check_if_boundary_f(face_mi->to_lf,face->element->f_ve,face_mi->ve_inds,mesh_vert));
	const_cast_bool(&face->curved,
	                check_if_curved_f(face_mi->to_lf,sim->domain_type,face->element->f_ve,face_mi->ve_inds,mesh_vert));

	return face;
}

static void destructor_Face (struct Face* face)
{
	UNUSED(face);
}

// Level 1 ********************************************************************************************************** //

static bool check_if_boundary_f
	(const int to_lf, const struct const_Multiarray_Vector_i*const f_ve, const struct const_Vector_i*const ve_inds,
	 const struct Mesh_Vertices*const mesh_vert)
{
	if (to_lf > BC_STEP_SC)
		return true;

	return check_ve_condition(f_ve,ve_inds,mesh_vert->ve_boundary,mesh_vert->ve_bc,false);
}

static bool check_if_curved_f
	(const int to_lf, const int domain_type, const struct const_Multiarray_Vector_i*const f_ve,
	 const struct const_Vector_i*const ve_inds, const struct Mesh_Vertices*const mesh_vert)
{
	if (domain_type == DOM_MAPPED)
		return true;

// likely redundant {
	if (to_lf > 2*BC_STEP_SC)
		return true;
// }

	return check_ve_condition(f_ve,ve_inds,mesh_vert->ve_boundary,mesh_vert->ve_bc,true);
}
