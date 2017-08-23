// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "Volume.h"
#include "Intrusive.h"
#include "Matrix.h"
#include "Simulation.h"
#include "Mesh.h"
#include "constants_elements.h"

#include <stdlib.h>
#include <stdbool.h>

#include "Macros.h"

#include "Multiarray.h"
#include "Vector.h"
#include "const_cast.h"

#include "constants_mesh.h"
#include "constants_bc.h"

// Static function declarations ************************************************************************************* //

/// \brief Container for \ref Volume related mesh information.
struct Volume_mesh_info {
	int elem_type;                        ///< The type of \ref Element associated with the volume.
	const struct const_Vector_i* ve_inds; ///< The indices of the vertices of the volume.

	const struct const_Vector_i* to_lf;   ///< The relevant row of \ref Mesh_Connectivity::v_to_lf.
};

/// \brief Constructor for an individual \ref Volume.
static struct Volume* constructor_Volume
	(const struct Simulation*const sim,          ///< The \ref Simulation.
	 const struct Volume_mesh_info*const vol_mi, ///< The \ref Volume_mesh_info.
	 const struct const_Matrix_d*const nodes,    ///< \ref Mesh_Data::nodes.
	 const struct Mesh_Vertices*const mesh_vert  ///< \ref Mesh_Vertices.
	);

/// \brief Destructor for an individual \ref Volume.
static void destructor_Volume
	(struct Volume* volume ///< Standard.
	);

// Interface functions ********************************************************************************************** //

struct Intrusive_List* constructor_Volume_List (struct Simulation*const sim, const struct Mesh*const mesh)
{
	struct Intrusive_List* volumes = constructor_empty_IL();

	const struct const_Matrix_d*const nodes = mesh->mesh_data->nodes;

	const struct const_Vector_i*const            elem_types = mesh->mesh_data->elem_types;
	const struct const_Multiarray_Vector_i*const node_nums = mesh->mesh_data->node_nums;

	const struct const_Multiarray_Vector_i*const v_to_lf = mesh->mesh_conn->v_to_lf;

	const ptrdiff_t n_v = v_to_lf->extents[0];
	for (ptrdiff_t v = 0; v < n_v; ++v) {
		const ptrdiff_t ind_v = v + mesh->mesh_data->ind_v;
		struct Volume_mesh_info vol_mi =
			{ .elem_type = elem_types->data[ind_v],
			  .ve_inds   = node_nums->data[ind_v],
			  .to_lf     = v_to_lf->data[v],
			};

		push_back_IL(volumes,(struct Intrusive_Link*) constructor_Volume(sim,&vol_mi,nodes,mesh->mesh_vert));
	}
	sim->n_v = n_v;

EXIT_UNSUPPORTED;
/// \bug The volume of index 6 is being marked as being on the boundary...

	return volumes;
}

void destructor_Volumes (struct Intrusive_List* volumes)
{
	for (const struct Intrusive_Link* curr = volumes->first; curr; ) {
		struct Intrusive_Link* next = curr->next;
		destructor_Volume((struct Volume*) curr);
		curr = next;
	}
	destructor_IL(volumes);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Check if the current volume is on a boundary.
 *
 *	This function works similarly to \ref check_if_curved_v.
 *
 *	\return True if on a boundary. */
static bool check_if_boundary_v
	(const struct const_Vector_i*const to_lf,           ///< Current volume component of \ref Mesh_Connectivity::v_to_lf.
	 const struct const_Multiarray_Vector_i*const f_ve, ///< Defined in \ref Element.
	 const struct const_Vector_i*const ve_inds,         ///< The vertex indices for the volume.
	 const struct const_Vector_i*const ve_boundary      ///< Defined in \ref Mesh_Vertices.
	);

/** \brief Check if the current volume is curved.
 *
 *	This function returns `true` if:
 *		1. The domain is mapped (and all elements are curved); or
 *		2. The volume contains a face which has at least two vertices which were marked as curved.
 *
 *	It is not enough to simply check if any of the adjacent faces lie on a curved boundary as this does not account for
 *	3D volumes which only have a curved edge.
 *
 *	\return True if curved. */
static bool check_if_curved_v
	(const int domain_type,                             ///< Defined in \ref Simulation.
	 const struct const_Multiarray_Vector_i*const f_ve, ///< Defined in \ref Element.
	 const struct const_Vector_i*const ve_inds,         ///< The vertex indices for the volume.
	 const struct const_Vector_i*const ve_curved        ///< Defined in \ref Mesh_Vertices.
	);

/** \brief Constructor for the xyz coordinates of the volume vertices.
 *	\return See brief. */
static struct Matrix_d* constructor_volume_vertices
	(const struct const_Vector_i*const ve_inds, ///< The vertex indices.
	 const struct const_Matrix_d*const nodes    ///< \ref Mesh_Data::nodes.
	);

static struct Volume* constructor_Volume
	(const struct Simulation*const sim, const struct Volume_mesh_info*const vol_mi,
	 const struct const_Matrix_d*const nodes, const struct Mesh_Vertices*const mesh_vert)
{
	struct Volume* volume = malloc(sizeof *volume); // returned


	const_constructor_move_Matrix_d(&volume->xyz_ve,constructor_volume_vertices(vol_mi->ve_inds,nodes));

	for (int i = 0; i < NFMAX;    ++i) {
	for (int j = 0; j < NSUBFMAX; ++j) {
		const_cast_Face(&volume->faces[i][j],NULL);
	}}

	const_cast_const_Element(&volume->element,get_element_by_type(sim->elements,vol_mi->elem_type));

	const_cast_bool(&volume->boundary,
	                check_if_boundary_v(vol_mi->to_lf,volume->element->f_ve,vol_mi->ve_inds,mesh_vert->ve_boundary));

	const_cast_bool(&volume->curved,
	                check_if_curved_v(sim->domain_type,volume->element->f_ve,vol_mi->ve_inds,mesh_vert->ve_curved));
printf("%d %d %d\n\n\n",volume->element->type,volume->boundary,volume->curved);

	return volume;
}

static void destructor_Volume (struct Volume* volume)
{
	UNUSED(volume);
}

// Level 1 ********************************************************************************************************** //

/** \brief Check if a sufficient number of vertices satisfy the condition.
 *
 *	Currently used to check:
 *		- curved;
 *		- boundary.
 *
 *	\return See brief. */
static bool check_for_ve_condition
	(const struct const_Multiarray_Vector_i*const f_ve, ///< Defined in \ref Element.
	 const struct const_Vector_i*const ve_inds,         ///< The vertex indices for the volume.
	 const struct const_Vector_i*const ve_condition     ///< The vertex condition to check for.
	);

static bool check_if_boundary_v
	(const struct const_Vector_i*const to_lf, const struct const_Multiarray_Vector_i*const f_ve,
	 const struct const_Vector_i*const ve_inds, const struct const_Vector_i*const ve_boundary)
{
	const ptrdiff_t i_max = to_lf->extents[0];
	for (ptrdiff_t i = 0; i < i_max; ++i) {
		if (to_lf->data[i] > BC_STEP_SC)
			return true;
	}

	return check_for_ve_condition(f_ve,ve_inds,ve_boundary);
}

static bool check_if_curved_v
	(const int domain_type, const struct const_Multiarray_Vector_i*const f_ve,
	 const struct const_Vector_i*const ve_inds, const struct const_Vector_i*const ve_curved)
{
	if (domain_type == DOM_MAPPED)
		return true;

	return check_for_ve_condition(f_ve,ve_inds,ve_curved);
}

static struct Matrix_d* constructor_volume_vertices
	(const struct const_Vector_i*const ve_inds, const struct const_Matrix_d*const nodes)
{
	struct Matrix_d* dest = constructor_empty_Matrix_d('R',ve_inds->extents[0],nodes->extents[1]);

	const ptrdiff_t i_max = dest->extents[0];
	for (ptrdiff_t i = 0; i < i_max; ++i)
		set_row_Matrix_d(i,dest,get_row_const_Matrix_d(ve_inds->data[i],nodes));

print_Matrix_d(dest);
	return dest;
}

// Level 2 ********************************************************************************************************** //

static bool check_for_ve_condition
	(const struct const_Multiarray_Vector_i*const f_ve, const struct const_Vector_i*const ve_inds,
	 const struct const_Vector_i*const ve_condition)
{
	const int n_lf = f_ve->extents[0];
	for (int lf = 0; lf < n_lf; ++lf) {
		const struct const_Vector_i*const f_ve_f = f_ve->data[lf];

		int count_ve_curved = 0;

		const int n_ve_f = f_ve_f->extents[0];
		for (int ve = 0; ve < n_ve_f; ++ve) {
			const ptrdiff_t ind_ve = ve_inds->data[f_ve_f->data[ve]];
			if (ve_condition->data[ind_ve]) {
				if (++count_ve_curved == 2)
					return true;
			}
		}
	}
	return false;
}
