// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "volume.h"
#include "constants_elements.h"
#include "intrusive.h"

#include <stdlib.h>
#include <stdbool.h>

#include "macros.h"
#include "constants_mesh.h"
#include "constants_bc.h"
#include "constants_intrusive.h"

#include "multiarray.h"
#include "matrix.h"
#include "vector.h"

#include "simulation.h"
#include "mesh.h"
#include "mesh_readers.h"
#include "mesh_connectivity.h"
#include "mesh_vertices.h"
#include "const_cast.h"
#include "face.h"
#include "element.h"
#include "geometry.h"

// Static function declarations ************************************************************************************* //

/// \brief Container for \ref Volume related mesh information.
struct Volume_mesh_info {
	int elem_type;                        ///< The type of \ref Element associated with the volume.
	const struct const_Vector_i* ve_inds; ///< The indices of the vertices of the volume.

	const struct const_Vector_i* to_lf;   ///< The relevant row of \ref Mesh_Connectivity::v_to_lf.
};

/// \brief Constructor for an individual \ref Volume.
static struct Volume* constructor_Volume
	(const struct Simulation*const sim,          ///< \ref Simulation.
	 const struct Mesh*const mesh,               ///< \ref Mesh.
	 const struct Volume_mesh_info*const vol_mi, ///< \ref Volume_mesh_info.
	 const int index                             ///< The volume index.
	);

/** \brief Check for a match of boundaries relating to the two input vertex bc vectors.
 *
 *	The two input vectors of vertex boundary conditions are searched for a matching boundary condition possibly subject
 *	to the constraint of looking only at curved boundary conditions.
 *
 *	\return `true` if a matching boundary condition is found for the two vertices; `false` otherwise. */
static bool find_bc_match
	(const struct const_Vector_i*const ve_bc_0, ///< The vector of boundary conditions associated with the first vertex.
	 const struct const_Vector_i*const ve_bc_1, ///< The vector of boundary conditions associated with the second vertex.
	 const bool curved_only                     /**< Flag indicating whether only curved boundary conditions should be
	                                             *   considered. */
	);

// Interface functions ********************************************************************************************** //

struct Intrusive_List* constructor_Volumes (struct Simulation*const sim, const struct Mesh*const mesh)
{
	struct Intrusive_List* volumes = constructor_empty_IL(IL_VOLUME);

	const struct const_Vector_i*const            elem_types = mesh->mesh_data->elem_types;
	const struct const_Multiarray_Vector_i*const node_nums  = mesh->mesh_data->node_nums;

	const struct const_Multiarray_Vector_i*const v_to_lf = mesh->mesh_conn->v_to_lf;

	const ptrdiff_t n_v = compute_size(v_to_lf->order,v_to_lf->extents);
	for (ptrdiff_t v = 0; v < n_v; ++v) {
		const ptrdiff_t ind_v = v + mesh->mesh_data->ind_v;
		struct Volume_mesh_info vol_mi =
			{ .elem_type = elem_types->data[ind_v],
			  .ve_inds   = node_nums->data[ind_v],
			  .to_lf     = v_to_lf->data[v],
			};

		push_back_IL(volumes,(struct Intrusive_Link*) constructor_Volume(sim,mesh,&vol_mi,v));
	}
	sim->n_v = n_v;

	set_up_geometry(sim,volumes);

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

void destructor_Volume (struct Volume* volume)
{
	(volume->xyz_ve ? destructor_Matrix_d((struct Matrix_d*)volume->xyz_ve) : EXIT_DESTRUCTOR);
	(volume->geom_coef ? destructor_Matrix_d((struct Matrix_d*)volume->geom_coef) : EXIT_DESTRUCTOR);
	(volume->sol_coef ? destructor_Multiarray_d((struct Multiarray_d*)volume->sol_coef) : EXIT_DESTRUCTOR);
}

bool check_ve_condition
	(const struct const_Multiarray_Vector_i*const f_ve, const struct const_Vector_i*const ve_inds,
	 const struct const_Vector_i*const ve_condition, const struct const_Multiarray_Vector_i*const ve_bc,
	 const bool curved_only)
{
	const int n_lf = compute_size(f_ve->order,f_ve->extents);
	for (int lf = 0; lf < n_lf; ++lf) {
		const struct const_Vector_i*const f_ve_f = f_ve->data[lf];

		int count_ve_condition = 0;

		const struct const_Vector_i* ve_bc_V = NULL;

		const int n_ve_f = f_ve_f->ext_0;
		for (int ve = 0; ve < n_ve_f; ++ve) {
			const int ind_ve = ve_inds->data[f_ve_f->data[ve]];
			if (ve_condition->data[ind_ve]) {
				if (count_ve_condition == 0)
					ve_bc_V = ve_bc->data[ind_ve];

				// Ensure that the the vertices are on the same curved surface
				if (find_bc_match(ve_bc_V,ve_bc->data[ind_ve],curved_only)) {
					if (++count_ve_condition == 2)
						return true;
				}
			}
		}
	}
	return false;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Check if the current volume is on a boundary.
 *
 *	This function works similarly to \ref check_if_curved_v.
 *
 *	\return `true` if on a boundary. */
static bool check_if_boundary_v
	(const struct const_Vector_i*const to_lf,           ///< Current volume component of \ref Mesh_Connectivity::v_to_lf.
	 const struct const_Multiarray_Vector_i*const f_ve, ///< Defined in \ref Element.
	 const struct const_Vector_i*const ve_inds,         ///< The vertex indices for the volume.
	 const struct Mesh_Vertices*const mesh_vert         ///< \ref Mesh_Vertices.
	);

/** \brief Check if the current volume is curved.
 *
 *	This function returns `true` if:
 *		1. The domain is mapped (and all elements are curved); or
 *		2. The volume contains an edge which has at least two vertices which were marked as curved.
 *
 *	It is not enough to simply check if any of the adjacent faces lie on a curved boundary as this does not account for
 *	3D volumes which only have a curved edge.
 *
 *	\return `true` if curved. */
static bool check_if_curved_v
	(const int domain_type,                             ///< Defined in \ref Simulation.
	 const struct const_Multiarray_Vector_i*const f_ve, ///< Defined in \ref Element.
	 const struct const_Vector_i*const ve_inds,         ///< The vertex indices for the volume.
	 const struct Mesh_Vertices*const mesh_vert         ///< \ref Mesh_Vertices.
	);

/** \brief Constructor for the xyz coordinates of the volume vertices.
 *	\return See brief. */
static struct Matrix_d* constructor_volume_vertices
	(const struct const_Vector_i*const ve_inds, ///< The vertex indices.
	 const struct const_Matrix_d*const nodes    ///< \ref Mesh_Data::nodes.
	);

static struct Volume* constructor_Volume
	(const struct Simulation*const sim, const struct Mesh*const mesh,
	 const struct Volume_mesh_info*const vol_mi, const int index)
{
	const struct const_Matrix_d*const nodes = mesh->mesh_data->nodes;
	const struct Mesh_Vertices*const mesh_vert = mesh->mesh_vert;

	struct Volume* volume = calloc(1,sizeof *volume); // returned
	const_cast_i(&volume->index,index);

	const_constructor_move_Matrix_d(&volume->xyz_ve,constructor_volume_vertices(vol_mi->ve_inds,nodes));

	for (int i = 0; i < NFMAX;    ++i) {
	for (int j = 0; j < NSUBFMAX; ++j) {
		const_cast_Face(&volume->faces[i][j],NULL);
	}}

	const_cast_const_Element(&volume->element,get_element_by_type(sim->elements,vol_mi->elem_type));

	const_cast_bool(&volume->boundary,
	                check_if_boundary_v(vol_mi->to_lf,volume->element->f_ve,vol_mi->ve_inds,mesh_vert));
	const_cast_bool(&volume->curved,
	                check_if_curved_v(sim->domain_type,volume->element->f_ve,vol_mi->ve_inds,mesh_vert));

	const_constructor_move_Matrix_d(&volume->geom_coef,constructor_default_Matrix_d());
	const_constructor_move_Multiarray_d(&volume->sol_coef,constructor_default_Multiarray_d());

	return volume;
}

static bool find_bc_match
	(const struct const_Vector_i*const ve_bc_0, const struct const_Vector_i*const ve_bc_1, const bool curved_only)
{
	const ptrdiff_t i_max = ve_bc_0->ext_0;
	for (ptrdiff_t i = 0; i < i_max; ++i) {
		const int bc_0 = ve_bc_0->data[i];

		if (curved_only && bc_0 < 2*BC_STEP_SC)
			continue;

		if (find_val_Vector_i(ve_bc_1,bc_0,false))
			return true;
	}
	return false;
}

// Level 1 ********************************************************************************************************** //

static bool check_if_boundary_v
	(const struct const_Vector_i*const to_lf, const struct const_Multiarray_Vector_i*const f_ve,
	 const struct const_Vector_i*const ve_inds, const struct Mesh_Vertices*const mesh_vert)
{
	// If the volume has a face on a domain boundary
	const ptrdiff_t i_max = to_lf->ext_0;
	for (ptrdiff_t i = 0; i < i_max; ++i) {
		if (to_lf->data[i] > BC_STEP_SC)
			return true;
	}

	// If the volume has 2 vertices on a domain boundary (i.e. a boundary edge)
	return check_ve_condition(f_ve,ve_inds,mesh_vert->ve_boundary,mesh_vert->ve_bc,false);
}

static bool check_if_curved_v
	(const int domain_type, const struct const_Multiarray_Vector_i*const f_ve,
	 const struct const_Vector_i*const ve_inds, const struct Mesh_Vertices*const mesh_vert)
{
	if (domain_type == DOM_PARAMETRIC)
		return true;

	// If the volume has 2 vertices on a curved domain boundary (i.e. a curved boundary edge)
	return check_ve_condition(f_ve,ve_inds,mesh_vert->ve_curved,mesh_vert->ve_bc,true);
}

static struct Matrix_d* constructor_volume_vertices
	(const struct const_Vector_i*const ve_inds, const struct const_Matrix_d*const nodes)
{
	struct Matrix_d* dest = constructor_empty_Matrix_d('R',ve_inds->ext_0,nodes->ext_1);

	const ptrdiff_t i_max = dest->ext_0;
	for (ptrdiff_t i = 0; i < i_max; ++i)
		set_row_Matrix_d(i,dest,get_row_const_Matrix_d(ve_inds->data[i],nodes));

	return dest;
}
