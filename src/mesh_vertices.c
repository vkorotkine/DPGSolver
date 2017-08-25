// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "mesh_vertices.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Macros.h"
#include "constants_mesh.h"
#include "constants_bc.h"

#include "Element.h"
#include "Multiarray.h"
#include "Vector.h"

#include "mesh_geometry_cylinder_hollow_section.h"

// Static function declarations ************************************************************************************* //

/// \brief Container for locally computed \ref Mesh_Vertices members.
struct Mesh_Vertices_l {
	struct Vector_i* ve_curved;        // Local version of variable defined in \ref Mesh_Vertices.
	struct Vector_i* ve_boundary;      // Local version of variable defined in \ref Mesh_Vertices.
	struct Multiarray_Vector_i* ve_bc; // Local version of variable defined in \ref Mesh_Vertices.
};

/** \brief Constructor for the \ref Mesh_Vertices.
 *	\return Standard. */
static struct Mesh_Vertices* constructor_Mesh_Vertices
	(const struct Mesh_Vertices_l*const mesh_vert_l ///< \ref Mesh_Vertices_l.
	);

/// \brief Correct the coordinates of the vertex nodes in curved meshes.
static void correct_mesh_vertices
	(const struct Mesh_Input*const mesh_input,   ///< \ref Mesh_Input.
	 const struct Mesh_Vertices*const mesh_vert, ///< \ref Mesh_Vertices.
	 const struct const_Matrix_d*const nodes     ///< Defined in \ref Mesh_Data.
	);

// Interface functions ********************************************************************************************** //

struct Mesh_Vertices* mesh_process_vertices
	(const struct Mesh*const mesh, const struct const_Intrusive_List* elements,
	 const struct Mesh_Input*const mesh_input)
{
	// Set up Mesh_Vertices
	struct Mesh_Vertices_l mesh_vert_l;

	const ptrdiff_t ind_v = mesh->mesh_data->ind_v,
	                n_ve = mesh->mesh_data->nodes->extents[0];

	const struct const_Vector_i*const       elem_types  = mesh->mesh_data->elem_types,
	                           *const*const volume_nums = &mesh->mesh_data->node_nums->data[ind_v];

	const struct const_Multiarray_Vector_i*const v_to_v  = mesh->mesh_conn->v_to_v,
	                                      *const v_to_lf = mesh->mesh_conn->v_to_lf;


	mesh_vert_l.ve_curved   = constructor_empty_Vector_i(n_ve); // keep
	mesh_vert_l.ve_boundary = constructor_empty_Vector_i(n_ve); // keep
	mesh_vert_l.ve_bc       = constructor_empty_Multiarray_Vector_i (1,n_ve); // keep
	set_to_zero_Vector_i(mesh_vert_l.ve_curved);
	set_to_zero_Vector_i(mesh_vert_l.ve_boundary);

	const ptrdiff_t n_v = v_to_lf->extents[0];
	for (ptrdiff_t v = 0; v < n_v; ++v) {
		const struct const_Vector_i*const v_to_v_V  = v_to_v->data[v],
		                           *const v_to_lf_V = v_to_lf->data[v];

		const int lf_max = v_to_v_V->extents[0];
		for (int lf = 0; lf < lf_max; ++lf) {
			const int v_to_lf_i = v_to_lf_V->data[lf];
			if (v_to_lf_i < BC_STEP_SC)
				continue;

			struct const_Element* element = get_element_by_type(elements,elem_types->data[ind_v+v]);
			const struct const_Vector_i*const f_ve = element->f_ve->data[lf];
			const ptrdiff_t n_ve_f = f_ve->extents[0];

			for (int ve = 0; ve < n_ve_f; ++ve) {
				const int ind_ve = volume_nums[v]->data[f_ve->data[ve]];

				mesh_vert_l.ve_boundary->data[ind_ve] = 1;
				push_back_Vector_i(mesh_vert_l.ve_bc->data[ind_ve],v_to_lf_i,true,true);
				if (v_to_lf_i > 2*BC_STEP_SC)
					mesh_vert_l.ve_curved->data[ind_ve] = 1;
			}
		}
	}

	struct Mesh_Vertices* mesh_vert = constructor_Mesh_Vertices(&mesh_vert_l);

	// Correct vertex coordinates if necessary.
	correct_mesh_vertices(mesh_input,mesh_vert,mesh->mesh_data->nodes);

	return mesh_vert;
}

void destructor_Mesh_Vertices (struct Mesh_Vertices* mesh_vert)
{
	destructor_Vector_i((struct Vector_i*)mesh_vert->ve_curved);
	destructor_Vector_i((struct Vector_i*)mesh_vert->ve_boundary);
	destructor_Multiarray_Vector_i((struct Multiarray_Vector_i*)mesh_vert->ve_bc);

	free(mesh_vert);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Function pointer to mesh_snap_to_boundary functions.
 *	\param input_path The path to the input file.
 *	\param ve_curved  Defined in \ref Mesh_Vertices.
 *	\param ve_bc      Defined in \ref Mesh_Data.
 */
typedef void (*mesh_snap_to_boundary_fptr)
	(const char*const input_path,
	 const struct const_Vector_i*const ve_curved,
	 const struct Matrix_d*const nodes
	);

static struct Mesh_Vertices* constructor_Mesh_Vertices (const struct Mesh_Vertices_l*const mesh_vert_l)
{
	struct Mesh_Vertices* mesh_vert = malloc(sizeof *mesh_vert); // returned

	const_constructor_move_Vector_i(&mesh_vert->ve_curved,mesh_vert_l->ve_curved);
	const_constructor_move_Vector_i(&mesh_vert->ve_boundary,mesh_vert_l->ve_boundary);
	const_constructor_move_Multiarray_Vector_i(&mesh_vert->ve_bc,mesh_vert_l->ve_bc);

	return mesh_vert;
}

static mesh_snap_to_boundary_fptr set_fptr_mesh_snap (const char*const geom_name)
{
	if (strstr(geom_name,"n-Cylinder_HollowSection"))
		return mesh_snap_to_cylinder__hollow_section;

	EXIT_ADD_SUPPORT;
}

static void correct_mesh_vertices
	(const struct Mesh_Input*const mesh_input, const struct Mesh_Vertices*const mesh_vert,
	 const struct const_Matrix_d*const nodes)
{
	if (mesh_input->domain_type != DOM_CURVED)
		return;

	mesh_snap_to_boundary_fptr mesh_snap_to_boundary = set_fptr_mesh_snap(mesh_input->geom_name);

	mesh_snap_to_boundary(mesh_input->input_path,mesh_vert->ve_curved,(struct Matrix_d*)nodes);
}
