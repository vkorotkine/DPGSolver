// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "Volume.h"
#include "Intrusive.h"
#include "Simulation.h"
#include "Mesh.h"

#include <stdlib.h>
#include <stdbool.h>

#include "Macros.h"

#include "const_cast.h"

#include "constants_mesh.h"
#include "constants_bc.h"

// Static function declarations ************************************************************************************* //

struct Volume_mesh_info {
	unsigned int elem_type;
	const struct const_Vector_ui* ve_inds;

	const struct const_Vector_ui* to_v;
	const struct const_Vector_ui* to_lf;
};

/// \brief Constructor for an individual \ref Volume.
static struct Volume* constructor_Volume
	(const struct Simulation*const sim,          ///< The \ref Simulation.
	 const struct Volume_mesh_info*const vol_mi, ///< The \ref Volume_mesh_info.
	 const struct const_Matrix_d*const nodes     ///< \ref Mesh_Data::nodes.
	);

/// \brief Destructor for an individual \ref Volume.
static void destructor_Volume
	(struct Volume* volume ///< Standard.
	);

// Interface functions ********************************************************************************************** //

struct Intrusive_List* constructor_Volume_List (const struct Simulation*const sim, const struct Mesh*const mesh)
{
	struct Intrusive_List* Volumes = constructor_empty_IL();

	const struct const_Matrix_d*const nodes = mesh->mesh_data->nodes;

	const struct const_Vector_ui*const            elem_types = mesh->mesh_data->elem_types;
	const struct const_Multiarray_Vector_ui*const node_nums = mesh->mesh_data->node_nums;

	const struct const_Multiarray_Vector_ui*const v_to_v  = mesh->mesh_conn->v_to_v,
	                                       *const v_to_lf = mesh->mesh_conn->v_to_lf;

	const unsigned int n_v = v_to_v->extents[0];
	for (size_t v = 0; v < n_v; ++v) {
		struct Volume_mesh_info vol_mi =
			{ .elem_type = elem_types->data[v],
			  .ve_inds   = node_nums->data[v],
			  .to_lf     = v_to_lf->data[v],
			};

		push_back_IL(Volumes,(struct Intrusive_Link*) constructor_Volume(sim,&vol_mi,nodes));
	}

	return Volumes;
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

/** \brief Check if the current volume is on a boundary of the domain by searching for a boundary condition in `to_lf`.
 *	\return Positive if on a boundary. */
static bool check_if_boundary
	(const struct const_Vector_ui*const to_lf ///< Current volume component of \ref Mesh_Connectivity::v_to_lf.
	);

/** \brief Check if the current volume is curved.
 *	\return Positive if curved. */
static bool check_if_curved
	(const struct const_Vector_ui*const to_lf, ///< Current volume component of \ref Mesh_Connectivity::v_to_lf.
	 const unsigned int domain_type            ///< \ref Simulation::domain_type.
	);

/** \brief Constructor for the xyz coordinates of the volume vertices.
 *	\return See brief. */
static struct Matrix_d* constructor_volume_vertices
	(const struct Vector_ui*const ve_inds,   ///< The vertex indices.
	 const struct const_Matrix_d*const nodes ///< \ref Mesh_Data::nodes.
	);

static struct Volume* constructor_Volume
	(const struct Simulation*const sim, const struct Volume_mesh_info*const vol_mi,
	 const struct const_Matrix_d*const nodes)
{
	struct Volume* volume = malloc(sizeof *volume); // returned

	const_cast_bool(&volume->boundary,check_if_boundary(vol_mi->to_lf));
	const_cast_bool(&volume->curved,  check_if_curved(vol_mi->to_lf,sim->domain_type));

	const_constructor_move_Matrix_d(&volume->xyz_ve,constructor_volume_vertices(vol_mi->ve_inds,nodes));

	const_cast_const_Element(&volume->element,get_element_by_type(sim->elements,vol_mi->elem_type));

	return volume;
}

static void destructor_Volume (struct Volume* volume)
{
	UNUSED(volume);
}

// Level 1 ********************************************************************************************************** //

static bool check_if_boundary (const struct const_Vector_ui*const to_lf)
{
	const size_t i_max = to_lf->extents[0];
	for (size_t i = 0; i < i_max; ++i) {
		if (to_lf->data[i] > BC_STEP_SC)
			return true;
	}
	return false;
}

static bool check_if_curved (const struct const_Vector_ui*const to_lf, const unsigned int domain_type)
{
	if (domain_type == DOM_MAPPED)
		return true;

	const size_t i_max = to_lf->extents[0];
	for (size_t i = 0; i < i_max; ++i) {
		if (to_lf->data[i] > 2*BC_STEP_SC)
			return true;
	}
	return false;
}

static struct Matrix_d* constructor_volume_vertices
	(const struct Vector_ui*const ve_inds, const struct const_Matrix_d*const nodes)
{
	struct Matrix_d* dest = constructor_empty_Matrix_d('R',ve_inds->extents[0],nodes->extents[1]);

	const size_t i_max = dest->extents[0];
	for (size_t i = 0; i < i_max; ++i) {
		set_row
	}
}
