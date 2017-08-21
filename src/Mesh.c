// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 *	\brief Reads a mesh file and returns relevant data as part of

 \todo Complete this.
 */

#include "Mesh.h"
#include "Intrusive.h"
#include "mesh_readers.h"
#include "mesh_connectivity.h"

#include <limits.h>
#include <string.h>

#include "Macros.h"
#include "Element.h"

#include "file_processing.h"

// Static function declarations ************************************************************************************* //

/** \brief Constructor for a \ref Mesh.
 *	\return Standard. */
struct Mesh* constructor_Mesh ();

// Interface functions ********************************************************************************************** //

struct Mesh* set_up_mesh
	(const char*const mesh_name_full, const unsigned int d, const struct const_Intrusive_List* elements)
{
	struct Mesh* mesh = constructor_Mesh();

	*(struct Mesh_Data**)&         mesh->mesh_data = mesh_reader(mesh_name_full,d);
	*(struct Mesh_Connectivity**)& mesh->mesh_conn = mesh_connect(mesh->mesh_data,elements);

	return mesh;
}

void destructor_Mesh (struct Mesh* mesh)
{
	destructor_Mesh_Data((struct Mesh_Data*)mesh->mesh_data);
	destructor_Mesh_Connectivity((struct Mesh_Connectivity*)mesh->mesh_conn);
	free(mesh);
}

size_t get_first_volume_index (const struct const_Vector_ui*const elem_per_dim, const unsigned int d)
{
	size_t ind = 0;
	for (unsigned int dim = 0; dim < d; dim++)
		ind += elem_per_dim->data[dim];
	return ind;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

struct Mesh* constructor_Mesh ()
{
	struct Mesh* mesh = malloc(1 * sizeof *mesh); // returned
	return mesh;
}
