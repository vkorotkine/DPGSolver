// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 *	\brief Reads a mesh file and returns relevant data as part of
 *
 *	General comments.
 */

#include "set_up_mesh.h"

#include <limits.h>
#include <string.h>

#include "Macros.h"
#include "Intrusive.h"
#include "Element.h"

#include "file_processing.h"
#include "mesh_readers.h"
#include "mesh_connectivity.h"


struct Mesh* constructor_Mesh ()
{
	struct Mesh* mesh = malloc(1 * sizeof *mesh); // returned
	return mesh;
}

void destructor_Mesh (struct Mesh* mesh)
{
	destructor_Mesh_Data((struct Mesh_Data*)mesh->mesh_data);
//	destructor_Mesh_Connectivity((struct Mesh_Connectivity*)mesh->mesh_connectivity);
	free(mesh);
}

struct Mesh* set_up_mesh (const char*const mesh_name_full, const unsigned int d)
{
	struct Mesh* mesh = constructor_Mesh();

// Set up Elements?
	struct Intrusive_List* Elements = constructor_Element_List(d);
UNUSED(Elements);

	*(struct Mesh_Data**)&         mesh->mesh_data = mesh_reader(mesh_name_full,d);
	*(struct Mesh_Connectivity**)& mesh->mesh_conn = mesh_connect(mesh->mesh_data);

//destructor_Mesh(mesh);

	return mesh;
}

