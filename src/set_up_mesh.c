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

#include "file_processing.h"
#include "mesh_readers.h"


struct Mesh* constructor_Mesh ()
{
	struct Mesh* mesh = malloc(1 * sizeof *mesh); // returned
	return mesh;
}

void destructor_Mesh (struct Mesh* mesh_data)
{
	free(mesh_data);
}

struct Mesh* set_up_mesh (const char*const mesh_name_full, const unsigned int d)
{
	struct Mesh* mesh = constructor_Mesh ();
	mesh_reader(mesh_name_full,d);
DO_NOTHING_P(mesh_name_full);
	return mesh;
}

