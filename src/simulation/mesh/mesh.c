/* {{{
This file is part of DPGSolver.

DPGSolver is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or any later version.

DPGSolver is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along with DPGSolver.  If not, see
<http://www.gnu.org/licenses/>.
}}} */
/** \file
 */

#include "mesh.h"

#include <limits.h>
#include <string.h>

#include "macros.h"

#include "vector.h"

#include "intrusive.h"
#include "mesh_readers.h"
#include "mesh_connectivity.h"
#include "mesh_vertices.h"
#include "element.h"
#include "file_processing.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

struct Mesh* constructor_Mesh (const struct Mesh_Input* mesh_input, const struct const_Intrusive_List* elements)
{
	bool elements_provided = true;
	if (elements == NULL) {
		elements_provided = false;
		elements = constructor_Elements(mesh_input->d); // destructed
	}

	struct Mesh* mesh = calloc(1,sizeof *mesh); // returned

	*(struct Mesh_Data**)&         mesh->mesh_data = constructor_Mesh_Data(mesh_input->mesh_name_full,mesh_input->d);
	*(struct Mesh_Connectivity**)& mesh->mesh_conn = constructor_Mesh_Connectivity(mesh->mesh_data,elements);
	*(struct Mesh_Vertices**)&     mesh->mesh_vert = constructor_Mesh_Vertices(mesh,elements,mesh_input);

	if (!elements_provided)
		destructor_const_Elements(elements);
	else
		set_elements_present(elements,mesh->mesh_data->elem_types);

	return mesh;
}

void destructor_Mesh (struct Mesh* mesh)
{
	destructor_Mesh_Data((struct Mesh_Data*)mesh->mesh_data);
	destructor_Mesh_Connectivity((struct Mesh_Connectivity*)mesh->mesh_conn);
	destructor_Mesh_Vertices((struct Mesh_Vertices*)mesh->mesh_vert);
	free(mesh);
}

ptrdiff_t get_first_volume_index (const struct const_Vector_i*const elem_per_dim, const int d)
{
	ptrdiff_t ind = 0;
	for (int dim = 0; dim < d; dim++)
		ind += elem_per_dim->data[dim];
	return ind;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
