// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
///	\file

#include "mesh_connectivity.h"
#include "mesh_readers.h"

#include "Parameters.h"
#include "Macros.h"

#include "Multiarray.h"
#include "Vector.h"

static struct Multiarray_Vector_ui* compute_f_ve (const struct Mesh_Data*const mesh_data);

struct Mesh_Connectivity* mesh_connect (const struct Mesh_Data*const mesh_data)
{
	const size_t d = mesh_data->nodes->extents[1];
UNUSED(d);

	struct Multiarray_Vector_ui* f_ve = compute_f_ve(mesh_data);
UNUSED(f_ve);

	struct Mesh_Connectivity* mesh_connectivity = NULL;
	return mesh_connectivity;
}

// Static functions ************************************************************************************************* //

static struct Vector_ui* count_elements_per_dim (const struct const_Vector_ui*const elem_types)
{
	struct Vector_ui* count = constructor_empty_Vector_ui(DMAX+1);
	set_to_zero_Vector_ui(count);
	for (size_t i = 0; i < elem_types->extents[0]; i++) {
		const unsigned int elem_type = elem_types->data[i];

		switch (elem_type) {
		case POINT:
			count->data[0]++;
			break;
		case LINE:
			count->data[1]++;
			break;
		case TRI: case QUAD:
			count->data[2]++;
			break;
		case TET: case HEX: case WEDGE: case PYR:
			count->data[3]++;
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
		}
	}

	return count;
}

//static compute_index_

/*static size_t compute_number_of_volume_faces
	(const struct const_Vector_ui*const elem_types, const size_t*const elem_per_dim)
{
UNUSED(elem_types);
UNUSED(elem_per_dim);
return 0;
}*/

/// \brief Compute the list of (f)ace (ve)rtices for each face.
static struct Multiarray_Vector_ui* compute_f_ve (const struct Mesh_Data*const mesh_data)
{
	const struct const_Vector_ui*const            elem_types = mesh_data->elem_types;
	const struct const_Multiarray_Vector_ui*const node_nums  = mesh_data->node_nums;

	struct Vector_ui* elem_per_dim = count_elements_per_dim(elem_types);
	print_Vector_ui(elem_per_dim);
	EXIT_UNSUPPORTED;
UNUSED(node_nums);
}
