// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
///	\file

#include "mesh_connectivity.h"
#include "mesh_readers.h"
#include "Intrusive.h"

#include "Parameters.h"
#include "Macros.h"

#include "Multiarray.h"
#include "Vector.h"
#include "Element.h"
#include "const_cast.h"

// Static function declarations ************************************************************************************* //

/// \brief Container for local connectivity related information.
struct Conn_info {
	const unsigned int d;                 ///< The dimension.
	struct Vector_ui* elem_per_dim;       ///< The number of elements of each dimension.
	struct const_Vector_ui* volume_types; ///< Pointer to the first volume entry in \ref Mesh_Data::elem_types.
	struct Vector_ui* v_n_f;              ///< The number of faces for each volume.
};

/// \brief Constructor for \ref Conn_info.
static struct Conn_info* constructor_Conn_info
	(const struct Mesh_Data*const mesh_data, ///< Standard.
	 const struct Intrusive_List* elements   ///< Standard.
	);

/// \brief Destructor for \ref Conn_info.
static void destructor_Conn_info
	(struct Conn_info* conn_info ///< Standard.
	);

/// \brief Compute the list of (f)ace (ve)rtices for each face.
static struct Multiarray_Vector_ui* compute_f_ve
	(const struct Mesh_Data*const mesh_data, ///< Standard.
	 const struct Intrusive_List* elements,  ///< Standard.
	 const struct Conn_info* conn_info       ///< The \ref Conn_info.
	);

// Interface functions ********************************************************************************************** //

struct Mesh_Connectivity* mesh_connect (const struct Mesh_Data*const mesh_data, const struct Intrusive_List* elements)
{
	if (elements == NULL)
		EXIT_ADD_SUPPORT; // Add constructor. Don't forget to destruct.

	struct Conn_info* conn_info = constructor_Conn_info(mesh_data,elements); // destructed


	struct Multiarray_Vector_ui* f_ve = compute_f_ve(mesh_data,elements,conn_info);
//print_Multiarray_Vector_ui(f_ve);
	struct Vector_ui* indices_f_ve = sort_Multiarray_Vector_ui(f_ve,true); // destructed
print_Multiarray_Vector_ui(f_ve);
print_Vector_ui(indices_f_ve);


	const unsigned int d   = conn_info->d,
	                   n_v = conn_info->elem_per_dim->data[d],
	                   n_f = sum_Vector_ui(conn_info->v_n_f);

	struct Vector_ui* v_n_f = conn_info->v_n_f;

	unsigned int v_to_v_ui[n_f];
	unsigned int v_to_f_ui[n_f];

	size_t ind_vf = 0;
	for (size_t v = 0; v < n_v; ++v) {
		for (size_t f = 0, f_max = v_n_f->data[v]; f < f_max; ++f) {
			v_to_v_ui[ind_vf] = v;
			v_to_f_ui[ind_vf] = f;
			++ind_vf;
		}
	}

	struct Multiarray_Vector_ui* v_to_v =
		constructor_copy_Multiarray_Vector_ui(v_to_v_ui,conn_info->v_n_f->data,1,n_v); // keep
	struct Multiarray_Vector_ui* v_to_f =
		constructor_copy_Multiarray_Vector_ui(v_to_f_ui,conn_info->v_n_f->data,1,n_v); // keep
print_Multiarray_Vector_ui(v_to_v);
print_Multiarray_Vector_ui(v_to_f);

	destructor_Vector_ui(indices_f_ve);
	destructor_Conn_info(conn_info);

	struct Mesh_Connectivity* mesh_connectivity = NULL;
	return mesh_connectivity;
}

// Static functions ************************************************************************************************* //

/** \brief See return.
 *	\return The number of elements of each dimension.
 */
static struct Vector_ui* count_elements_per_dim
	(const struct const_Vector_ui*const elem_types ///< Defined in \ref Conn_info.
	);

/** \brief See return.
 *	\return The index of the first volume.
 */
static size_t get_first_volume_index
	(const struct Vector_ui*const elem_per_dim, ///< Defined in \ref Conn_info.
	 const unsigned int d                       ///< Defined in \ref Conn_info.
	);

/** \brief See return.
 *	\return The sum of the number of faces of all volumes.
 */
static size_t compute_sum_n_f
	(const struct Intrusive_List* elements,          ///< Standard.
	 const struct const_Vector_ui*const volume_types ///< Defined in \ref Conn_info.
	);

static struct Conn_info* constructor_Conn_info
	(const struct Mesh_Data*const mesh_data, const struct Intrusive_List* elements)
{
	const size_t d = mesh_data->nodes->extents[1];

	struct Vector_ui* elem_per_dim = count_elements_per_dim(mesh_data->elem_types);

	const unsigned int n_v   = elem_per_dim->data[d],
	                   ind_v = get_first_volume_index(elem_per_dim,d);

	struct const_Vector_ui* volume_types =
		constructor_move_const_Vector_ui_ui(n_v,false,&mesh_data->elem_types->data[ind_v]);

	struct Vector_ui* v_n_f = constructor_empty_Vector_ui(n_v);
	for (unsigned int v = 0; v < n_v; ++v) {
		const struct Element*const element = get_element_by_type(elements,volume_types->data[v]);
		v_n_f->data[v] = element->n_f;
	}

	struct Conn_info* conn_info = malloc(sizeof *conn_info); // returned;
	const_cast_ui(&conn_info->d,d);
	conn_info->elem_per_dim = elem_per_dim;
	conn_info->volume_types = volume_types;
	conn_info->v_n_f        = v_n_f;

	return conn_info;
}

static void destructor_Conn_info (struct Conn_info* conn_info)
{
	destructor_Vector_ui(conn_info->elem_per_dim);
	free(conn_info);
}

static struct Multiarray_Vector_ui* compute_f_ve
	(const struct Mesh_Data*const mesh_data, const struct Intrusive_List* elements, const struct Conn_info* conn_info)
{
	const unsigned int d     = conn_info->d,
	                   ind_v = get_first_volume_index(conn_info->elem_per_dim,d),
	                   num_v = conn_info->elem_per_dim->data[d];

	struct const_Vector_ui* volume_types = conn_info->volume_types;

	const unsigned int sum_n_f = compute_sum_n_f(elements,volume_types);
	struct Multiarray_Vector_ui* f_ve = constructor_empty_Multiarray_Vector_ui(1,sum_n_f); // returned

	const struct const_Vector_ui*const*const volume_nums = &mesh_data->node_nums->data[ind_v];
	for (unsigned int v = 0, ind_f = 0; v < num_v; ++v) {
		const struct Element*const element = get_element_by_type(elements,volume_types->data[v]);
		for (size_t f = 0, f_max = conn_info->v_n_f->data[v]; f < f_max; ++f) {
			const struct const_Vector_ui*const f_ve_f = element->f_ve->data[f];
			const size_t n_n = f_ve_f->extents[0];

			struct Vector_ui* f_ve_curr = f_ve->data[ind_f];
			reserve_Vector_ui(f_ve_curr,n_n);
			for (size_t n = 0; n < n_n; ++n)
				f_ve_curr->data[n] = volume_nums[v]->data[f_ve_f->data[n]];
			++ind_f;
		}
	}

	return f_ve;
}

static struct Vector_ui* count_elements_per_dim (const struct const_Vector_ui*const elem_types)
{
	struct Vector_ui* count = constructor_empty_Vector_ui(DMAX+1); // returned
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

static size_t get_first_volume_index (const struct Vector_ui*const elem_per_dim, const unsigned int d)
{
	size_t ind = 0;
	for (unsigned int dim = 0; dim < d; dim++)
		ind += elem_per_dim->data[dim];
	return ind;
}

static size_t compute_sum_n_f (const struct Intrusive_List* elements, const struct const_Vector_ui*const volume_types)
{
	size_t sum_n_f = 0;
	for (size_t v = 0, v_max = volume_types->extents[0]; v < v_max; ++v) {
		const struct Element*const element = get_element_by_type(elements,volume_types->data[v]);
		sum_n_f += element->n_f;
	}
	return sum_n_f;
}
