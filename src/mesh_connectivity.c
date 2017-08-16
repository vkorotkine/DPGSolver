// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
///	\file

#include "mesh_connectivity.h"
#include "mesh_readers.h"
#include "Intrusive.h"

#include <limits.h>

#include "Parameters.h"
#include "Macros.h"

#include "Multiarray.h"
#include "Vector.h"
#include "Element.h"

#include "const_cast.h"
#include "allocators.h"

// Static function declarations ************************************************************************************* //

/// \brief Container for locally computed \ref Mesh_Connectivity members.
struct Mesh_Connectivity_l {
	struct Multiarray_Vector_ui* v_to_v;  ///< Defined in \ref Mesh_Connectivity.
	struct Multiarray_Vector_ui* v_to_lf; ///< Defined in \ref Mesh_Connectivity.
};

/// \brief Container for local connectivity related information.
struct Conn_info {
	// Available from mesh_data:
	const unsigned int d;                 ///< The dimension.
	struct Vector_ui* elem_per_dim;       ///< The number of elements of each dimension.
	struct const_Vector_ui* volume_types; ///< Pointer to the first volume entry in \ref Mesh_Data::elem_types.
	struct Vector_ui* v_n_lf;             ///< The number of faces for each volume.

	// Computed here:
	struct Multiarray_Vector_ui* f_ve;     ///< Global face to vertex correspondence.
	struct Vector_ui*            ind_f_ve; ///< Indices of \ref f_ve after sorting.
};

/** \brief Constructor for \ref Conn_info.
 *	\return Standard. */
static struct Conn_info* constructor_Conn_info
	(const struct Mesh_Data*const mesh_data, ///< Standard.
	 const struct Intrusive_List* elements   ///< Standard.
	);

/// \brief Destructor for \ref Conn_info.
static void destructor_Conn_info
	(struct Conn_info* conn_info ///< Standard.
	);

/** \brief Compute the list of (f)ace (ve)rtices for each face.
 *	\return See brief. */
static void compute_f_ve
	(const struct Mesh_Data*const mesh_data, ///< Standard.
	 const struct Intrusive_List* elements,  ///< Standard.
	 struct Conn_info* conn_info             ///< The \ref Conn_info.
	);

/** \brief Compute the volume to (volume, local face) correspondence.
 *
 *	This function works by finding the volume and local face indices corresponding to the faces and then establishing
 *	connections by checking for adjacent matching face vertices in the sorted list of face vertices. If no connection
 *	is found, the entries are set to UINT_MAX.
 */
static void compute_v_to__v_lf
	(const struct Conn_info*const conn_info,      ///< The \ref Conn_info.
	 struct Mesh_Connectivity_l*const mesh_conn_l ///< The \ref Mesh_Connectivity_l.
	);

static void compute_bc_info
	()
{
//	const size_t d    = conn_info->d,
//	             n_bf = conn_info->elem_per_dim->data[d-1];
}

// Interface functions ********************************************************************************************** //

struct Mesh_Connectivity* mesh_connect (const struct Mesh_Data*const mesh_data, const struct Intrusive_List* elements)
{
	if (elements == NULL)
		EXIT_ADD_SUPPORT; // Add constructor. Don't forget to destruct.

	struct Mesh_Connectivity_l mesh_conn_l;
	struct Conn_info* conn_info = constructor_Conn_info(mesh_data,elements); // destructed

	compute_f_ve(mesh_data,elements,conn_info);

	compute_v_to__v_lf(conn_info,&mesh_conn_l);
print_Multiarray_Vector_ui(mesh_conn_l.v_to_v);
print_Multiarray_Vector_ui(mesh_conn_l.v_to_lf);

	compute_bc_info();


	destructor_Multiarray_Vector_ui(conn_info->f_ve);
	destructor_Vector_ui(conn_info->ind_f_ve);
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

	struct Vector_ui* v_n_lf = constructor_empty_Vector_ui(n_v);
	for (unsigned int v = 0; v < n_v; ++v) {
		const struct Element*const element = get_element_by_type(elements,volume_types->data[v]);
		v_n_lf->data[v] = element->n_f;
	}

	struct Conn_info* conn_info = malloc(sizeof *conn_info); // returned;
	const_cast_ui(&conn_info->d,d);
	conn_info->elem_per_dim = elem_per_dim;
	conn_info->volume_types = volume_types;
	conn_info->v_n_lf       = v_n_lf;

	return conn_info;
}

static void destructor_Conn_info (struct Conn_info* conn_info)
{
	destructor_Vector_ui(conn_info->elem_per_dim);
	free(conn_info);
}



/// \brief Container for periodic node information.
struct Periodic_Node {
	char   dir;           ///< The direction of the periodicity. Options: 'x', 'y', 'z'.
	double coord[DMAX-1]; ///< The node coordinates in directions other than the periodic one.
};

static bool check_bf_periodic
	(const char sm,        ///< Indicator for whether it is desired to check for 'M'aster or 'S'lave.
	 const unsigned int bc ///< The value of the boundary condition.
	)
{
	const unsigned int bc_base = bc % BC_STEP_SC;
	switch (sm) {
	case 'M':
		return ((bc_base == PERIODIC_XL) || (bc_base == PERIODIC_YL) || (bc_base == PERIODIC_ZL));
		break;
	case 'S':
		return ((bc_base == PERIODIC_XR) || (bc_base == PERIODIC_YR) || (bc_base == PERIODIC_ZR));
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
}

static size_t count_periodic_faces
	(const size_t ind_bf, const size_t n_bf, const struct const_Matrix_ui*const elem_tags)
{
	size_t count = 0;

	const size_t n_max = ind_bf+n_bf;
	for (size_t n = ind_bf; n < n_max; ++n) {
		if (check_bf_periodic('M',get_val_const_Matrix_ui(n,0,elem_tags)))
			++count;
	}
	return count;
}

/// \brief Container for periodic face information.
struct Periodic_Face_Info {
	struct Periodic_Node** per_nodes;       ///< The list of \ref Periodic_Node entities.
	struct Vector_ui**     bface_node_nums; ///< The boundary face node numbers.
	struct Vector_ui*      ind_sorted;      ///< The indices corresponding to the sorted boundary face node numbers.
};

static struct Periodic_Face_Info* constructor_Periodic_Face_Info (size_t const n_pf)
{
	struct Periodic_Face_Info* p_f_info = malloc(sizeof *p_f_info); // returned

	p_f_info->per_nodes       = malloc(n_pf * sizeof *(p_f_info->per_nodes));       // destructed
	p_f_info->bface_node_nums = malloc(n_pf * sizeof *(p_f_info->bface_node_nums)); // destructed

	return p_f_info;
}

static void destructor_Periodic_Face_Info (struct Periodic_Face_Info* p_f_info)
{
	free(p_f_info->per_nodes);
	free(p_f_info->bface_node_nums);

	free(p_f_info);
}

static void set_p_nodes_and_nums
	(struct Periodic_Face_Info* p_info[2], const size_t ind_bf, const size_t n_bf,
	 const struct Mesh_Data*const mesh_data)
{
	const struct const_Matrix_d*const             nodes     = mesh_data->nodes;
	const struct const_Matrix_ui*const            elem_tags = mesh_data->elem_tags;
	const struct const_Multiarray_Vector_ui*const node_nums = mesh_data->node_nums;

	const size_t n_max = ind_bf+n_bf;
	for (size_t n = ind_bf; n < n_max; ++n) {
		if (check_bf_periodic('M',get_val_const_Matrix_ui(n,0,elem_tags)))
			; // Fill master [0]
		else if (check_bf_periodic('S',get_val_const_Matrix_ui(n,0,elem_tags)))
			; // Fill slave [1]
	}
}

/** \brief Correct the face vertex correspondence if the mesh is periodic.
 *
 *	The correction modifies the vertices of the "slave" faces to those of the "master" faces. The correspondence is
 *	established by checking the equivalence of the mesh node coordinates in the appropriate directions.
 */
void correct_f_ve_for_periodic
	(struct Multiarray_Vector_ui*const f_ve, ///< Defined in \ref Conn_info.
	 const struct Mesh_Data*const mesh_data, ///< Defined in \ref Mesh_Data.
	 const struct Conn_info*const conn_info  ///< The \ref Conn_info.
	)
{
	if (mesh_data->periodic_corr == NULL)
		return;

	// 1) Make lists of master and slave [direction, centroid] containers (sorted).
	// 2) Make lists of master and slave boundary face vertices (sorted Vector data only).
	// 2.1) Find face indices of master faces (binary search).
	// 2.2) Replace face indices of slave faces (binary search) with those of corresponding master faces.

	const size_t d      = conn_info->d,
	             ind_bf = get_first_volume_index(conn_info->elem_per_dim,d-1),
	             n_bf   = conn_info->elem_per_dim->data[d-1];

	const size_t n_pf = count_periodic_faces(ind_bf,n_bf,mesh_data->elem_tags);

	const size_t n_ms = 2;
	struct Periodic_Face_Info* p_info[n_ms];

	for (size_t i = 0; i < n_ms; i++)
		p_info[i] = constructor_Periodic_Face_Info(n_pf);



	for (size_t i = 0; i < n_ms; i++)
		destructor_Periodic_Face_Info(p_info[i]);


UNUSED(p_info);

printf("%zu\n",n_pf);
UNUSED(n_pf);

print_const_Matrix_d(mesh_data->nodes);
print_const_Matrix_ui(mesh_data->elem_tags);
print_const_Multiarray_Vector_ui(mesh_data->node_nums);
//print_const_Matrix_ui(mesh_data->periodic_corr);
//print_Multiarray_Vector_ui(conn_info->f_ve);
	const struct const_Vector_ui*const*const bface_nums = &mesh_data->node_nums->data[ind_bf];
	for (size_t bf = 0; bf < n_bf; ++bf) {
	}
// void * std::bsearch (const void *key, const void *array, size_t count, size_t size, comparison_fn_t compare)
UNUSED(ind_bf);
UNUSED(bface_nums);


UNUSED(f_ve);
EXIT_UNSUPPORTED;
}


static void compute_f_ve
	(const struct Mesh_Data*const mesh_data, const struct Intrusive_List* elements, struct Conn_info* conn_info)
{
	const unsigned int d     = conn_info->d,
	                   ind_v = get_first_volume_index(conn_info->elem_per_dim,d),
	                   n_v   = conn_info->elem_per_dim->data[d];

	struct const_Vector_ui* volume_types = conn_info->volume_types;

	const unsigned int sum_n_f = compute_sum_n_f(elements,volume_types);
	struct Multiarray_Vector_ui* f_ve = constructor_empty_Multiarray_Vector_ui(1,sum_n_f); // returned

	const struct const_Vector_ui*const*const volume_nums = &mesh_data->node_nums->data[ind_v];
	for (unsigned int v = 0, ind_f = 0; v < n_v; ++v) {
		const struct Element*const element = get_element_by_type(elements,volume_types->data[v]);
		for (size_t f = 0, f_max = conn_info->v_n_lf->data[v]; f < f_max; ++f) {
			const struct const_Vector_ui*const f_ve_f = element->f_ve->data[f];
			const size_t n_n = f_ve_f->extents[0];

			struct Vector_ui* f_ve_curr = f_ve->data[ind_f];
			reserve_Vector_ui(f_ve_curr,n_n);
			for (size_t n = 0; n < n_n; ++n)
				f_ve_curr->data[n] = volume_nums[v]->data[f_ve_f->data[n]];
			++ind_f;
		}
	}

	conn_info->f_ve     = f_ve; // keep
	conn_info->ind_f_ve = sort_Multiarray_Vector_ui(f_ve,true); // keep

	correct_f_ve_for_periodic(f_ve,mesh_data,conn_info);

}

static void compute_v_to__v_lf (const struct Conn_info*const conn_info, struct Mesh_Connectivity_l*const mesh_conn_l)
{
	struct Vector_ui* ind_f_ve_V = conn_info->ind_f_ve;

	const size_t d   = conn_info->d,
	             n_v = conn_info->elem_per_dim->data[d],
	             n_f = ind_f_ve_V->extents[0];

	struct Vector_ui* v_n_lf = conn_info->v_n_lf;

	// Store global volume and local face indices corresponding to each global face (reordered).
	unsigned int* ind_v_ui  = mallocator(UINT_T,1,n_f); // moved
	unsigned int* ind_lf_ui = mallocator(UINT_T,1,n_f); // moved

	for (size_t ind_vf = 0, v = 0; v < n_v; ++v) {
		const size_t lf_max = v_n_lf->data[v];
		for (size_t lf = 0; lf < lf_max; ++lf) {
			ind_v_ui[ind_vf]  = v;
			ind_lf_ui[ind_vf] = lf;
			++ind_vf;
		}
	}

	struct Vector_ui* ind_v_V  = constructor_move_Vector_ui_ui(n_f,true,ind_v_ui);  // destructed
	struct Vector_ui* ind_lf_V = constructor_move_Vector_ui_ui(n_f,true,ind_lf_ui); // destructed

	reorder_Vector_ui(ind_v_V,ind_f_ve_V->data);
	reorder_Vector_ui(ind_lf_V,ind_f_ve_V->data);

	// Compute v_to_v and v_to_lf
	unsigned int* v_to_v_ui  = mallocator(UINT_T,1,n_f); // moved
	unsigned int* v_to_lf_ui = mallocator(UINT_T,1,n_f); // moved

	struct Multiarray_Vector_ui* f_ve = conn_info->f_ve;
	const unsigned int*const ind_f_ve_ui = ind_f_ve_V->data;

	bool found_match = false;
	for (size_t f = 0; f < n_f; f++) {
		const size_t ind_0 = ind_f_ve_ui[f];
		if ((f+1) < n_f && check_equal_Vector_ui(f_ve->data[f],f_ve->data[f+1])) {
			const size_t ind_1 = ind_f_ve_ui[f+1];

			v_to_v_ui[ind_0]  = ind_v_ui[f+1];
			v_to_lf_ui[ind_0] = ind_lf_ui[f+1];
			v_to_v_ui[ind_1]  = ind_v_ui[f];
			v_to_lf_ui[ind_1] = ind_lf_ui[f];

			found_match = true;
		} else {
			if (!found_match) {
				v_to_v_ui[ind_0]  = UINT_MAX;
				v_to_lf_ui[ind_0] = UINT_MAX;
			}
			found_match = false;
		}
	}
	destructor_Vector_ui(ind_v_V);
	destructor_Vector_ui(ind_lf_V);

	mesh_conn_l->v_to_v  = constructor_copy_Multiarray_Vector_ui(v_to_v_ui,conn_info->v_n_lf->data,1,n_v);  // keep
	mesh_conn_l->v_to_lf = constructor_copy_Multiarray_Vector_ui(v_to_lf_ui,conn_info->v_n_lf->data,1,n_v); // keep
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
