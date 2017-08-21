// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
///	\file

#include "mesh_connectivity.h"
#include "mesh_readers.h"
#include "Intrusive.h"

#include <limits.h>

#include "constants_core.h"
#include "constants_elements.h"
#include "constants_bc.h"
#include "Macros.h"

#include "Multiarray.h"
#include "Vector.h"
#include "Element.h"

#include "mesh_periodic.h"
#include "const_cast.h"
#include "allocators.h"

// Static function declarations ************************************************************************************* //

/// \brief Container for locally computed \ref Mesh_Connectivity members.
struct Mesh_Connectivity_l {
	struct Multiarray_Vector_ui* v_to_v;  ///< Local version of \ref Mesh_Connectivity::v_to_v.
	struct Multiarray_Vector_ui* v_to_lf; ///< Local version of \ref Mesh_Connectivity::v_to_lf.
};

/** \brief Constructor for \ref Conn_info.
 *	\return Standard. */
static struct Conn_info* constructor_Conn_info
	(const struct Mesh_Data*const mesh_data,          ///< Standard.
	 const struct const_Intrusive_List*const elements ///< Standard.
	);

/// \brief Destructor for \ref Conn_info.
static void destructor_Conn_info
	(struct Conn_info* conn_info ///< Standard.
	);

/** \brief Compute the list of (f)ace (ve)rtices for each face.
 *	\return See brief. */
static void compute_f_ve
	(const struct Mesh_Data*const mesh_data,           ///< Standard.
	 const struct const_Intrusive_List*const elements, ///< Standard.
	 struct Conn_info* conn_info                       ///< The \ref Conn_info.
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

/** \brief Add the boundary condition information to \ref v_to_lf.
 *
 *	This function replaces all invalid entries in the volume to local face container with the number corresponding to
 *	the appropriate boundary condition.
 */
static void add_bc_info
	(const struct Mesh_Data*const mesh_data,      ///< Standard.
	 const struct Conn_info*const conn_info,      ///< The \ref Conn_info.
	 struct Mesh_Connectivity_l*const mesh_conn_l ///< The \ref Mesh_Connectivity_l.
	);

/** \brief Constructor for the \ref Mesh_Connectivity.
 *	\return Standard. */
static struct Mesh_Connectivity* constructor_Mesh_Connectivity
	(const struct Mesh_Connectivity_l*const mesh_conn_l ///< The \ref Mesh_Connectivity_l.
	);

// Interface functions ********************************************************************************************** //

struct Mesh_Connectivity* mesh_connect
	(const struct Mesh_Data*const mesh_data, const struct const_Intrusive_List* elements)
{
	if (elements == NULL)
		EXIT_ADD_SUPPORT; // Add constructor. Don't forget to destruct.

	struct Mesh_Connectivity_l mesh_conn_l;
	struct Conn_info* conn_info = constructor_Conn_info(mesh_data,elements); // destructed

	compute_f_ve(mesh_data,elements,conn_info);
	compute_v_to__v_lf(conn_info,&mesh_conn_l);
	add_bc_info(mesh_data,conn_info,&mesh_conn_l);

	destructor_Multiarray_Vector_ui(conn_info->f_ve);
	destructor_Vector_ui(conn_info->ind_f_ve);
	destructor_Conn_info(conn_info);

	struct Mesh_Connectivity* mesh_connectivity = constructor_Mesh_Connectivity(&mesh_conn_l); // keep

print_const_Multiarray_Vector_ui(mesh_connectivity->v_to_v);

	return mesh_connectivity;
}

void destructor_Mesh_Connectivity (struct Mesh_Connectivity* mesh_conn)
{
	destructor_Multiarray_Vector_ui((struct Multiarray_Vector_ui*)mesh_conn->v_to_v);
	destructor_Multiarray_Vector_ui((struct Multiarray_Vector_ui*)mesh_conn->v_to_lf);

	free(mesh_conn);
}

size_t get_first_volume_index (const struct Vector_ui*const elem_per_dim, const unsigned int d)
{
	size_t ind = 0;
	for (unsigned int dim = 0; dim < d; dim++)
		ind += elem_per_dim->data[dim];
	return ind;
}

void set_f_node_nums (struct Vector_ui**const f_node_nums, const struct const_Vector_ui*const node_nums)
{
	*f_node_nums = constructor_copy_Vector_ui_ui(node_nums->extents[0],node_nums->data);
	sort_Vector_ui(*f_node_nums);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Container for the list of boundary face information.
struct Boundary_Face_Info {
	size_t n_pfe;                   ///< The number of physical face elements.
	size_t n_bf;                    /**< The number of boundary faces (n_bf = n_pfe - n_pf, where n_pf: number of
	                                 *   periodic faces). */
	struct Boundary_Face** b_faces; ///< The list of \ref Boundary_Face entities.
};

/// \brief Container for boundary face information.
struct Boundary_Face {
	unsigned int bc;             ///< The value of the boundary condition.
	struct Vector_ui* node_nums; ///< The node numbers of the face vertices.
};

/** \brief Constructor for \ref Boundary_Face_Info.
 *	\return Standard. */
static struct Boundary_Face_Info* constructor_Boundary_Face_Info
	(size_t const n_pfe, ///< \ref Boundary_Face_Info::n_pfe.
	 size_t const n_bf   ///< \ref Boundary_Face_Info::n_bf.
	);

/// \brief Destructor for \ref Boundary_Face_Info.
static void destructor_Boundary_Face_Info
	(struct Boundary_Face_Info* bf_info ///< Standard.
	);

/** \brief See return.
 *	\return The number of elements of each dimension.
 */
static struct Vector_ui* count_elements_per_dim
	(const struct const_Vector_ui*const elem_types ///< Defined in \ref Conn_info.
	);

/** \brief See return.
 *	\return The sum of the number of faces of all volumes.
 */
static size_t compute_sum_n_f
	(const struct const_Intrusive_List*const elements, ///< Standard.
	 const struct const_Vector_ui*const volume_types   ///< Defined in \ref Conn_info.
	);

/// \brief Set boundary face info for all entries in the list.
static void set_bf_info
	(struct Boundary_Face_Info* bf_info,    ///< The \ref Boundary_Face_Info.
	 const size_t ind_pfe,                  ///< Index of the first physical face element in the mesh element list.
	 const struct Mesh_Data*const mesh_data ///< The \ref Mesh_Data.
	);

/** \brief Count the number of boundary faces.
 *	\return See brief. */
static size_t count_boundary_faces
	(const size_t ind_pfe,                        /**< Index of the first physical face element in the mesh element
	                                               *   list. */
	 const size_t n_pfe,                          ///< The number of physical face elements.
	 const struct const_Matrix_ui*const elem_tags ///< \ref Mesh_Data::elem_tags.
	);

/** \brief Reorder the boundary face entries in the list according to the node numbering.
 *	\warning This is not currently done in place.
 */
static void reorder_b_faces
	(struct Boundary_Face**const b_faces, ///< \ref Boundary_Face_Info::b_faces.
	 struct Vector_ui* ordering           ///< The new ordering
	);

/** \brief Comparison function for std::bsearch between \ref Boundary_Face\*\* `a` and `b`.
 *	\return The \ref cmp_Vector_ui of the `node_nums` of `a` and `b`.
 *
 *	\note Input Vectors must be have sorted data.
 */
static int cmp_Boundary_Face
	(const void *a, ///< Variable 1.
	 const void *b  ///< Variable 2.
	);

/// \brief Update v_to_lf replacing unused boundary entries with the boundary values.
static void update_v_to_lf_bc
	(struct Multiarray_Vector_ui*const v_to_lf, ///< \ref Mesh_Connectivity_l::v_to_lf.
	 const unsigned int*const v_to_lf_ui        ///< The `v_to_lf` data with bc information included.
	);

static struct Conn_info* constructor_Conn_info
	(const struct Mesh_Data*const mesh_data, const struct const_Intrusive_List*const elements)
{
	const size_t d = mesh_data->nodes->extents[1];

	struct Vector_ui* elem_per_dim = count_elements_per_dim(mesh_data->elem_types);

	const unsigned int n_v   = elem_per_dim->data[d],
	                   ind_v = get_first_volume_index(elem_per_dim,d);

	struct const_Vector_ui* volume_types =
		constructor_move_const_Vector_ui_ui(n_v,false,&mesh_data->elem_types->data[ind_v]);

	struct Vector_ui* v_n_lf = constructor_empty_Vector_ui(n_v);
	for (unsigned int v = 0; v < n_v; ++v) {
		const struct const_Element*const element = get_element_by_type(elements,volume_types->data[v]);
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
	destructor_Vector_ui((struct Vector_ui*)conn_info->volume_types);
	destructor_Vector_ui(conn_info->v_n_lf);
	free(conn_info);
}

static void compute_f_ve
	(const struct Mesh_Data*const mesh_data, const struct const_Intrusive_List*const elements,
	 struct Conn_info* conn_info)
{
	const unsigned int d     = conn_info->d,
	                   ind_v = get_first_volume_index(conn_info->elem_per_dim,d),
	                   n_v   = conn_info->elem_per_dim->data[d];

	struct const_Vector_ui* volume_types = conn_info->volume_types;

	const unsigned int sum_n_f = compute_sum_n_f(elements,volume_types);
	struct Multiarray_Vector_ui* f_ve = constructor_empty_Multiarray_Vector_ui(1,sum_n_f); // returned

	const struct const_Vector_ui*const*const volume_nums = &mesh_data->node_nums->data[ind_v];
	for (unsigned int v = 0, ind_f = 0; v < n_v; ++v) {
		const struct const_Element*const element = get_element_by_type(elements,volume_types->data[v]);
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

	correct_f_ve_for_periodic(mesh_data,conn_info);
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
	unsigned int* v_to_v_ui  = mallocator(UINT_T,1,n_f); // free
	unsigned int* v_to_lf_ui = mallocator(UINT_T,1,n_f); // free

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

	free(v_to_v_ui);
	free(v_to_lf_ui);
}

static void add_bc_info
	(const struct Mesh_Data*const mesh_data, const struct Conn_info*const conn_info,
	 struct Mesh_Connectivity_l*const mesh_conn_l)
{
	const size_t d       = conn_info->d,
	             ind_pfe = get_first_volume_index(conn_info->elem_per_dim,d-1),
	             n_pfe   = conn_info->elem_per_dim->data[d-1],
	             n_bf    = count_boundary_faces(ind_pfe,n_pfe,mesh_data->elem_tags);

	if (n_bf == 0)
		return;

	// Set the boundary face info and sort it by node_nums.
	struct Boundary_Face_Info* bf_info = constructor_Boundary_Face_Info(n_pfe,n_bf); // destructed
	set_bf_info(bf_info,ind_pfe,mesh_data);

	// Copy the pointers to the node_nums into a Multiarray_Vector_ui (for sorting).
	struct Multiarray_Vector_ui* bf_ve = constructor_empty_Multiarray_Vector_ui(1,n_bf); // destructed
	bf_ve->owns_data = false;
	for (size_t i = 0; i < n_bf; ++i) {
		destructor_Vector_ui(bf_ve->data[i]);
		bf_ve->data[i] = bf_info->b_faces[i]->node_nums;
	}

	struct Vector_ui* ind_bf_ve = sort_Multiarray_Vector_ui(bf_ve,true); // destructed
	destructor_Multiarray_Vector_ui(bf_ve);

	reorder_b_faces(bf_info->b_faces,ind_bf_ve);
	destructor_Vector_ui(ind_bf_ve);

	// Set the unused entries in v_to_lf to the values of the associated boundary conditions.
	struct Vector_ui*const v_to_lf_V = collapse_Multiarray_Vector_ui(mesh_conn_l->v_to_lf); // destructed

	unsigned int*const v_to_lf_ui  = v_to_lf_V->data,
	            *const ind_f_ve_ui = conn_info->ind_f_ve->data;

	struct Vector_ui*const*const f_ve_V = conn_info->f_ve->data;

	struct Boundary_Face*const f_curr = malloc(sizeof *f_curr); // free
	f_curr->bc = UINT_MAX;

	const size_t n_f = v_to_lf_V->extents[0];
	for (size_t n = 0; n < n_f; ++n) {
		const size_t f = ind_f_ve_ui[n];
		if (v_to_lf_ui[f] != UINT_MAX)
			continue;

		f_curr->node_nums = f_ve_V[n];
		struct Boundary_Face*const*const bf_curr =
			bsearch(&f_curr,bf_info->b_faces,n_bf,sizeof(bf_info->b_faces[0]),cmp_Boundary_Face);

		v_to_lf_ui[f] = (*bf_curr)->bc;
	}
	free(f_curr);

	update_v_to_lf_bc(mesh_conn_l->v_to_lf,v_to_lf_ui);

	destructor_Vector_ui(v_to_lf_V);
	destructor_Boundary_Face_Info(bf_info);
}

static struct Mesh_Connectivity* constructor_Mesh_Connectivity (const struct Mesh_Connectivity_l*const mesh_conn_l)
{
	struct Mesh_Connectivity* mesh_conn = malloc(sizeof *mesh_conn); // returned

	const_constructor_move_Multiarray_Vector_ui(&mesh_conn->v_to_v,mesh_conn_l->v_to_v);
	const_constructor_move_Multiarray_Vector_ui(&mesh_conn->v_to_lf,mesh_conn_l->v_to_lf);

	return mesh_conn;
}

// Level 1 ********************************************************************************************************** //

/** \brief Constructor for \ref Boundary_Face.
 *	\return Standard. */
static struct Boundary_Face* constructor_Boundary_Face ( );

/** \brief Copy constructor for \ref Boundary_Face_Info::b_faces.
 *	\return Standard. */
static struct Boundary_Face** constructor_b_faces
	(const size_t n_bf,               ///< \ref Boundary_Face_Info::n_bf.
	 struct Boundary_Face**const src  ///< The source \ref Boundary_Face_Info::b_faces data.
	);

/// Destructor for \ref Boundary_Face_Info::b_faces.
static void destructor_b_faces
	(const size_t n_bf,               ///< \ref Boundary_Face_Info::n_bf.
	 struct Boundary_Face**const src  ///< Standard.
	);

/// \brief Destructor for \ref Boundary_Face.
static void destructor_Boundary_Face
	(struct Boundary_Face* bf ///< Standard.
	);

/** \brief Check if the physical face element is a boundary which is not periodic.
 *	\return See brief. */
static bool check_pfe_boundary
	(const unsigned int bc ///< The value of the boundary condition.
	);

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

static size_t compute_sum_n_f
	(const struct const_Intrusive_List*const elements, const struct const_Vector_ui*const volume_types)
{
	size_t sum_n_f = 0;
	for (size_t v = 0, v_max = volume_types->extents[0]; v < v_max; ++v) {
		const struct const_Element*const element = get_element_by_type(elements,volume_types->data[v]);
		sum_n_f += element->n_f;
	}
	return sum_n_f;
}

static struct Boundary_Face_Info* constructor_Boundary_Face_Info (const size_t n_pfe, const size_t n_bf)
{
	struct Boundary_Face_Info* bf_info = malloc(sizeof *bf_info); // returned

	struct Boundary_Face** b_faces = malloc(n_bf * sizeof *b_faces); // keep
	for (size_t i = 0; i < n_bf; ++i)
		b_faces[i] = constructor_Boundary_Face(); // destructed

	bf_info->n_pfe   = n_pfe;
	bf_info->n_bf    = n_bf;
	bf_info->b_faces = b_faces;

	return bf_info;
}

static void destructor_Boundary_Face_Info (struct Boundary_Face_Info* bf_info)
{
	const size_t n_bf = bf_info->n_bf;
	for (size_t i = 0; i < n_bf; ++i)
		destructor_Boundary_Face(bf_info->b_faces[i]);
	free(bf_info->b_faces);

	free(bf_info);
}

static void set_bf_info
	(struct Boundary_Face_Info* bf_info, const size_t ind_pfe, const struct Mesh_Data*const mesh_data)
{
	const struct const_Matrix_ui*const            elem_tags = mesh_data->elem_tags;
	const struct const_Multiarray_Vector_ui*const node_nums = mesh_data->node_nums;

	size_t count_bf = 0;

	const size_t n_max = ind_pfe+bf_info->n_bf;
	for (size_t n = ind_pfe; n < n_max; ++n) {
		const unsigned int bc = get_val_const_Matrix_ui(n,0,elem_tags);

		if (!check_pfe_boundary(bc))
			continue;

		struct Boundary_Face*const bf = bf_info->b_faces[count_bf];
		bf->bc = bc;

		set_f_node_nums(&bf->node_nums,node_nums->data[n]);

		++count_bf;
	}

	if (count_bf != bf_info->n_bf)
		EXIT_ERROR("Did not find the correct number of boundary face entities");
}

static size_t count_boundary_faces
	(const size_t ind_pfe, const size_t n_pfe, const struct const_Matrix_ui*const elem_tags)
{
	size_t count = 0;

	const size_t n_max = ind_pfe+n_pfe;
	for (size_t n = ind_pfe; n < n_max; ++n) {
		if (check_pfe_boundary(get_val_const_Matrix_ui(n,0,elem_tags)))
			++count;
	}
	return count;
}

static void reorder_b_faces (struct Boundary_Face**const b_faces, struct Vector_ui* ordering)
{
	const size_t n_bf = compute_size(1,ordering->extents);

	struct Boundary_Face** b_faces_cpy = constructor_b_faces(n_bf,b_faces); // destructed

	const unsigned int*const ordering_ui = ordering->data;
	for (size_t i = 0; i < n_bf; ++i) {
		const unsigned int ind_bf = ordering_ui[i];
		b_faces_cpy[i]->bc = b_faces[ind_bf]->bc;
		set_to_data_Vector_ui(b_faces_cpy[i]->node_nums,b_faces[ind_bf]->node_nums->data);
	}

	for (size_t i = 0; i < n_bf; ++i) {
		b_faces[i]->bc = b_faces_cpy[i]->bc;
		set_to_data_Vector_ui(b_faces[i]->node_nums,b_faces_cpy[i]->node_nums->data);
	}

	destructor_b_faces(n_bf,b_faces_cpy);
}

static int cmp_Boundary_Face (const void *a, const void *b)
{
	const struct Boundary_Face*const*const ia = (const struct Boundary_Face*const*const) a,
	                          *const*const ib = (const struct Boundary_Face*const*const) b;

	return cmp_Vector_ui(&(*ia)->node_nums,&(*ib)->node_nums);
}

static void update_v_to_lf_bc (struct Multiarray_Vector_ui*const v_to_lf, const unsigned int*const v_to_lf_ui)
{
	size_t count = 0;
	const size_t i_max = compute_size(v_to_lf->order,v_to_lf->extents);
	for (size_t i = 0; i < i_max; ++i) {
		struct Vector_ui* v_to_lf_curr = v_to_lf->data[i];
		unsigned int*const data = v_to_lf_curr->data;
		const size_t j_max = compute_size(1,v_to_lf_curr->extents);
		for (size_t j = 0; j < j_max; ++j) {
			if (data[j] == UINT_MAX)
				data[j] = v_to_lf_ui[count];
			count++;
		}
	}
}

// Level 2 ********************************************************************************************************** //

/** \brief Constructor for \ref Boundary_Face where the memory is allocated for the node_nums but not set.
 *	\return Standard. */
static struct Boundary_Face* constructor_empty_Boundary_Face
	(struct Boundary_Face* src ///< The source data.
	);

static struct Boundary_Face* constructor_Boundary_Face ( )
{
	struct Boundary_Face* bf = malloc(sizeof *bf); // returned

	return bf;
}

static void destructor_Boundary_Face (struct Boundary_Face* bf)
{
	destructor_Vector_ui(bf->node_nums);
	free(bf);
}

static struct Boundary_Face** constructor_b_faces (const size_t n_bf, struct Boundary_Face**const src)
{
	struct Boundary_Face** dest = malloc(n_bf * sizeof *dest); // free
	for (size_t i = 0; i < n_bf; ++i)
		dest[i] = constructor_empty_Boundary_Face(src[i]); // destructed

	return dest;
}

static void destructor_b_faces (const size_t n_bf, struct Boundary_Face** src)
{
	for (size_t i = 0; i < n_bf; ++i)
		destructor_Boundary_Face(src[i]);
	free(src);
}

static bool check_pfe_boundary (const unsigned int bc)
{
	const unsigned int bc_base = bc % BC_STEP_SC;
	switch (bc_base) {
		case BC_INFLOW: // Advection
		case BC_OUTFLOW:
		case BC_DIRICHLET: // Poisson
		case BC_NEUMANN:
		case BC_RIEMANN: // Euler
		case BC_SLIPWALL:
		case BC_BACKPRESSURE:
		case BC_TOTAL_TP:
		case BC_SUPERSONIC_IN:
		case BC_SUPERSONIC_OUT:
		case BC_NOSLIP_T: // Navier-Stokes
		case BC_NOSLIP_ADIABATIC:
			return true;
			break;
		default:
			return false;
			break;
	}
}

// Level 3 ********************************************************************************************************** //

static struct Boundary_Face* constructor_empty_Boundary_Face (struct Boundary_Face* src)
{
	struct Boundary_Face* bf = malloc(sizeof *bf); // returned

	bf->bc        = src->bc;
	bf->node_nums = constructor_empty_Vector_ui(src->node_nums->extents[0]);

	return bf;
}
