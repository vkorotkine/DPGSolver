// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "mesh_periodic.h"
#include "mesh_connectivity.h"

#include <limits.h>

#include "constants_core.h"
#include "constants_bc.h"
#include "Macros.h"

#include "multiarray.h"
#include "matrix.h"
#include "vector.h"
#include "mesh.h"

// Static function declarations ************************************************************************************* //

/// \brief Container for a list of periodic face information.
struct Periodic_Face_Info {
	ptrdiff_t n_pf;                 ///< The number of periodic faces.
	struct Periodic_Face** p_faces; ///< The list of \ref Periodic_Face entities.
};

/** \brief Constructor for \ref Periodic_Face_Info.
 *	\return Standard. */
static struct Periodic_Face_Info* constructor_Periodic_Face_Info
	(const ptrdiff_t n_pf ///< \ref Periodic_Face_Info::n_pf.
	);

/// \brief Destructor for \ref Periodic_Face_Info.
static void destructor_Periodic_Face_Info
	(struct Periodic_Face_Info* pf_info ///< Standard.
	);

/// \brief Container for periodic face information.
struct Periodic_Face {
	char   dir;                 ///< The direction of the periodicity. Options: 'x', 'y', 'z'.
	double centr[DMAX-1];       ///< The centroid coordinates in directions other than the periodic one.
	struct Vector_i* node_nums; ///< The node numbers of the face vertices.
};

/** \brief Constructor for \ref Periodic_Face.
 *	\return Standard. */
struct Periodic_Face* constructor_Periodic_Face ();

/// \brief Destructor for \ref Periodic_Face.
void destructor_Periodic_Face
	(struct Periodic_Face* pf ///< Standard.
	);

/** \brief Count the number of periodic faces.
 *	\return See brief. */
static ptrdiff_t count_periodic_faces
	(const ptrdiff_t ind_pfe,                    ///< Index of the first physical face element in the mesh element list.
	 const ptrdiff_t n_pfe,                      ///< The number of physical face elements.
	 const struct const_Matrix_i*const elem_tags ///< \ref Mesh_Data::elem_tags.
	);

/** \brief Set the master and slave periodic face info for all entries in the list.
 *
 *	\note The final lists are sorted in order of the [`dir`,`centr`] pairs.
 */
static void set_pf_info
	(struct Periodic_Face_Info* pf_info[N_MS], ///< The \ref Periodic_Face_Info for the master [0] and slave [1] faces.
	 const ptrdiff_t ind_pfe,                  ///< Index of the first physical face element in the mesh element list.
	 const ptrdiff_t n_pfe,                    ///< The number of physical face elements.
	 const struct Mesh_Data*const mesh_data    ///< The \ref Mesh_Data.
	);

/** \brief Substitute the node numbers of the master faces for those of the slave faces.
 *
 *	\warning Both `f_ve` and `pf_info` must be sorted before entering this function.
 */
static void substitute_m_for_s
	(struct Multiarray_Vector_i*const f_ve,   ///< \ref Conn_info::f_ve.
	 struct Periodic_Face_Info* pf_info[N_MS] ///< The \ref Periodic_Face_Info for the master [0] and slave [1] faces.
	);

// Interface functions ********************************************************************************************** //

void correct_f_ve_for_periodic (const struct Mesh_Data*const mesh_data, struct Conn_info*const conn_info)
{
	if (mesh_data->periodic_corr == NULL)
		return;

	const ptrdiff_t d       = conn_info->d,
	                ind_pfe = get_first_volume_index(conn_info->elem_per_dim,d-1),
	                n_pfe   = conn_info->elem_per_dim->data[d-1];

	const ptrdiff_t n_pf = count_periodic_faces(ind_pfe,n_pfe,mesh_data->elem_tags);

	struct Periodic_Face_Info* pf_info[N_MS];

	for (ptrdiff_t i = 0; i < N_MS; i++)
		pf_info[i] = constructor_Periodic_Face_Info(n_pf);

	set_pf_info(pf_info,ind_pfe,n_pfe,mesh_data);

	struct Multiarray_Vector_i*const f_ve = conn_info->f_ve;
	substitute_m_for_s(f_ve,pf_info);

	for (ptrdiff_t i = 0; i < N_MS; i++)
		destructor_Periodic_Face_Info(pf_info[i]);

	// re-sort f_ve.
	struct Vector_i* ind_f_ve_final = sort_Multiarray_Vector_i(f_ve,true); // destructed
	reorder_Vector_i(conn_info->ind_f_ve,ind_f_ve_final->data);

	destructor_Vector_i(ind_f_ve_final);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Check if the physical face element is a master or slave periodic boundary.
 *	\return See brief. */
static bool check_pfe_periodic
	(const char sm, ///< Indicator for whether it is desired to check for 'M'aster or 'S'lave.
	 const int bc   ///< The value of the boundary condition.
	);

/// \brief Set the periodic face direction.
static void set_pf_dir
	(struct Periodic_Face*const pf, ///< The \ref Periodic_Face.
	 const int bc                   ///< The boundary condition.
	);

/// \brief Set the periodic face centroid.
static void set_pf_centr
	(struct Periodic_Face*const pf,               ///< The \ref Periodic_Face.
	 const struct const_Vector_i*const node_nums, ///< The entry of \ref Mesh_Data::node_nums for the current face.
	 const struct const_Matrix_d*const nodes      ///< The \ref Mesh_Data::nodes.
	);

/// \brief Sort the master and slave periodic face lists based on the [`dir`,`centr`] pairs.
static void sort_pf_info
	(struct Periodic_Face_Info* pf_info[N_MS] ///< The \ref Periodic_Face_Info for the master [0] and slave [1] faces.
	);

static struct Periodic_Face_Info* constructor_Periodic_Face_Info (const ptrdiff_t n_pf)
{
	struct Periodic_Face_Info* pf_info = malloc(sizeof *pf_info); // returned

	pf_info->n_pf = n_pf;

	struct Periodic_Face** p_faces = malloc(n_pf * sizeof *p_faces); // keep
	for (ptrdiff_t i = 0; i < n_pf; ++i)
		p_faces[i] = constructor_Periodic_Face();

	pf_info->p_faces = p_faces; // destructed

	return pf_info;
}

static void destructor_Periodic_Face_Info (struct Periodic_Face_Info* pf_info)
{
	const ptrdiff_t n_pf = pf_info->n_pf;

	for (ptrdiff_t i = 0; i < n_pf; ++i)
		destructor_Periodic_Face(pf_info->p_faces[i]);
	free(pf_info->p_faces);

	free(pf_info);
}

struct Periodic_Face* constructor_Periodic_Face ()
{
	struct Periodic_Face* pf = malloc(sizeof *pf); // returned;

	// `node_nums` must be subsequently copy contructed.
	pf->node_nums = NULL;

	return pf;
}

void destructor_Periodic_Face (struct Periodic_Face* pf)
{
	if (pf->node_nums == NULL)
		EXIT_ERROR("A Periodic Face was not fully constructed");

	destructor_Vector_i(pf->node_nums);
	free(pf);
}

static ptrdiff_t count_periodic_faces
	(const ptrdiff_t ind_pfe, const ptrdiff_t n_pfe, const struct const_Matrix_i*const elem_tags)
{
	ptrdiff_t count = 0;

	const ptrdiff_t n_max = ind_pfe+n_pfe;
	for (ptrdiff_t n = ind_pfe; n < n_max; ++n) {
		if (check_pfe_periodic('M',get_val_const_Matrix_i(n,0,elem_tags)))
			++count;
	}
	return count;
}

static void set_pf_info
	(struct Periodic_Face_Info* pf_info[N_MS], const ptrdiff_t ind_pfe, const ptrdiff_t n_pfe,
	 const struct Mesh_Data*const mesh_data)
{
	const struct const_Matrix_d*const            nodes     = mesh_data->nodes;
	const struct const_Matrix_i*const            elem_tags = mesh_data->elem_tags;
	const struct const_Multiarray_Vector_i*const node_nums = mesh_data->node_nums;

	ptrdiff_t count_ms[N_MS] = {0};

	const ptrdiff_t n_max = ind_pfe+n_pfe;
	for (ptrdiff_t n = ind_pfe; n < n_max; ++n) {
		const int bc = get_val_const_Matrix_i(n,0,elem_tags);

		int ind_ms = -1;
		if (check_pfe_periodic('M',bc))
			ind_ms = 0;
		else if (check_pfe_periodic('S',bc))
			ind_ms = 1;
		else
			continue;

		struct Periodic_Face*const pf = pf_info[ind_ms]->p_faces[count_ms[ind_ms]];
		set_pf_dir(pf,bc);
		set_pf_centr(pf,node_nums->data[n],nodes);
		set_f_node_nums(&pf->node_nums,node_nums->data[n]);

		++count_ms[ind_ms];
	}

	sort_pf_info(pf_info);

	for (ptrdiff_t i = 0; i < N_MS; i++) {
		if (count_ms[i] != pf_info[i]->n_pf)
			EXIT_ERROR("Did not find the correct number of periodic entities");
	}
}

static void substitute_m_for_s (struct Multiarray_Vector_i*const f_ve, struct Periodic_Face_Info* pf_info[N_MS])
{
	const ptrdiff_t n_pf = pf_info[0]->n_pf;

	struct Periodic_Face_Info* pf_info_m = pf_info[0],
	                         * pf_info_s = pf_info[1];

	// Form the list of pointers to the slave node_nums.
	struct Vector_i** f_ve_slave = malloc(n_pf * sizeof *f_ve_slave); // free

	const ptrdiff_t nitems = compute_size(f_ve->order,f_ve->extents);
	for (ptrdiff_t i = 0; i < n_pf; ++i) {
		struct Vector_i*const*const slave =
			bsearch(&pf_info_s->p_faces[i]->node_nums,f_ve->data,nitems,sizeof(f_ve->data[0]),cmp_Vector_i);

		if (slave == NULL)
			EXIT_ERROR("Did not find the slave periodic face in the list");

		f_ve_slave[i] = *slave;
	}

	// Substitute into f_ve.
	// Note: this must not be done until all entries are found so that the list remains sorted.
	for (ptrdiff_t i = 0; i < n_pf; ++i)
		copy_data_Vector_i_Vector_i(pf_info_m->p_faces[i]->node_nums,f_ve_slave[i]);

	free(f_ve_slave);
}

// Level 1 ********************************************************************************************************** //

/** \brief Comparison function for std::qsort between \ref Periodic_Face\* `a` and `b`.
 *	\return The comparison of `a` and `b` based first on the direction and then on the centroid.
 */
static int cmp_Periodic_Face
	(const void *a, ///< Variable 1.
	 const void *b  ///< Variable 2.
	);

static bool check_pfe_periodic (const char sm, const int bc)
{
	const int bc_base = bc % BC_STEP_SC;
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

static void set_pf_dir (struct Periodic_Face*const pf, const int bc)
{
	const int bc_base = bc % BC_STEP_SC;
	switch (bc_base) {
		case PERIODIC_XL: case PERIODIC_XR: pf->dir = 'x'; break;
		case PERIODIC_YL: case PERIODIC_YR: pf->dir = 'y'; break;
		case PERIODIC_ZL: case PERIODIC_ZR: pf->dir = 'z'; break;
		default:
			EXIT_UNSUPPORTED;
			break;
	}
}

static void set_pf_centr
	(struct Periodic_Face*const pf, const struct const_Vector_i*const node_nums,
	 const struct const_Matrix_d*const nodes)
{
	ptrdiff_t skip = DMAX;
	switch (pf->dir) {
		case 'x': skip = 0; break;
		case 'y': skip = 1; break;
		case 'z': skip = 2; break;
		default:
			EXIT_UNSUPPORTED;
			break;
	}

	double*const centr = pf->centr;
	for (int i = 0; i < DMAX-1; ++i)
		centr[i] = 0.0;

	const ptrdiff_t d     = nodes->extents[1],
	             n_max = node_nums->extents[0];
	for (ptrdiff_t n = 0; n < n_max; ++n) {
		const double*const nodes_r = get_row_const_Matrix_d(node_nums->data[n],nodes);

		ptrdiff_t ind_c = 0;
		for (ptrdiff_t dim = 0; dim < d; ++dim) {
			if (dim != skip) {
				centr[ind_c] += nodes_r[dim];
				++ind_c;
			}
		}
	}

	for (int i = 0; i < DMAX-1; ++i)
		centr[i] /= n_max;
}

static void sort_pf_info (struct Periodic_Face_Info* pf_info[N_MS])
{
	for (int i = 0; i < N_MS; ++i)
		qsort(pf_info[i]->p_faces,pf_info[i]->n_pf,sizeof(pf_info[i]->p_faces[0]),cmp_Periodic_Face);
}

// Level 2 ********************************************************************************************************** //

static int cmp_Periodic_Face (const void *a, const void *b)
{
	const struct Periodic_Face*const*const ia = (const struct Periodic_Face*const*const) a,
	                          *const*const ib = (const struct Periodic_Face*const*const) b;

	const char dir_a = (*ia)->dir,
	           dir_b = (*ib)->dir;

	if (dir_a > dir_b)
		return 1;
	else if (dir_a < dir_b)
		return -1;

	const double*const centr_a = (*ia)->centr,
	            *const centr_b = (*ib)->centr;

	for (int i = 0; i < DMAX-1; i++) {
		if (centr_a[i] > centr_b[i])
			return 1;
		else if (centr_a[i] < centr_b[i])
			return -1;
	}
	return 0;
}
