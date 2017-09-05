// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "face.h"
#include "intrusive.h"



#include "simulation.h"
#include "mesh.h"

#include <stdlib.h>

#include "macros.h"
#include "constants_mesh.h"
#include "constants_bc.h"
#include "constants_tol.h"

#include "volume.h"
#include "element.h"
#include "const_cast.h"
#include "math_functions.h"

// Static function declarations ************************************************************************************* //

/// \brief Container for \ref Face related mesh information.
struct Face_mesh_info {
	struct const_Element* element;        ///< The pointer to the \ref Element corresponding to the face.

	const struct const_Vector_i* ve_inds; ///< The indices of the vertices of the face.

	/// \brief Container for neighbouring info.
	struct Neigh_info_mi {
		int ind_lf;            ///< \ref Face::Neighbour_Info::ind_lf.
		struct Volume* volume; ///< \ref Face::Neighbour_Info::volume.
	} neigh_info[2]; ///< \ref Neigh_info_mi.
};

/// \brief Constructor for an individual \ref Face.
static struct Face* constructor_Face
	(const struct Simulation*const sim,         ///< \ref Simulation.
	 const struct Mesh*const mesh,              ///< \ref Mesh.
	 const struct Face_mesh_info*const face_mi, ///< \ref Face_mesh_info.
	 const int index                            ///< The face index.
	);

/// \brief Destructor for an individual \ref Face.
static void destructor_Face
	(struct Face* face ///< Standard.
	);

/** \brief Compute the vector of vertex indices for the given local face of the input volume.
 *	\return See brief. */
const struct const_Vector_i* compute_ve_inds_f
	(const struct Volume*const volume,            ///< The input \ref Volume.
	 const struct const_Vector_i*const ve_inds_v, ///< The vertex indices of the volume.
	 const int lf                                 ///< The local face under consideration.
	);

// Interface functions ********************************************************************************************** //

struct Intrusive_List* constructor_Face_List (struct Simulation*const sim, const struct Mesh*const mesh)
{
	struct Intrusive_List* faces = constructor_empty_IL();

	const struct const_Multiarray_Vector_i*const node_nums = mesh->mesh_data->node_nums;

	const struct const_Multiarray_Vector_i*const v_to_v  = mesh->mesh_conn->v_to_v,
	                                      *const v_to_lf = mesh->mesh_conn->v_to_lf;

	struct Volume** volume_array = malloc(sim->n_v * sizeof *volume_array); // free

	ptrdiff_t ind_v = 0;
	for (const struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		volume_array[ind_v] = (struct Volume*) curr;
		++ind_v;
	}

	ptrdiff_t n_f = 0;

	const ptrdiff_t v_max = sim->n_v;
	for (ptrdiff_t v = 0; v < v_max; ++v) {
		struct Volume*const volume = volume_array[v];

		const struct const_Vector_i*const v_to_v_V = v_to_v->data[v];

		const ptrdiff_t ind_v = v + mesh->mesh_data->ind_v;
		const int lf_max = v_to_v_V->ext_0;
		for (int lf = 0; lf < lf_max; ++lf) {
			if (volume->faces[lf][0] != NULL) // Already found this face.
				continue;

			const int v_neigh = v_to_v_V->data[lf];
			const struct const_Vector_i*const ve_inds =
				compute_ve_inds_f(volume,node_nums->data[ind_v],lf); // destructed

			const int ind_lf = v_to_lf->data[v]->data[lf];
			struct Volume*const volume_n = ( v_neigh == -1 ? NULL : volume_array[v_neigh] );

			struct Face_mesh_info face_mi =
				{ .element              = get_element_by_face(volume->element,lf),
				  .ve_inds              = ve_inds,
				  .neigh_info[0] = { .ind_lf = lf,
				                     .volume = volume, },
				  .neigh_info[1] = { .ind_lf = ind_lf,
				                     .volume = volume_n, },
				};

			push_back_IL(faces,(struct Intrusive_Link*) constructor_Face(sim,mesh,&face_mi,n_f));

			destructor_Vector_i((struct Vector_i*)ve_inds);

			struct Face* face = (struct Face*) faces->last;
			const_cast_Face(&volume->faces[lf][0],face);
			if (v_neigh != -1)
				const_cast_Face(&face_mi.neigh_info[1].volume->faces[ind_lf][0],face);

			++n_f;
		}
	}
	free(volume_array);

	sim->n_f = n_f;

	return faces;
}

void destructor_Faces (struct Intrusive_List* faces)
{
	for (const struct Intrusive_Link* curr = faces->first; curr; ) {
		struct Intrusive_Link* next = curr->next;
		destructor_Face((struct Face*) curr);
		curr = next;
	}
	destructor_IL(faces);
}

void const_cast_Face (const struct Face*const* dest, const struct Face*const src)
{
	*(struct Face**) dest = (struct Face*) src;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Check if the current face is on a boundary.
 *
 *	This function works similarly to \ref check_if_boundary_v.
 *
 *	\return True if on a boundary. */
static bool check_if_boundary_f
	(const int ind_lf,                                  ///< Current face component of \ref Mesh_Connectivity::v_to_lf.
	 const struct const_Multiarray_Vector_i*const f_ve, ///< Defined in \ref Element.
	 const struct const_Vector_i*const ve_inds,         ///< The vertex indices for the volume.
	 const struct Mesh_Vertices*const mesh_vert         ///< \ref Mesh_Vertices.
	);

/** \brief Check if the current face is curved.
 *
 *	This function works similarly to \ref check_if_curved_v.
 *
 *	\return True if curved. */
static bool check_if_curved_f
	(const int ind_lf,                                  ///< Current face component of \ref Mesh_Connectivity::v_to_lf.
	 const int domain_type,                             ///< \ref Simulation::domain_type.
	 const struct const_Multiarray_Vector_i*const f_ve, ///< Defined in \ref Element.
	 const struct const_Vector_i*const ve_inds,         ///< The vertex indices for the volume.
	 const struct Mesh_Vertices*const mesh_vert         ///< \ref Mesh_Vertices.
	);

/// \brief Set up the \ref Face::Neigh_Info::ind_ord
static void set_up_Face__Neigh_Info__ind_ord
	(struct Face* face ///< \ref Face.
	);

static struct Face* constructor_Face
	(const struct Simulation*const sim, const struct Mesh*const mesh, const struct Face_mesh_info*const face_mi,
	 const int index)
{
	const struct Mesh_Vertices*const mesh_vert = mesh->mesh_vert;

	struct Face* face = malloc(sizeof *face); // returned
	const_cast_i(&face->index,index);

	for (int i = 0; i < 2; ++i) {
		if (face_mi->neigh_info[i].volume) {
			face->neigh_info[i].ind_lf   = face_mi->neigh_info[i].ind_lf;
			face->neigh_info[i].ind_href = 0;
			face->neigh_info[i].ind_sref = 0;
			face->neigh_info[i].ind_ord  = -1; // Set to invalid value, to be subsequently determined.
			face->neigh_info[i].volume   = face_mi->neigh_info[i].volume;
		} else {
			face->neigh_info[i].ind_lf   = -1;
			face->neigh_info[i].ind_href = -1;
			face->neigh_info[i].ind_sref = -1;
			face->neigh_info[i].ind_ord  = -1;
			face->neigh_info[i].volume   = NULL;
		}
	}
	const_cast_const_Element(&face->element,face_mi->element);

	const int ind_lf = face_mi->neigh_info[1].ind_lf;

	const_cast_bool(&face->boundary,check_if_boundary_f(ind_lf,face->element->f_ve,face_mi->ve_inds,mesh_vert));
	const_cast_bool(&face->curved,
	                check_if_curved_f(ind_lf,sim->domain_type,face->element->f_ve,face_mi->ve_inds,mesh_vert));
	const_cast_i(&face->bc,( ind_lf > BC_STEP_SC ? ind_lf : -1 ));

	set_up_Face__Neigh_Info__ind_ord(face);

	return face;
}

static void destructor_Face (struct Face* face)
{
	UNUSED(face); // freed as part of the list.
}

const struct const_Vector_i* compute_ve_inds_f
	(const struct Volume*const volume, const struct const_Vector_i*const ve_inds_v, const int lf)
{
	const struct const_Vector_i*const f_ve_f = volume->element->f_ve->data[lf];

	const ptrdiff_t i_max = f_ve_f->ext_0;
	struct Vector_i*const ve_inds_f = constructor_empty_Vector_i(i_max); // moved
	for (ptrdiff_t i = 0; i < i_max; ++i)
		ve_inds_f->data[i] = ve_inds_v->data[f_ve_f->data[i]];

	const struct const_Vector_i*const dest = NULL;
	const_constructor_move_Vector_i(&dest,ve_inds_f); // returned
	return dest;
}

// Level 1 ********************************************************************************************************** //

/// \brief Set Face::Neigh_Info::ind_ord based on the xyz coordinates of the face vertices.
static void set_ind_ord
	(struct Neigh_Info neigh_info[2],            ///< \ref Face::Neigh_Info.
	 const struct const_Matrix_d*const xyz_ve[2] /**< The xyz coordinates of the vertices obtained from
	                                              *   \ref Volume::xyz_ve on either side. */
	);

static bool check_if_boundary_f
	(const int ind_lf, const struct const_Multiarray_Vector_i*const f_ve, const struct const_Vector_i*const ve_inds,
	 const struct Mesh_Vertices*const mesh_vert)
{
	if (ind_lf > BC_STEP_SC)
		return true;

	return check_ve_condition(f_ve,ve_inds,mesh_vert->ve_boundary,mesh_vert->ve_bc,false);
}

static bool check_if_curved_f
	(const int ind_lf, const int domain_type, const struct const_Multiarray_Vector_i*const f_ve,
	 const struct const_Vector_i*const ve_inds, const struct Mesh_Vertices*const mesh_vert)
{
	if (domain_type == DOM_PARAMETRIC)
		return true;

	if (ind_lf > 2*BC_STEP_SC)
		return true;

	return check_ve_condition(f_ve,ve_inds,mesh_vert->ve_boundary,mesh_vert->ve_bc,true);
}

static void set_up_Face__Neigh_Info__ind_ord (struct Face* face)
{
	if (face->boundary) {
		face->neigh_info[0].ind_ord = -1;
		face->neigh_info[1].ind_ord = -1;
		return;
	}

	const int d = (face->element->d)+1;

	const struct const_Matrix_d* xyz_ve[2] = { NULL };
	for (int i = 0; i < 2; ++i) {
		struct Neigh_Info neigh_info = face->neigh_info[i];
		const int ind_lf = neigh_info.ind_lf;

		const struct Volume*const volume = neigh_info.volume;
		const struct const_Element*const element = volume->element;

		const struct const_Vector_i*const f_ve_lf = element->f_ve->data[ind_lf];

		const struct const_Matrix_d* xyz_ve_v = volume->xyz_ve;

		xyz_ve[i] = constructor_copy_extract_const_Matrix_d(xyz_ve_v,f_ve_lf); // destructed
	}

	set_ind_ord(face->neigh_info,xyz_ve);
EXIT_UNSUPPORTED;

	for (int i = 0; i < 2; ++i)
		destructor_Matrix_d((struct Matrix_d*)xyz_ve[i]);
}

// Level 2 ********************************************************************************************************** //

/** Container holding information relating to the face vertex correspondence.
 *
 *	`matches_R_to_L` provides the index reordering of the vertices of the face as interpolated from the 'R'ight volume
 *	required to match the vertices of the face as interpolated from the 'L'eft volume. The opposite holds for
 *	`matches_L_to_R`.
 */
struct Vertex_Correspondence {
	struct Vector_i* matches_R_to_L, ///< Matching indices for the vertices of the 'R'ight as seen from the 'L'eft.
	               * matches_L_to_R, ///< Matching indices for the vertices of the 'L'eft as seen from the 'R'ight.
};

static struct Vertex_Correspondence* constructor_Vertex_Correspondence (const struct const_Matrix_d*const xyz_ve[2])

static void set_ind_ord (struct Neigh_Info neigh_info[2], const struct const_Matrix_d*const xyz_ve[2])
{
	const ptrdiff_t d = xyz_ve[0]->ext_1;
	if (d == 1) {
		neigh_info[0] = 0;
		neigh_info[1] = 0;
		return;
	}

	struct Vertex_Correspondence* vert_corr = constructor_Vertex_Correspondence(xyz_ve); // destructed

	const int n_ve = vert_corr->matches_R_to_L->ext_0;
	if (n_ve == 2) { // LINE
		const int n_possible = 2;
// Possibly use static const int below if compiler complains
		const int matches_possible[n_possible][n_ve] =
			{ {0, 1,},
			  {1, 0,}, };

		for (int i = 0; i < n_possible; ++i) {
			if (check_equal_Vector_i_i(matches_R_to_L,matches_possible[i]))
				neigh_info[0] = i;
			if (check_equal_Vector_i_i(matches_L_to_R,matches_possible[i]))
				neigh_info[0] = i;
		}
	}

	destructor_Vertex_Correspondence(vert_corr);
}

// Level 3 ********************************************************************************************************** //

/**	\brief Check if the face is periodic based on the coordinates of the centroid of the face vertex coordinates.
 *	\return 0 if not periodic; the value associated with the periodic boundary of the given direction otherwise. */
static int check_face_for_periodicity
	(const struct const_Matrix_d*const xyz_ve[2] ///< Defined in \ref set_ind_ord.
	);

static struct Vertex_Correspondence* constructor_Vertex_Correspondence (const struct const_Matrix_d*const xyz_ve[2])
{
	const int bc_periodic = check_face_for_periodicity(xyz_ve);

	struct Vertex_Correspondence* vert_corr = malloc(sizeof *vert_corr); // returned

	return vert_corr;
}

// Level 4 ********************************************************************************************************** //

static int check_face_for_periodicity (const struct const_Matrix_d*const xyz_ve[2])
{
	if (xyz_ve[0]->layout != 'R')
		EXIT_ERROR("Expected row major.");

	const ptrdiff_t n_ve = xyz_ve[0]->ext_0,
	                d    = xyz_ve[0]->ext_1;

	struct Vector_d* centroid[2] = { NULL };
	for (int i = 0; i < 2; ++i) {
		centroid[i] = constructor_sum_Vector_d_const_Matrix_d('R',xyz_ve[i]); // destructed

		for (int dim = 0; dim < d; ++dim)
			centroid[i]->data[dim] /= n_ve;
	}

	int bc = 0;
	for (int dim = 0; dim < d; ++dim) {
		if (!equal_d(centroid[0]->data[dim],centroid[1]->data[dim],EPS)) {
			switch (dim) {
				case 0: bc = PERIODIC_XL; break;
				case 1: bc = PERIODIC_YL; break;
				case 2: bc = PERIODIC_ZL; break;
			}
		}
	}

	for (int i = 0; i < 2; ++i)
		destructor_Vector_d(centroid[i]);

	return bc;
}
