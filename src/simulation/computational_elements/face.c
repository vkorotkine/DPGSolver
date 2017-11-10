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

#include "face.h"

#include <stdlib.h>
#include <math.h>

#include "macros.h"
#include "definitions_mesh.h"
#include "definitions_bc.h"
#include "definitions_tol.h"
#include "definitions_intrusive.h"

#include "multiarray.h"
#include "matrix.h"
#include "vector.h"

#include "const_cast.h"
#include "element.h"
#include "math_functions.h"
#include "mesh.h"
#include "mesh_readers.h"
#include "mesh_connectivity.h"
#include "mesh_vertices.h"
#include "nodes_correspondence.h"
#include "simulation.h"
#include "volume.h"

// Static function declarations ************************************************************************************* //

/// \brief Container for \ref Face related mesh information.
struct Face_mesh_info {
	const struct const_Element* element; ///< The pointer to the \ref Element corresponding to the face.

	const struct const_Vector_i* ve_inds; ///< The indices of the vertices of the face.

	/// \brief Container for neighbouring info.
	struct Neigh_info_mi {
		int ind_lf;            ///< \ref Face::Neigh_Info::ind_lf.
		struct Volume* volume; ///< \ref Face::Neigh_Info::volume.
	} neigh_info[2]; ///< \ref Neigh_info_mi.
};

/** \brief Constructor for an individual \ref Face.
 *  \return Standard. */
static struct Face* constructor_Face
	(const struct Simulation*const sim,         ///< \ref Simulation.
	 const struct Mesh*const mesh,              ///< \ref Mesh.
	 const struct Face_mesh_info*const face_mi, ///< \ref Face_mesh_info.
	 const int index                            ///< The face index.
	);

/// \brief Destructor for a \ref Face in a \ref Intrusive_List, excluding the memory for the link itself.
static void destructor_Face_link
	(struct Face* face ///< Standard.
	);

/** \brief Compute the vector of vertex indices for the given local face of the input volume.
 *  \return See brief. */
const struct const_Vector_i* compute_ve_inds_f
	(const struct Volume*const volume,            ///< The input \ref Volume.
	 const struct const_Vector_i*const ve_inds_v, ///< The vertex indices of the volume.
	 const int lf                                 ///< The local face under consideration.
	);

// Interface functions ********************************************************************************************** //

struct Intrusive_List* constructor_Faces (struct Simulation*const sim, const struct Mesh*const mesh)
{
	struct Intrusive_List* faces = constructor_empty_IL(IL_FACE,NULL);

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
				{ .element            = get_element_by_face(volume->element,lf),
				  .ve_inds            = ve_inds,
				  .neigh_info[0] = { .ind_lf = lf,
				                     .volume = volume, },
				  .neigh_info[1] = { .ind_lf = ind_lf,
				                     .volume = volume_n, },
				};

			push_back_IL(faces,(struct Intrusive_Link*) constructor_Face(sim,mesh,&face_mi,n_f));

			destructor_const_Vector_i(ve_inds);

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
		destructor_Face_link((struct Face*) curr);
		curr = next;
	}
	destructor_IL(faces);
}

void const_cast_Face (const struct Face*const* dest, const struct Face*const src)
{
	*(struct Face**) dest = (struct Face*) src;
}

int get_face_element_index (const struct Face*const face)
{
	return get_face_element_index_by_ind_lf(face->neigh_info[0].volume->element,face->neigh_info[0].ind_lf);
}

int compute_side_index_face (const struct Face* face, const struct Volume* vol)
{
	if (face->neigh_info[0].volume == vol) {
		return 0;
	} else {
		assert(face->neigh_info[1].volume == vol);
		return 1;
	}
	EXIT_ERROR("Did not find the volume.\n");
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Check if the current face is on a boundary.
 *
 *	This function works similarly to \ref check_if_boundary_v.
 *
 *  \return True if on a boundary. */
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
 *  \return True if curved. */
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

	struct Face* face = calloc(1,sizeof *face); // returned
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

static void destructor_Face_link (struct Face* face)
{
	// Do nothing.
	UNUSED(face);
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

	const struct const_Matrix_d* xyz_ve[2] = { NULL };
	for (int i = 0; i < 2; ++i) {
		struct Neigh_Info neigh_info = face->neigh_info[i];
		const int ind_lf = neigh_info.ind_lf;

		const struct Volume*const volume = neigh_info.volume;
		const struct const_Element*const element = volume->element;

		const struct const_Vector_i*const f_ve_lf = element->f_ve->data[ind_lf];

		const struct const_Matrix_d* xyz_ve_v = constructor_default_const_Matrix_d();
		set_const_Matrix_from_Multiarray_d(xyz_ve_v,volume->xyz_ve,(ptrdiff_t[]){});

		xyz_ve[i] = constructor_copy_extract_const_Matrix_d(xyz_ve_v,f_ve_lf); // destructed
		destructor_const_Matrix_d(xyz_ve_v);
	}

	set_ind_ord(face->neigh_info,xyz_ve);

	for (int i = 0; i < 2; ++i)
		destructor_const_Matrix_d(xyz_ve[i]);
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
	               * matches_L_to_R; ///< Matching indices for the vertices of the 'L'eft as seen from the 'R'ight.
};

/** \brief Constructor for \ref Vertex_Correspondence.
 *  \return Standard. */
static struct Vertex_Correspondence* constructor_Vertex_Correspondence
	(const struct const_Matrix_d*const xyz_ve[2] ///< Defined in \ref set_ind_ord.
	);

/// \brief Destructor for \ref Vertex_Correspondence.
static void destructor_Vertex_Correspondence
	(struct Vertex_Correspondence* vert_corr ///< \ref Vertex_Correspondence.
	);

static void set_ind_ord (struct Neigh_Info neigh_info[2], const struct const_Matrix_d*const xyz_ve[2])
{
	const ptrdiff_t d = xyz_ve[0]->ext_1;
	if (d == 1) {
		neigh_info[0].ind_ord = 0;
		neigh_info[1].ind_ord = 0;
		return;
	}

	struct Vertex_Correspondence* vert_corr = constructor_Vertex_Correspondence(xyz_ve); // destructed

	struct Vector_i* matches_R_to_L = vert_corr->matches_R_to_L,
	               * matches_L_to_R = vert_corr->matches_L_to_R;

	int n_possible = 0;
	int* matches_possible_i;

	const int n_ve = vert_corr->matches_R_to_L->ext_0;
	if (n_ve == 2) { // LINE
		n_possible = LINE_N_PERM;
		matches_possible_i = (int[]) { 0,1, 1,0, };
	} else if (n_ve == 3) { // TRI
		n_possible = TRI_N_PERM;
		matches_possible_i = (int[]) { 0,1,2, 1,2,0, 2,0,1, 0,2,1, 2,1,0, 1,0,2 };
	} else if (n_ve == 4) { // QUAD
		n_possible = QUAD_N_PERM;
		matches_possible_i = (int[]) { 0,1,2,3, 1,0,3,2, 2,3,0,1, 3,2,1,0, 0,2,1,3, 2,0,3,1, 1,3,0,2, 3,1,2,0};
	} else {
		EXIT_UNSUPPORTED;
	}

	struct Matrix_i* matches_possible =
		constructor_copy_Matrix_i_i('R',n_possible,n_ve,matches_possible_i); // destructed

	for (int i = 0; i < n_possible; ++i) {
		if (check_equal_Vector_i_i(matches_R_to_L,get_row_Matrix_i(i,matches_possible)))
			neigh_info[1].ind_ord = i;
		if (check_equal_Vector_i_i(matches_L_to_R,get_row_Matrix_i(i,matches_possible)))
			neigh_info[0].ind_ord = i;
	}

	destructor_Matrix_i(matches_possible);

	destructor_Vertex_Correspondence(vert_corr);
}

// Level 3 ********************************************************************************************************** //

/** \brief Check if the face is periodic based on the coordinates of the centroid of the face vertex coordinates.
 *  \return 0 if not periodic; the value associated with the periodic boundary of the given direction otherwise. */
static int check_face_for_periodicity
	(const struct const_Matrix_d*const xyz_ve[2] ///< Defined in \ref set_ind_ord.
	);

/** \brief Constructor for a \ref Vector_i\* holding the indices of the matches between the input xyz coordinate
 *         matrices.
 *  \return Standard.
 *
 *	The matrices must have `layout = 'R'`.
 */
struct Vector_i* constructor_matches_Vector_i_Matrix_d
	(const struct const_Matrix_d*const xyz_m, ///< The master xyz coordinates.
	 const struct const_Matrix_d*const xyz_s, ///< The slave xyz coordinates.
	 const struct Vector_i*const ind_skip     ///< The indices to skip (if ext_0 > 0).
	);

static struct Vertex_Correspondence* constructor_Vertex_Correspondence (const struct const_Matrix_d*const xyz_ve[2])
{
	struct Vertex_Correspondence* vert_corr = calloc(1,sizeof *vert_corr); // returned

	const int bc_periodic = check_face_for_periodicity(xyz_ve);

	int n_skip = 0;
	int ind_skip_i[1] = {-1};
	if (bc_periodic) {
		n_skip = 1;
		switch (bc_periodic) {
		case PERIODIC_XL: ind_skip_i[0] = 0; break;
		case PERIODIC_YL: ind_skip_i[0] = 1; break;
		case PERIODIC_ZL: ind_skip_i[0] = 2; break;
		default:
			EXIT_UNSUPPORTED;
			break;
		}
	}
	struct Vector_i* ind_skip = constructor_copy_Vector_i_i(n_skip,ind_skip_i); // destructed;

	vert_corr->matches_R_to_L = constructor_matches_Vector_i_Matrix_d(xyz_ve[0],xyz_ve[1],ind_skip), // destructed
	vert_corr->matches_L_to_R = constructor_matches_Vector_i_Matrix_d(xyz_ve[1],xyz_ve[0],ind_skip); // destructed

	destructor_Vector_i(ind_skip);

	return vert_corr;
}

static void destructor_Vertex_Correspondence (struct Vertex_Correspondence* vert_corr)
{
	destructor_Vector_i(vert_corr->matches_R_to_L);
	destructor_Vector_i(vert_corr->matches_L_to_R);

	free(vert_corr);
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
		centroid[i] = constructor_sum_Vector_d_const_Matrix_d('C',xyz_ve[i]); // destructed
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

struct Vector_i* constructor_matches_Vector_i_Matrix_d
	(const struct const_Matrix_d*const xyz_m, const struct const_Matrix_d*const xyz_s,
	 const struct Vector_i*const ind_skip)
{
	const ptrdiff_t ext_1  = xyz_m->ext_1,
	                n_skip = ind_skip->ext_0;
	if (n_skip >= ext_1)
		EXIT_ERROR("Too many entries in `ind_skip` (Unsupported: %td >= %td).",n_skip,ext_1);

	const ptrdiff_t ext_0 = xyz_m->ext_0;
	struct Vector_i* dest = constructor_empty_Vector_i(ext_0); // returned

	for (ptrdiff_t i = 0; i < ext_0; ++i) {
		dest->data[i] = -1;
		const double*const data_m = get_row_const_Matrix_d(i,xyz_m);
		for (ptrdiff_t i2 = 0; i2 < ext_0; ++i2) {
			const double*const data_s = get_row_const_Matrix_d(i2,xyz_s);

			double diff = 0.0;
			ptrdiff_t ind = 0;
			for (ptrdiff_t j = 0; j < ext_1; ++j) {
				if (ind < n_skip && j == ind_skip->data[ind]) {
					++ind;
					continue;
				}
				diff += fabs(data_m[j]-data_s[j]);
			}

			if (diff < EPS) {
				dest->data[i] = i2;
				break;
			}
		}
		if (dest->data[i] == -1) {
			print_Vector_i(ind_skip);
			print_const_Matrix_d(xyz_m);
			print_const_Matrix_d(xyz_s);
			EXIT_ERROR("Did not find the matching index from the slave xyz coordinates.");
		}
	}

	return dest;
}
