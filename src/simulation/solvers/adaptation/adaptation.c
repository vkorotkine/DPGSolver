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

#include "adaptation.h"

#include <assert.h>
#include <stdlib.h>
#include "gsl/gsl_math.h"

#include "macros.h"
#include "definitions_adaptation.h"
#include "definitions_bc.h"
#include "definitions_elements.h"
#include "definitions_h_ref.h"
#include "definitions_intrusive.h"
#include "definitions_mesh.h"

#include "face_solver_adaptive.h"
#include "element_adaptation.h"
#include "element_solver.h"
#include "volume_solver_adaptive.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "computational_elements.h"
#include "const_cast.h"
#include "geometry.h"
#include "intrusive.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "solve.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Mark volumes for adaptation based on the input strategy.
static void mark_volumes_to_adapt
	(const struct Simulation* sim, ///< Defined for \ref adapt_hp.
	 const int adapt_strategy      ///< Defined for \ref adapt_hp.
	);

/// \brief Perform hp adaptation of \ref Simulation::volumes.
static void adapt_hp_volumes
	(struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Perform hp adaptation of \ref Simulation::faces.
static void adapt_hp_faces
	(struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destruct any unused computational elements which have spawned an h-adapted parent or children.
static void destruct_unused_computational_elements
	(struct Simulation* sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void adapt_hp (struct Simulation* sim, const int adapt_strategy)
{
	assert(sim->volumes->name == IL_VOLUME_SOLVER);
	assert(sim->faces->name   == IL_FACE_SOLVER);
	assert(list_is_derived_from("solver",'e',sim));

	constructor_derived_computational_elements(sim,IL_SOLVER_ADAPTIVE); // destructed

	mark_volumes_to_adapt(sim,adapt_strategy);
	adapt_hp_volumes(sim);
	adapt_hp_faces(sim);
	destruct_unused_computational_elements(sim);
	// Don't forget to update ind_dof.

	destructor_derived_computational_elements(sim,IL_SOLVER);

	update_ind_dof(sim);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Update hp adaptation related members of \ref Adaptive_Solver_Volume\*s.
static void update_hp_members_volumes
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Update hp adaptation related members of \ref Adaptive_Solver_Face\*s.
static void update_hp_members_faces
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Update \ref Simulation::volumes for each volume where a pointer to either a child **or** parent is present.
static void update_list_volumes
	(struct Simulation*const sim ///< \ref Simulation.
	);

/// \brief Update \ref Simulation::faces for each face where a pointer to either a child **or** parent is present.
static void update_list_existing_faces
	(struct Simulation*const sim ///< \ref Simulation.
	);

/// \brief Update \ref Simulation::faces adding new faces in any h-refined volumes.
static void update_list_new_faces
	(struct Simulation*const sim ///< \ref Simulation.
	);

/// \brief Update geometry for \ref Solver_Volume_T\*s.
static void update_geometry_volumes
	(struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Update geometry for \ref Solver_Face_T\*s.
static void update_geometry_faces
	(struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Project the solution for \ref Solver_Volume_T\*s.
static void project_solution_volumes
	(struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Project the solution for \ref Solver_Face_T\*s.
static void project_solution_faces
	(struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Update the pointers of \ref Volume::faces for any volumes which are adjacent to an h-adapted face.
static void update_volume_face_pointers
	(struct Simulation* sim ///< \ref Simulation.
	);

static void mark_volumes_to_adapt (const struct Simulation* sim, const int adapt_strategy)
{
	switch (adapt_strategy) {
	case ADAPT_S_P_REFINE:
		for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
			struct Adaptive_Solver_Volume* a_s_vol = (struct Adaptive_Solver_Volume*) curr;
			a_s_vol->adapt_type = ADAPT_P_REFINE;
		}
		break;
	case ADAPT_S_P_COARSE:
		for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
			struct Adaptive_Solver_Volume* a_s_vol = (struct Adaptive_Solver_Volume*) curr;
			a_s_vol->adapt_type = ADAPT_P_COARSE;
		}
		break;
	case ADAPT_S_H_REFINE:
		for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
			struct Adaptive_Solver_Volume* a_s_vol = (struct Adaptive_Solver_Volume*) curr;
			a_s_vol->adapt_type = ADAPT_H_REFINE;
		}
		break;
	case ADAPT_S_H_COARSE:
		for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
			struct Adaptive_Solver_Volume* a_s_vol = (struct Adaptive_Solver_Volume*) curr;
			a_s_vol->adapt_type = ADAPT_H_COARSE;
		}
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",adapt_strategy);
		break;
	}
}

static void adapt_hp_volumes (struct Simulation* sim)
{
	update_hp_members_volumes(sim);
	update_list_volumes(sim);
	update_geometry_volumes(sim);
	project_solution_volumes(sim);
}

static void adapt_hp_faces (struct Simulation* sim)
{
	update_hp_members_faces(sim);
	update_list_existing_faces(sim);
	project_solution_faces(sim);
	update_volume_face_pointers(sim);
	update_list_new_faces(sim);
	update_geometry_faces(sim);
}

static void destruct_unused_computational_elements (struct Simulation* sim)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Adaptive_Solver_Volume* a_s_vol = (struct Adaptive_Solver_Volume*) curr;

		switch (a_s_vol->adapt_type) {
		case ADAPT_NONE:     // fallthrough
		case ADAPT_P_REFINE: // fallthrough
		case ADAPT_P_COARSE: // fallthrough
			break; // Do nothing.
		case ADAPT_H_COARSE:
			EXIT_ADD_SUPPORT; // Remove children
			break; // Do nothing.
		case ADAPT_H_REFINE:
			EXIT_ADD_SUPPORT; // Remove parent
			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",a_s_vol->adapt_type);
			break;
		}
	}


	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Adaptive_Solver_Face* a_s_face = (struct Adaptive_Solver_Face*) curr;

		switch (a_s_face->adapt_type) {
		case ADAPT_NONE:     // fallthrough
		case ADAPT_P_REFINE: // fallthrough
		case ADAPT_P_COARSE: // fallthrough
			break; // Do nothing.
		case ADAPT_H_COARSE:
			EXIT_ADD_SUPPORT; // Remove children
			break; // Do nothing.
		case ADAPT_H_REFINE:
			EXIT_ADD_SUPPORT; // Remove parent
			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",a_s_face->adapt_type);
			break;
		}
	}
}

// Level 1 ********************************************************************************************************** //

/// \brief Update \ref Solver_Volume_T::p_ref.
static void update_volume_p_ref
	(struct Adaptive_Solver_Volume* a_s_vol, ///< \ref Adaptive_Solver_Volume.
	 const struct Simulation* sim            ///< \ref Simulation.
	);

/** \brief Update h-refinement related members of the current \ref Adaptive_Solver_Volume.
 *
 *  \note New volumes are constructed here (either the children or the parent).
 */
static void update_volume_h
	(struct Adaptive_Solver_Volume*const a_s_vol, ///< Pointer to current \ref Adaptive_Solver_Volume.
	 const struct Simulation*const sim            ///< \ref Simulation.
	);

/** \brief Update h-refinement related members of the current \ref Adaptive_Solver_Face.
 *
 *  \note New faces are constructed here (either the children or the parent).
 */
static void update_face_h
	(struct Adaptive_Solver_Face*const a_s_face, ///< Pointer to current \ref Adaptive_Solver_Face.
	 const struct Simulation*const sim           ///< \ref Simulation.
	);

/// \brief Update \ref Solver_Face_T::p_ref.
static void update_face_p_ref
	(struct Adaptive_Solver_Face* a_s_face, ///< \ref Adaptive_Solver_Face.
	 const struct Simulation* sim           ///< \ref Simulation.
	);

/** \brief Get the index of the dominant volume.
 *  \return See brief. */
static int get_dominant_volume_index
	(const struct Adaptive_Solver_Face* a_s_face ///< \ref Adaptive_Solver_Face.
	);

/// \brief Swap \ref Face::neigh_info members.
static void swap_face_neigh
	(struct Adaptive_Solver_Face* a_s_face ///< \ref Adaptive_Solver_Face.
	);

/** \brief Return whether one of the volumes neighbouring the current face was updated.
 *  \return See brief. */
static bool neighbours_were_updated
	(const struct Adaptive_Solver_Face* a_s_face ///< \ref Adaptive_Solver_Face.
	);

/// \brief Set \ref Adaptive_Solver_Face::adapt_type.
static void set_face_adapt_type
	(struct Adaptive_Solver_Face* a_s_face ///< \ref Adaptive_Solver_Face.
	);

/** \brief If the dominant volume on the current mesh is changed as a result of the h-adaptation, update required
 *         members of the current face before adapting. */
static void swap_dominant_volume_if_necessary
	(struct Adaptive_Solver_Face* a_s_face ///< \ref Adaptive_Solver_Face.
	);

/// \brief Project all \ref Solver_Volume_T solution-related data to the p-adaptive volume.
static void compute_projection_p_volume
	(struct Adaptive_Solver_Volume* a_s_vol, ///< \ref Adaptive_Solver_Volume.
	 const struct Simulation* sim            ///< \ref Simulation.
	);

/// \brief Project all \ref Solver_Volume_T solution-related data to the h-adaptive volume.
static void compute_projection_h_volume
	(struct Adaptive_Solver_Volume* a_s_vol, ///< \ref Adaptive_Solver_Volume.
	 const struct Simulation* sim            ///< \ref Simulation.
	);

/// \brief Project all \ref Solver_Face_T solution-related data to the p-adaptive face.
static void compute_projection_p_face
	(struct Adaptive_Solver_Face* a_s_face, ///< \ref Adaptive_Solver_Face.
	 const struct Simulation* sim           ///< \ref Simulation.
	);

/// \brief Project all \ref Solver_Face_T solution-related data to the h-refined face.
static void compute_projection_h_refine_face
	(struct Adaptive_Solver_Face* a_s_face, ///< \ref Adaptive_Solver_Face.
	 const struct Simulation* sim           ///< \ref Simulation.
	);

/// \brief Insert newly created faces (internal faces of h-refined volumes) into the face list.
static void insert_new_faces
	(const struct Adaptive_Solver_Volume*const a_s_vol, ///< The h-refined volume.
	 const struct Simulation*const sim                  ///< \ref Simulation.
	);

/** \brief Return whether the face has just been newly created from h-refinement or h-coarsening.
 *  \return See brief. */
static bool is_face_h_adapted
	(const struct Adaptive_Solver_Face*const a_s_face ///< \ref Adaptive_Solver_Face.
	);

/// \brief Set \ref Volume::faces from the data for for the current side of the input \ref Face.
static void set_volume_face_ptr
	(const struct Face*const face, ///< \ref Face.
	 const int side_index          ///< The index of the side of the face under consideration.
	);

static void update_hp_members_volumes (const struct Simulation* sim)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Adaptive_Solver_Volume* a_s_vol = (struct Adaptive_Solver_Volume*) curr;

		const int adapt_type = a_s_vol->adapt_type;
		if (adapt_type == ADAPT_NONE)
			continue;

		a_s_vol->updated = true;
		switch (adapt_type) {
		case ADAPT_P_REFINE: // fallthrough
		case ADAPT_P_COARSE:
			update_volume_p_ref(a_s_vol,sim);
			break;
		case ADAPT_H_REFINE: // fallthrough
		case ADAPT_H_COARSE:
			update_volume_h(a_s_vol,sim);
			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",adapt_type);
			break;
		}
	}
}

static void update_hp_members_faces (const struct Simulation* sim)
{
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Adaptive_Solver_Face* a_s_face = (struct Adaptive_Solver_Face*) curr;

		if (!neighbours_were_updated(a_s_face))
			continue;

		set_face_adapt_type(a_s_face);
		const int adapt_type = a_s_face->adapt_type;
		if (adapt_type == ADAPT_NONE)
			continue;

		swap_dominant_volume_if_necessary(a_s_face);
		switch (adapt_type) {
		case ADAPT_P_REFINE: // fallthrough
		case ADAPT_P_COARSE:
			update_face_p_ref(a_s_face,sim);
			break;
		case ADAPT_H_REFINE: // fallthrough
		case ADAPT_H_COARSE:
			update_face_h(a_s_face,sim);
			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",adapt_type);
			break;
		}

		const int ind_dom_vol = get_dominant_volume_index(a_s_face);
		if (ind_dom_vol == 1)
			swap_face_neigh(a_s_face);
	}
}

static void update_list_volumes (struct Simulation*const sim)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Adaptive_Solver_Volume* a_s_vol = (struct Adaptive_Solver_Volume*) curr;

		if (a_s_vol->child_0) {
			assert(a_s_vol->parent == NULL);
			curr = replace_IL(sim->volumes,1,curr,a_s_vol->child_0);
		} else if (a_s_vol->parent) {
			assert(a_s_vol->child_0 == NULL);
			EXIT_ADD_SUPPORT; // curr = replace_IL(sim->volumes,n_children,curr,a_s_vol->parent);
		}
	}
}

static void update_list_existing_faces (struct Simulation*const sim)
{
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Adaptive_Solver_Face* a_s_face = (struct Adaptive_Solver_Face*) curr;

		if (a_s_face->child_0) {
			assert(a_s_face->parent == NULL);
			curr = replace_IL(sim->faces,1,curr,a_s_face->child_0);
		} else if (a_s_face->parent) {
			assert(a_s_face->child_0 == NULL);
			EXIT_ADD_SUPPORT; // curr = replace_IL(sim->faces,n_children,curr,a_s_face->parent);
		}
	}
}

static void update_list_new_faces (struct Simulation*const sim)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Adaptive_Solver_Volume* a_s_vol = (struct Adaptive_Solver_Volume*) curr;

		switch (a_s_vol->adapt_type) {
		case ADAPT_NONE:     // fallthrough
		case ADAPT_P_REFINE: // fallthrough
		case ADAPT_P_COARSE: // fallthrough
		case ADAPT_H_COARSE:
			break; // Do nothing.
		case ADAPT_H_REFINE:
			insert_new_faces(a_s_vol,sim);
			EXIT_ADD_SUPPORT;
			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",a_s_vol->adapt_type);
			break;
		}
	}
}

static void update_geometry_volumes (struct Simulation* sim)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Adaptive_Solver_Volume* a_s_vol = (struct Adaptive_Solver_Volume*) curr;
		assert(a_s_vol->adapt_type != ADAPT_NONE);

		struct Solver_Volume* s_vol = (struct Solver_Volume*) curr;
		compute_geometry_volume(s_vol,sim);
	}
}

static void update_geometry_faces (struct Simulation* sim)
{
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Adaptive_Solver_Face* a_s_face = (struct Adaptive_Solver_Face*) curr;
		struct Solver_Face* s_face            = (struct Solver_Face*) curr;

		const int adapt_type = a_s_face->adapt_type;
		switch (adapt_type) {
		case ADAPT_NONE:
			break; // Do nothing.
		case ADAPT_P_REFINE: // fallthrough
		case ADAPT_P_COARSE:
			compute_geometry_face(s_face,sim);
			break;
		case ADAPT_H_REFINE:
			EXIT_ADD_SUPPORT; // Loop over child faces and set geometry.
		default:
			EXIT_ERROR("Unsupported: %d\n",adapt_type);
			break;
		}
	}
}

static void project_solution_volumes (struct Simulation* sim)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Adaptive_Solver_Volume* a_s_vol = (struct Adaptive_Solver_Volume*) curr;

		const int adapt_type = a_s_vol->adapt_type;
		switch (adapt_type) {
		case ADAPT_NONE:
			break; // Do nothing.
		case ADAPT_P_REFINE: // fallthrough
		case ADAPT_P_COARSE:
			compute_projection_p_volume(a_s_vol,sim);
			break;
		case ADAPT_H_REFINE: // fallthrough
		case ADAPT_H_COARSE:
			compute_projection_h_volume(a_s_vol,sim);
			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",adapt_type);
			break;
		}

	}
}

static void project_solution_faces (struct Simulation* sim)
{
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Adaptive_Solver_Face* a_s_face = (struct Adaptive_Solver_Face*) curr;

		const int adapt_type = a_s_face->adapt_type;
		switch (adapt_type) {
		case ADAPT_NONE:
			break; // Do nothing.
		case ADAPT_P_REFINE: // fallthrough
		case ADAPT_P_COARSE:
			compute_projection_p_face(a_s_face,sim);
			break;
		case ADAPT_H_REFINE: // fallthrough
			compute_projection_h_refine_face(a_s_face,sim);
			break;
		case ADAPT_H_COARSE:
			EXIT_ADD_SUPPORT;
			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",adapt_type);
			break;
		}
	}
}

static void update_volume_face_pointers (struct Simulation* sim)
{
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		const struct Adaptive_Solver_Face*const a_s_face = (struct Adaptive_Solver_Face*) curr;
		if (!is_face_h_adapted(a_s_face))
			continue;

		const struct Face*const face = (struct Face*) curr;
		set_volume_face_ptr(face,0);
		if (!face->boundary)
			set_volume_face_ptr(face,0);
	}
}

// Level 2 ********************************************************************************************************** //

/// \brief Constructs the children of the current volume.
static void constructor_volumes_h_refine
	(struct Adaptive_Solver_Volume*const a_s_vol_p, ///< Current \ref Adaptive_Solver_Volume which will be a parent.
	 const struct Simulation*const sim              ///< \ref Simulation.
	);

/// \brief Constructs the children of the current face.
static void constructor_faces_h_refine
	(struct Adaptive_Solver_Face*const a_s_face_p, ///< Current \ref Adaptive_Solver_Face which will be a parent.
	 const struct Simulation*const sim              ///< \ref Simulation.
	);

/// \brief Constructs the parent of the current faces to be coarsened.
static void constructor_faces_h_coarse
	(struct Adaptive_Solver_Face*const a_s_face_c0, ///< Current \ref Adaptive_Solver_Face which will be child_0.
	 const struct Simulation*const sim              ///< \ref Simulation.
	);

/** \brief Return the maximum mesh level of any volume associated with the input (i.e. including its children if
 *         present).
 *  \return See brief. */
int find_maximum_mesh_level
	(const struct Adaptive_Solver_Volume*const a_s_vol ///< The volume.
	);

/** \brief Get the pointer to the appropriate \ref Adaptation_Element::nc_ff \ref const_Vector_T\*.
 *  \return See brief. */
const struct const_Vector_i* get_operator__nc_ff
	(const int side_index_dest,       ///< The side index of the destination.
	 const struct Solver_Face* s_face ///< \ref Solver_Face_T.
	);

/** \brief Return the number of children associated with the isotropically h-refined volume of input type.
 *  \return See brief. */
static int get_n_children
	(const struct const_Element*const element ///< The \ref Element.
	);

/** \brief Return the "ind"ex of the 's'ub face of an h-"ref"ined face.
 *  \return See brief. */
static int compute_ind_sref
	(const struct Face*const face, ///< \ref Face.
	 const int side_index          ///< The index of the side of the face under consideration.
	);

/** \brief Constructor for a new internal \ref Adaptive_Solver_Face for an h-refined \ref Adaptive_Solver_Volume.
 *  \return See brief. */
static struct Adaptive_Solver_Face* constructor_Adaptive_Solver_Face_i_new
	(const int ind_vh,                                    ///< The 'ind'ex of the child 'h'-refined 'v'olume.
	 const int ind_lf,                                    ///< Local face index.
	 const struct Adaptive_Solver_Volume*const a_s_vol_p, ///< The parent volume.
	 const struct Simulation*const sim                    ///< \ref Simulation.
	);

static void update_volume_p_ref (struct Adaptive_Solver_Volume* a_s_vol, const struct Simulation* sim)
{
	const struct Solver_Volume* s_vol = (struct Solver_Volume*) a_s_vol;
	const int p_ref = s_vol->p_ref;

	const int adapt_type = a_s_vol->adapt_type;
	switch (adapt_type) {
	case ADAPT_P_REFINE:
		assert(p_ref < sim->p_ref[1]);
		const_cast_i(&s_vol->p_ref,p_ref+1);
		break;
	case ADAPT_P_COARSE:
		assert(p_ref > sim->p_ref[0]);
		const_cast_i(&s_vol->p_ref,p_ref-1);
		break;
	case ADAPT_H_REFINE: // fallthrough
	case ADAPT_H_COARSE: // fallthrough
	default:
		EXIT_ERROR("Unsupported: %d\n",adapt_type);
		break;
	}

	a_s_vol->child_0 = NULL;
	a_s_vol->parent  = NULL;
}

static void update_volume_h (struct Adaptive_Solver_Volume*const a_s_vol, const struct Simulation*const sim)
{
	const struct Solver_Volume* s_vol = (struct Solver_Volume*) a_s_vol;
	const int ml = s_vol->ml;

	const int adapt_type = a_s_vol->adapt_type;
	switch (adapt_type) {
	case ADAPT_H_REFINE:
		assert(ml < sim->ml[1]);
		constructor_volumes_h_refine(a_s_vol,sim);
		break;
	case ADAPT_H_COARSE:
		assert(ml > sim->ml[0]);
		EXIT_ADD_SUPPORT;
		break;
	case ADAPT_P_REFINE: // fallthrough
	case ADAPT_P_COARSE: // fallthrough
	default:
		EXIT_ERROR("Unsupported: %d\n",adapt_type);
		break;
	}
}

static void update_face_h (struct Adaptive_Solver_Face*const a_s_face, const struct Simulation*const sim)
{
	const struct Solver_Face* s_face = (struct Solver_Face*) a_s_face;
	const int ml = s_face->ml;

	const int adapt_type = a_s_face->adapt_type;
	switch (adapt_type) {
	case ADAPT_H_REFINE:
		assert(ml < sim->ml[1]);
		constructor_faces_h_refine(a_s_face,sim);
		break;
	case ADAPT_H_COARSE:
		assert(ml > sim->ml[0]);
		constructor_faces_h_coarse(a_s_face,sim);
		break;
	case ADAPT_P_REFINE: // fallthrough
	case ADAPT_P_COARSE: // fallthrough
	default:
		EXIT_ERROR("Unsupported: %d\n",adapt_type);
		break;
	}
}

static void update_face_p_ref (struct Adaptive_Solver_Face* a_s_face, const struct Simulation* sim)
{
	const struct Face* face = (struct Face*) a_s_face;
	const struct Solver_Face* s_face = (struct Solver_Face*) a_s_face;

	const struct Solver_Volume* s_vol_l = (struct Solver_Volume*) face->neigh_info[0].volume;
	if (face->boundary) {
		const_cast_i(&s_face->p_ref,s_vol_l->p_ref);
	} else {
		const struct Solver_Volume* s_vol_r = (struct Solver_Volume*) face->neigh_info[1].volume;
		const_cast_i(&s_face->p_ref,GSL_MAX(s_vol_l->p_ref,s_vol_r->p_ref));
	}
	assert(s_face->p_ref >= sim->p_ref[0]);
	assert(s_face->p_ref <= sim->p_ref[1]);
}

static int get_dominant_volume_index (const struct Adaptive_Solver_Face* a_s_face)
{
	const struct Face* face = (struct Face*) a_s_face;

	if (face->boundary) {
		return 0;
	} else {
		const struct Solver_Volume* s_vol[2] = { (struct Solver_Volume*) face->neigh_info[0].volume,
		                                         (struct Solver_Volume*) face->neigh_info[1].volume, };
		const int max_ml[2] = { find_maximum_mesh_level((struct Adaptive_Solver_Volume*)s_vol[0]),
		                        find_maximum_mesh_level((struct Adaptive_Solver_Volume*)s_vol[1]), };

		if ( (max_ml[1] > max_ml[0]) ||
		    ((max_ml[1] == max_ml[0]) && (s_vol[1]->p_ref > s_vol[0]->p_ref)))
			return 1;
		else
			return 0;
	}
	EXIT_ERROR("Did not find the dominant volume index.");
}

static void swap_face_neigh (struct Adaptive_Solver_Face* a_s_face)
{
	UNUSED(a_s_face);
	EXIT_ADD_SUPPORT;
}

static bool neighbours_were_updated (const struct Adaptive_Solver_Face* a_s_face)
{
	const struct Face* face = (struct Face*) a_s_face;

	const struct Adaptive_Solver_Volume* a_s_vol_l = (struct Adaptive_Solver_Volume*) face->neigh_info[0].volume;

	if (a_s_vol_l->updated == true)
		return true;

	if (!face->boundary) {
		const struct Adaptive_Solver_Volume* a_s_vol_r =
			(struct Adaptive_Solver_Volume*) face->neigh_info[1].volume;
		if (a_s_vol_r->updated == true)
			return true;
	}
	return false;
}

static void set_face_adapt_type (struct Adaptive_Solver_Face* a_s_face)
{
	const struct Face* face = (struct Face*) a_s_face;

	const int ind_dom_vol = get_dominant_volume_index(a_s_face);
	const struct Adaptive_Solver_Volume* a_s_vol_d =
		(struct Adaptive_Solver_Volume*) face->neigh_info[ind_dom_vol].volume;

	if (!a_s_vol_d->updated) {
		a_s_face->adapt_type = ADAPT_NONE;
		return;
	}

	const int max_ml_d = find_maximum_mesh_level(a_s_vol_d);

	const struct Solver_Face* s_face    = (struct Solver_Face*) a_s_face;
	const struct Solver_Volume* s_vol_d = (struct Solver_Volume*) a_s_vol_d;
	if (s_vol_d->p_ref > s_face->p_ref)
		a_s_face->adapt_type = ADAPT_P_REFINE;
	else if (s_vol_d->p_ref < s_face->p_ref)
		a_s_face->adapt_type = ADAPT_P_COARSE;
	else if (max_ml_d > s_face->ml)
		a_s_face->adapt_type = ADAPT_H_REFINE;
	else if (max_ml_d < s_face->ml)
		a_s_face->adapt_type = ADAPT_H_COARSE;
	else
		EXIT_ERROR("Did not find the adaptation type.");
}

static void swap_dominant_volume_if_necessary (struct Adaptive_Solver_Face* a_s_face)
{
	if (get_dominant_volume_index(a_s_face) == 0)
		return;

	struct Face*const face = (struct Face*) a_s_face;
	swap_neigh_info(face);

	struct Solver_Face*const s_face = (struct Solver_Face*) a_s_face;
	// Geometry related parameters need not be reordered as they will not be used before being discarded.
	if (s_face->nf_coef->extents[0] > 0) {
		const struct const_Vector_i*const nc_ff = get_operator__nc_ff(1,s_face);
		permute_Multiarray_d_V(s_face->nf_coef,nc_ff,'R');
	}

	if (s_face->s_coef->extents[0] > 0) {
		EXIT_ADD_SUPPORT;
	}

	EXIT_UNSUPPORTED; // Ensure that all is working as expected.
}

static void compute_projection_p_volume (struct Adaptive_Solver_Volume* a_s_vol, const struct Simulation* sim)
{
	const struct Volume* vol    = (struct Volume*) a_s_vol;
	struct Solver_Volume* s_vol = (struct Solver_Volume*) a_s_vol;
	const struct Adaptation_Element* a_e = &((struct Solver_Element*)vol->element)->a_e;

	const int p_i = a_s_vol->p_ref_prev,
	          p_o = s_vol->p_ref;
	const struct Operator* cc0_vs_vs = get_Multiarray_Operator(a_e->cc0_vs_vs,(ptrdiff_t[]){0,0,p_o,p_i});

	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	UNUSED(sim);
	const char op_format = 'd';

	const struct Multiarray_d* s_coef_p = s_vol->sol_coef;
	struct Multiarray_d* s_coef =
		constructor_mm_NN1_Operator_Multiarray_d(cc0_vs_vs,s_coef_p,'C',op_format,s_coef_p->order,NULL); // moved

	destructor_Multiarray_d((struct Multiarray_d*)s_coef_p);
	s_vol->sol_coef = s_coef;

	assert(compute_size(s_vol->grad_coef->order,s_vol->grad_coef->extents) == 0); // Add support for 2nd order.
}

static void compute_projection_h_volume (struct Adaptive_Solver_Volume* a_s_vol, const struct Simulation* sim)
{
	struct Solver_Volume*const s_vol         = (struct Solver_Volume*) a_s_vol;
	const struct Solver_Volume*const s_vol_p = (struct Solver_Volume*) a_s_vol->parent;
	const int p = s_vol->p_ref;
	if (a_s_vol->parent) {
		const struct Volume*const vol_p = (struct Volume*) a_s_vol->parent;

		const struct Solver_Element* s_e     = (struct Solver_Element*) vol_p->element;
		const struct Adaptation_Element* a_e = &s_e->a_e;

		assert(compute_size(s_vol_p->sol_coef->order,s_vol_p->sol_coef->extents) != 0);
		const int ind_h = a_s_vol->ind_h;
		const struct Operator*const cc0_vs_vs = get_Multiarray_Operator(a_e->cc0_vs_vs,(ptrdiff_t[]){ind_h,0,p,p});

		destructor_Multiarray_d(s_vol->sol_coef);
		s_vol->sol_coef = constructor_mm_NN1_Operator_Multiarray_d(
			cc0_vs_vs,s_vol_p->sol_coef,'C','d',2,NULL); // keep

		if (compute_size(s_vol_p->grad_coef->order,s_vol_p->grad_coef->extents) != 0) {
			EXIT_ADD_SUPPORT;
		}
	} else if (a_s_vol->child_0) {
		EXIT_ADD_SUPPORT;
	} else {
		EXIT_ERROR("Should not be entering here.");
	}

	if (test_case_explicitly_enforces_conservation(sim))
		set_Multiarray_d(s_vol->l_mult,(struct const_Multiarray_d*)s_vol_p->l_mult);
}

static void compute_projection_p_face (struct Adaptive_Solver_Face* a_s_face, const struct Simulation* sim)
{
	const struct Face* face    = (struct Face*) a_s_face;
	struct Solver_Face* s_face = (struct Solver_Face*) a_s_face;
	const struct Volume* vol   = face->neigh_info[0].volume;
	const struct Adaptation_Element* a_e = &((struct Solver_Element*)vol->element)->a_e;

	const struct Multiarray_d*const nf_coef_p = s_face->nf_coef;
	if (compute_size(nf_coef_p->order,nf_coef_p->extents) > 0) {
		const int ind_e = get_face_element_index(face),
		          p_i   = a_s_face->p_ref_prev,
		          p_o   = s_face->p_ref;
		const struct Operator* cc0_ff_ff =
			get_Multiarray_Operator(a_e->cc0_ff_ff,(ptrdiff_t[]){ind_e,ind_e,0,0,p_o,p_i});

		// sim may be used to store a parameter establishing which type of operator to use for the computation.
		UNUSED(sim);
		const char op_format = 'd';

		struct Multiarray_d* nf_coef = constructor_mm_NN1_Operator_Multiarray_d(
			cc0_ff_ff,nf_coef_p,'C',op_format,nf_coef_p->order,NULL); // moved

		destructor_Multiarray_d((struct Multiarray_d*)nf_coef_p);
		s_face->nf_coef = nf_coef;
	}

	assert(((struct Test_Case*)sim->test_case_rc->tc)->has_2nd_order == false); // Add support.
}

static void compute_projection_h_refine_face (struct Adaptive_Solver_Face* a_s_face, const struct Simulation* sim)
{
	const struct Face* face    = (struct Face*) a_s_face;
	struct Solver_Face* s_face = (struct Solver_Face*) a_s_face;
	const struct Volume* vol   = face->neigh_info[0].volume;
	const struct Adaptation_Element* a_e = &((struct Solver_Element*)vol->element)->a_e;

	const struct Multiarray_d* nf_coef_p = s_face->nf_coef;
	if (compute_size(nf_coef_p->order,nf_coef_p->extents) > 0) {
		const int ind_e = get_face_element_index(face),
		          p_i   = a_s_face->p_ref_prev,
		          p_o   = s_face->p_ref,
		          ind_h = a_s_face->ind_h;
		const struct Operator* cc0_ff_ff =
			get_Multiarray_Operator(a_e->cc0_ff_ff,(ptrdiff_t[]){ind_e,ind_e,ind_h,0,p_o,p_i});

		// sim may be used to store a parameter establishing which type of operator to use for the computation.
		UNUSED(sim);
		const char op_format = 'd';

		destructor_Multiarray_d(s_face->nf_coef);
		s_face->nf_coef = constructor_mm_NN1_Operator_Multiarray_d(
			cc0_ff_ff,nf_coef_p,'C',op_format,nf_coef_p->order,NULL); // moved
	}

	assert(((struct Test_Case*)sim->test_case_rc->tc)->has_2nd_order == false); // Add support.
}

static void insert_new_faces
	(const struct Adaptive_Solver_Volume*const a_s_vol, const struct Simulation*const sim)
{
	struct Intrusive_List*const faces_new = constructor_empty_IL(IL_INVALID,NULL); // destructed

	const struct Adaptive_Solver_Volume*const a_s_vol_p = (struct Adaptive_Solver_Volume*)a_s_vol->parent;
	const struct Volume*const vol_p                     = (struct Volume*) a_s_vol_p;

	struct Intrusive_Link* curr = a_s_vol_p->child_0;

	const int n_children = get_n_children(vol_p->element);
	for (int ind_vh = 0; ind_vh < n_children; ++ind_vh) {
		const struct Volume*const vol = (struct Volume*) curr;
		const int n_f = vol->element->n_f;
		for (int ind_lf = 0; ind_lf < n_f; ++ind_lf) {
			const struct Face* face = vol->faces[ind_lf][0];
			if (face != NULL)
				continue;

			push_back_IL(faces_new,(struct Intrusive_Link*)
				constructor_Adaptive_Solver_Face_i_new(ind_vh,ind_lf,a_s_vol_p,sim));

// push_back in faces_new. New face is always ind_href = 0. Add the python documentation to the new version and refer to
// it for required information.
printf("%d %d\n",ind_vh,ind_lf);
		}
		curr = curr->next;
	}
//	insert_List_into_List(faces_new,sim->faces,curr);
	destructor_IL(faces_new,false);
	EXIT_UNSUPPORTED; UNUSED(a_s_vol); UNUSED(sim);
}

static bool is_face_h_adapted (const struct Adaptive_Solver_Face*const a_s_face)
{
	const int adapt_type = a_s_face->adapt_type;
	if (adapt_type == ADAPT_H_REFINE || adapt_type == ADAPT_H_COARSE)
		return true;
	return false;
}

static void set_volume_face_ptr (const struct Face*const face, const int side_index)
{
	assert(!(side_index == 1 && face->boundary));

	const struct Neigh_Info*const neigh_info = &face->neigh_info[side_index];

	const int ind_lf   = neigh_info->ind_lf,
	          ind_sref = compute_ind_sref(face,side_index);
	const_cast_Face(&neigh_info->volume->faces[ind_lf][ind_sref],face);
}

// Level 3 ********************************************************************************************************** //

/** \brief Constructor for an h-refined \ref Adaptive_Solver_Volume.
 *  \return See brief.
 *
 *  Members dependent upon potentially refined face computational elements are left unspecified.
 */
static struct Adaptive_Solver_Volume* constructor_Adaptive_Solver_Volume_h_ref
	(const int ind_h,                                     ///< The index of the child h-refined element.
	 const struct Adaptive_Solver_Volume*const a_s_vol_p, ///< The parent volume.
	 const struct const_Multiarray_d*const xyz_ve_p2_i,   ///< \ref Volume::xyz_ve of polynomial degree 2.
	 const struct Simulation*const sim                    ///< \ref Simulation.
	);

/** \brief Constructor for an h-refined \ref Adaptive_Solver_Face.
 *  \return See brief. */
static struct Adaptive_Solver_Face* constructor_Adaptive_Solver_Face_h_ref
	(const int ind_h,                                    ///< The index of the child h-refined element.
	 const struct Adaptive_Solver_Face*const a_s_face_p, ///< The parent face.
	 const struct Simulation*const sim                   ///< \ref Simulation.
	);

/** \brief Constructor for the \ref Face base of the new internal \ref Adaptive_Solver_Face.
 *  \return See brief. */
static void constructor_Face_i_new
	(struct Face*const face,                              ///< The new face.
	 const int ind_vh,                                    ///< The 'ind'ex of the child 'h'-refined 'v'olume.
	 const int ind_lf,                                    ///< Local face index.
	 const struct Adaptive_Solver_Volume*const a_s_vol_p, ///< The parent volume.
	 const struct Simulation*const sim                    ///< \ref Simulation.
	);

static void constructor_volumes_h_refine
	(struct Adaptive_Solver_Volume*const a_s_vol_p, const struct Simulation*const sim)
{
	const struct Volume*const vol_p = (struct Volume*) a_s_vol_p;

	const struct Solver_Element* s_e     = (struct Solver_Element*) vol_p->element;
	const struct Adaptation_Element* a_e = &s_e->a_e;
	const struct Operator*const vv0_vv_vv = get_Multiarray_Operator(a_e->vv0_vv_vv,(ptrdiff_t[]){0,0,2,1});

	const struct const_Multiarray_d*const xyz_ve_p2 =
		constructor_mm_NN1_Operator_const_Multiarray_d(vv0_vv_vv,vol_p->xyz_ve,'C','d',2,NULL); // destructed

	struct Intrusive_List* volumes_c = constructor_empty_IL(IL_VOLUME_SOLVER_ADAPTIVE,NULL); // destructed

	const int n_children = get_n_children(vol_p->element);
	for (int n = 1; n <= n_children; ++n) {

		push_back_IL(volumes_c,(struct Intrusive_Link*)
			constructor_Adaptive_Solver_Volume_h_ref(n,a_s_vol_p,xyz_ve_p2,sim));

		if (n == 1)
			a_s_vol_p->child_0 = volumes_c->first;

	}
	destructor_const_Multiarray_d(xyz_ve_p2);
	destructor_IL(volumes_c,false);
}

static void constructor_faces_h_refine
	(struct Adaptive_Solver_Face*const a_s_face_p, const struct Simulation*const sim)
{
	const struct Face*const face_p = (struct Face*) a_s_face_p;

	struct Intrusive_List* faces_c = constructor_empty_IL(IL_FACE_SOLVER_ADAPTIVE,NULL); // destructed

	const int n_children = get_n_children(face_p->element);
	for (int n = 1; n <= n_children; ++n) {

		push_back_IL(faces_c,(struct Intrusive_Link*)constructor_Adaptive_Solver_Face_h_ref(n,a_s_face_p,sim));

		if (n == 1)
			a_s_face_p->child_0 = faces_c->first;

	}
	destructor_IL(faces_c,false);
}

static void constructor_faces_h_coarse
	(struct Adaptive_Solver_Face*const a_s_face_c0, ///< Current \ref Adaptive_Solver_Face which will be child_0.
	 const struct Simulation*const sim              ///< \ref Simulation.
	)
{
	// Ensure that all of the faces to become children have the same dominant volume. This is required for L2
	// projection of any solution members of \ref Solver_Face_T.
	EXIT_ADD_SUPPORT; UNUSED(a_s_face_c0); UNUSED(sim);
}

int find_maximum_mesh_level (const struct Adaptive_Solver_Volume*const a_s_vol)
{
	if (a_s_vol->child_0) {
		assert(a_s_vol->parent == NULL);
		return ((struct Solver_Volume*)a_s_vol->child_0)->ml;
	} else if (a_s_vol->parent) {
		assert(a_s_vol->child_0 == NULL);
		return ((struct Solver_Volume*)a_s_vol->parent)->ml;
	} else {
		assert(a_s_vol->parent == NULL);
		assert(a_s_vol->child_0 == NULL);
		return ((struct Solver_Volume*)a_s_vol)->ml;
	}
}

const struct const_Vector_i* get_operator__nc_ff (const int side_index_dest, const struct Solver_Face* s_face)
{
	const struct Neigh_Info* neigh_info = &((struct Face*)s_face)->neigh_info[side_index_dest];

	struct Volume* vol = neigh_info->volume;
	const struct Solver_Element* s_e     = (struct Solver_Element*) vol->element;
	const struct Adaptation_Element* a_e = &s_e->a_e;

	const int ind_ord = neigh_info->ind_ord,
	          ind_e   = get_face_element_index((struct Face*)s_face),
	          p_f     = s_face->p_ref;

	return get_const_Multiarray_Vector_i(a_e->nc_ff,(ptrdiff_t[]){ind_ord,ind_e,ind_e,0,0,p_f,p_f});
}

static int get_n_children (const struct const_Element*const element)
{
	switch (element->type) {
	case POINT: return 1;  break;
	case LINE:  return 2;  break;
	case TRI:   return 4;  break;
	case QUAD:  return 4;  break;
	case TET:   return 8;  break;
	case HEX:   return 8;  break;
	case WEDGE: return 8;  break;
	case PYR:   return 10; break;
	default:
		EXIT_ERROR("Unsupported: %d",element->type);
		break;
	}
}

static int compute_ind_sref (const struct Face*const face, const int side_index)
{
	const struct Neigh_Info*const neigh_info = &face->neigh_info[side_index];
	const int ind_href = neigh_info->ind_href;

	int ind_sref = -1;
	switch (face->element->type) {
	case POINT:
		assert(ind_href == H_POINT1_V0);
		ind_sref = ind_href;
		break;
	case LINE:
		switch (ind_href) {
		case H_LINE1_V0:
			ind_sref = ind_href;
			break;
		case H_LINE2_V0: case H_LINE2_V1:
			ind_sref = ind_href-1;
			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",ind_href);
			break;
		}
		break;
	case TRI:
		EXIT_ADD_SUPPORT;
		break;
	case QUAD:
		EXIT_ADD_SUPPORT;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",face->element->type);
		break;
	}
	return ind_sref;
}

static struct Adaptive_Solver_Face* constructor_Adaptive_Solver_Face_i_new
	(const int ind_vh, const int ind_lf, const struct Adaptive_Solver_Volume*const a_s_vol_p,
	 const struct Simulation*const sim)
{
	struct Adaptive_Solver_Face*const a_s_face = calloc(1,sizeof *a_s_face); // returned

	constructor_Face_i_new((struct Face*)a_s_face,ind_vh,ind_lf,a_s_vol_p,sim);
EXIT_ADD_SUPPORT; // Possibly compute any Solver_Face_T nf/sol here as there is nothing to project.

	a_s_face->adapt_type = ADAPT_NONE;
	a_s_face->p_ref_prev = -1;
	a_s_face->ind_h      = 0;

	a_s_face->updated    = true;
	a_s_face->child_0    = NULL;
	a_s_face->parent     = NULL;

	return a_s_face;
}

// Level 4 ********************************************************************************************************** //

/** \brief Constructor for the \ref Volume base of the h-refined \ref Adaptive_Solver_Volume.
 *  \return See brief. */
static void constructor_Volume_h_ref
	(struct Volume*const vol,                             ///< The child volume.
	 const int ind_h,                                     ///< The index of the child h-refined element.
	 const struct Adaptive_Solver_Volume*const a_s_vol_p, ///< The parent volume.
	 const struct const_Multiarray_d*const xyz_ve_p2_i,   ///< \ref Volume::xyz_ve of polynomial degree 2.
	 const struct Simulation*const sim                    ///< \ref Simulation.
	);

/** \brief Constructor for the \ref Face base of the h-refined \ref Adaptive_Solver_Face.
 *  \return See brief. */
static void constructor_Face_h_ref
	(struct Face*const face,                             ///< The child face.
	 const int ind_h,                                    ///< The index of the child h-refined element.
	 const struct Adaptive_Solver_Face*const a_s_face_p, ///< The parent volume.
	 const struct Simulation*const sim                   ///< \ref Simulation.
	);

/** \brief Constructor for the \ref Solver_Volume_T base of the h-refined \ref Adaptive_Solver_Volume.
 *  \return See brief. */
static void constructor_Solver_Volume_h_ref
	(struct Solver_Volume*const s_vol,                    ///< The child volume.
	 const struct Adaptive_Solver_Volume*const a_s_vol_p, ///< The parent volume.
	 const struct Simulation*const sim                    ///< \ref Simulation.
	);

/** \brief Constructor for the \ref Solver_Face_T base of the h-refined \ref Adaptive_Solver_Face.
 *  \return See brief. */
static void constructor_Solver_Face_h_ref
	(struct Solver_Face*const s_face,                    ///< The child face.
	 const struct Adaptive_Solver_Face*const a_s_face_p, ///< The parent volume.
	 const struct Simulation*const sim                   ///< \ref Simulation.
	);

/** \brief Check whether the internal face has an edge on a domain boundary.
 *  \return `true` if yes; `false` otherwise. */
static bool internal_face_has_boundary_edge
	(const int ind_vh,               ///< Defined for \ref constructor_Face_i_new.
	 const int ind_lf,               ///< Defined for \ref constructor_Face_i_new.
	 const struct Volume*const vol_p ///< The parent \ref Volume.
	);

/** \brief Check whether the internal face is curved.
 *  \return `true` if yes; `false` otherwise. */
static bool internal_face_is_curved
	(const int ind_vh,                 ///< Defined for \ref constructor_Face_i_new.
	 const int ind_lf,                 ///< Defined for \ref constructor_Face_i_new.
	 const struct Volume*const vol_p,  ///< The parent \ref Volume.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

/** \brief Set the members of \ref Face::Neigh_Info for the newly created internal face.
 *
 *  \note The values of the members can be generated for each of the supported element/h-refiment combinations using the
 *        [h_refinement_info.py] script.
 *
 *  <!-- References: -->
 *  [h_refinement_info.py]: h_refinement/h_refinement_info.py
 */
static void set_Neigh_Info_Face_i_new
	(struct Face*const face,                             ///< The new face.
	 const int ind_vh,                                   ///< The 'ind'ex of the child 'h'-refined 'v'olume.
	 const int ind_lf,                                   ///< Local face index.
	 const struct Adaptive_Solver_Volume*const a_s_vol_p ///< The parent volume.
	);


static struct Adaptive_Solver_Volume* constructor_Adaptive_Solver_Volume_h_ref
	(const int ind_h, const struct Adaptive_Solver_Volume*const a_s_vol_p,
	 const struct const_Multiarray_d*const xyz_ve_p2_i, const struct Simulation*const sim)
{
	struct Adaptive_Solver_Volume*const a_s_vol = calloc(1,sizeof *a_s_vol); // returned

	constructor_Volume_h_ref((struct Volume*)a_s_vol,ind_h,a_s_vol_p,xyz_ve_p2_i,sim);
	constructor_Solver_Volume_h_ref((struct Solver_Volume*)a_s_vol,a_s_vol_p,sim);

	const struct Solver_Volume*const s_vol_p = (struct Solver_Volume*) a_s_vol_p;
	a_s_vol->adapt_type = ADAPT_H_REFINE;
	a_s_vol->p_ref_prev = s_vol_p->p_ref;
	a_s_vol->ind_h      = ind_h;

	a_s_vol->updated    = true;
	a_s_vol->child_0    = NULL;
	a_s_vol->parent     = (struct Intrusive_Link*) a_s_vol_p;

	return a_s_vol;
}

static struct Adaptive_Solver_Face* constructor_Adaptive_Solver_Face_h_ref
	(const int ind_h, const struct Adaptive_Solver_Face*const a_s_face_p, const struct Simulation*const sim)
{
	struct Adaptive_Solver_Face*const a_s_face = calloc(1,sizeof *a_s_face); // returned

	constructor_Face_h_ref((struct Face*)a_s_face,ind_h,a_s_face_p,sim);
	constructor_Solver_Face_h_ref((struct Solver_Face*)a_s_face,a_s_face_p,sim);

	const struct Solver_Face*const s_face_p = (struct Solver_Face*) a_s_face_p;
	a_s_face->adapt_type = ADAPT_H_REFINE;
	a_s_face->p_ref_prev = s_face_p->p_ref;
	a_s_face->ind_h      = ind_h;

	a_s_face->updated    = true;
	a_s_face->child_0    = NULL;
	a_s_face->parent     = (struct Intrusive_Link*) a_s_face_p;

	return a_s_face;
}

static void constructor_Face_i_new
	(struct Face*const face, const int ind_vh, const int ind_lf,
	 const struct Adaptive_Solver_Volume*const a_s_vol_p, const struct Simulation*const sim)
{
	const_cast_i(&face->index,-1);

	const struct Volume*const vol_l = (struct Volume*) advance_Link(ind_vh,a_s_vol_p->child_0);
	const_cast_const_Element(&face->element,
		get_element_by_type(sim->elements,compute_elem_type_sub_ce(vol_l->element->type,'f',ind_lf)));

	const struct Volume*const vol_p = (struct Volume*) a_s_vol_p;
	const_cast_b(&face->boundary,internal_face_has_boundary_edge(ind_vh,ind_lf,vol_p));
	const_cast_b(&face->curved,internal_face_is_curved(ind_vh,ind_lf,vol_p,sim));
	const_cast_i(&face->bc,BC_INVALID);

	set_Neigh_Info_Face_i_new(face,ind_vh,ind_lf,a_s_vol_p);
}

// Level 5 ********************************************************************************************************** //

/// \brief Container for data used to set \ref Face::Neigh_Info for new internal faces.
struct Info_I {
	int ind_vh_1,         ///< Index of the h-refined child volume for side 1.
	    ind_lf_1;         ///< \ref Face::Neigh_Info::ind_lf for side 1.
	const int ind_ord[2]; ///< \ref Face::Neigh_Info::ind_ord for side 0/1.
};

/** \brief Constructor for \ref Volume::bc_faces for the h-refined child volume.
 *  \return See brief. */
static const struct const_Vector_i* constructor_bc_faces_h_ref
	(const int ind_h,                 ///< The h-refinement volume index.
	 const struct Volume*const vol_p, ///< The parent volume.
	 const struct Volume*const vol_c  ///< The child volume.
	);

/** \brief Constructor for \ref Volume::bc_edges for the h-refined child volume.
 *  \return See brief. */
static const struct const_Vector_i* constructor_bc_edges_h_ref
	(const int ind_h,                 ///< The h-refinement volume index.
	 const struct Volume*const vol_p, ///< The parent volume.
	 const struct Volume*const vol_c  ///< The child volume.
	);

/** \brief Return whether the child volume is on a boundary.
 *  \return See brief. */
static bool is_child_on_boundary
	(const struct Volume*const vol ///< The current volume.
	);

/** \brief Return whether the child volume is curved.
 *  \return See brief. */
static bool is_child_curved
	(const struct Volume*const vol,    ///< The current volume.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

/** \brief Get the index of the child volume associated with the current h-refined face.
 *  \return See brief. */
static int get_ind_child
	(const int ind_h,               ///< The index of the h-refinement of the face.
	 const int side_index,          ///< The index of the side of the face under consideration.
	 const struct Face*const face_p ///< The parent \ref Face.
	);

/** \brief Return the value of \ref Face::Neigh_Info::ind_lf for the dominant \ref Face::neigh_info for an h-refined
 *         face.
 *  \return See brief. */
static int get_ind_lf_h_ref
	(const int side_index,          ///< The index of the side of the face under consideration.
	 const struct Face*const face_p ///< The parent \ref Face.
	);

/** \brief Return the value of \ref Face::Neigh_Info::ind_ord for the dominant \ref Face::neigh_info for an h-refined
 *         face.
 *  \return See brief. */
static int get_ind_ord_h_ref
	(const int side_index,          ///< The index of the side of the face under consideration.
	 const struct Face*const face_p ///< The parent \ref Face.
	);

/** \brief Get the compound index associated with the current face in relation to its neighbouring volume.
 *  \return See brief.
 *
 *  The compound index provides a unique index associated with each possible face of an element (including all h-refined
 *  possibilities). */
static int get_ind_compound
	(const int ind_h,               ///< The index of the h-refinement of the face.
	 const int side_index,          ///< The index of the side of the face under consideration.
	 const struct Face*const face_p ///< The parent \ref Face.
	);

static void constructor_Volume_h_ref
	(struct Volume*const vol, const int ind_h, const struct Adaptive_Solver_Volume*const a_s_vol_p,
	 const struct const_Multiarray_d*const xyz_ve_p2_i, const struct Simulation*const sim)
{
	const_cast_i(&vol->index,-1);

	const struct Volume*const vol_p            = (struct Volume*) a_s_vol_p;
	const struct Solver_Volume*const s_vol_p   = (struct Solver_Volume*) a_s_vol_p;
	const struct const_Element*const element_p = vol_p->element;
	const_cast_const_Element(&vol->element,
		get_element_by_type(sim->elements,compute_elem_type_sub_ce(element_p->type,'v',ind_h)));
	vol->bc_faces = constructor_bc_faces_h_ref(ind_h,vol_p,vol);
	vol->bc_edges = constructor_bc_edges_h_ref(ind_h,vol_p,vol);

	const_cast_b(&vol->boundary,is_child_on_boundary(vol));
	const_cast_b(&vol->curved,is_child_curved(vol,sim));

	const struct const_Multiarray_d* xyz_ve_p2 = NULL;
	switch (sim->domain_type) {
	case DOM_STRAIGHT: // fallthrough
	case DOM_PARAMETRIC:
		xyz_ve_p2 = xyz_ve_p2_i;
		break;
	case DOM_BLENDED:
		xyz_ve_p2 = constructor_xyz_blended('v',xyz_ve_p2,s_vol_p,sim); // destructed
		break;
	default:
		EXIT_ERROR("Unsupported: %d",sim->domain_type);
		break;
	}

	const struct Solver_Element* s_e     = (struct Solver_Element*) vol_p->element;
	const struct Adaptation_Element* a_e = &s_e->a_e;
	const struct Operator*const vv0_vv_vv = get_Multiarray_Operator(a_e->vv0_vv_vv,(ptrdiff_t[]){ind_h,0,1,2});

	const_constructor_move_const_Multiarray_d(&vol->xyz_ve,
		constructor_mm_NN1_Operator_const_Multiarray_d(vv0_vv_vv,xyz_ve_p2,'C','d',2,NULL)); // keep

	if (xyz_ve_p2 != xyz_ve_p2_i)
		destructor_const_Multiarray_d(xyz_ve_p2);
}

static void constructor_Face_h_ref
	(struct Face*const face, const int ind_h, const struct Adaptive_Solver_Face*const a_s_face_p,
	 const struct Simulation*const sim)
{
	const_cast_i(&face->index,-1);

	const struct Face*const face_p             = (struct Face*) a_s_face_p;
	const struct const_Element*const element_p = face_p->element;
	const_cast_const_Element(&face->element,
		get_element_by_type(sim->elements,compute_elem_type_sub_ce(element_p->type,'v',ind_h)));

	const_cast_b(&face->boundary,face_p->boundary);
	const_cast_b(&face->curved,face_p->curved);
	const_cast_i(&face->bc,face_p->bc);

	const struct Volume*const vol_p_0 = face_p->neigh_info[0].volume;
	const struct Adaptive_Solver_Volume*const a_s_vol_p_0 = (struct Adaptive_Solver_Volume*) vol_p_0;

	const int ind_child_0 = get_ind_child(ind_h,0,face_p);
	face->neigh_info[0].ind_lf   = get_ind_lf_h_ref(0,face_p);
	face->neigh_info[0].ind_href = 0;
	face->neigh_info[0].ind_sref = -1;
	face->neigh_info[0].ind_ord  = get_ind_ord_h_ref(0,face_p);
	face->neigh_info[0].volume   = (struct Volume*) advance_Link(ind_child_0,a_s_vol_p_0->child_0);

	if (!face->boundary) {
		const struct Solver_Face*const s_face_p               = (struct Solver_Face*) a_s_face_p;
		const struct Volume*const vol_p_1                     = face_p->neigh_info[1].volume;
		const struct Adaptive_Solver_Volume*const a_s_vol_p_1 = (struct Adaptive_Solver_Volume*) vol_p_1;

		const int max_ml = find_maximum_mesh_level(a_s_vol_p_1);
		if (max_ml > s_face_p->ml) { // neighbour is also being h-refined
			const int ind_child_1 = get_ind_child(ind_h,1,face_p);
			face->neigh_info[1].ind_lf   = get_ind_lf_h_ref(1,face_p);
			face->neigh_info[1].ind_href = 0;
			face->neigh_info[1].ind_sref = -1;
			face->neigh_info[1].ind_ord  = get_ind_ord_h_ref(1,face_p);
			face->neigh_info[1].volume   = (struct Volume*) advance_Link(ind_child_1,a_s_vol_p_1->child_0);
		} else { // neighbour is not being refined
			const int ind_compound_1 = get_ind_compound(ind_h,1,face_p);

			face->neigh_info[1].ind_lf   = face_p->neigh_info[1].ind_lf,
			face->neigh_info[1].ind_href = ind_compound_1 % NFREFMAX;
			face->neigh_info[1].ind_sref = -1;
			face->neigh_info[1].ind_ord  = get_ind_ord_h_ref(1,face_p);
			face->neigh_info[1].volume   = face_p->neigh_info[1].volume;
		}
	}
}

static void constructor_Solver_Volume_h_ref
	(struct Solver_Volume*const s_vol, const struct Adaptive_Solver_Volume*const a_s_vol_p,
	 const struct Simulation*const sim)
{
	constructor_derived_Solver_Volume((struct Volume*)s_vol,sim);

	const struct Solver_Volume*const s_vol_p = (struct Solver_Volume*) a_s_vol_p;
	const_cast_i(&s_vol->p_ref,s_vol_p->p_ref);
	const_cast_i(&s_vol->ml,s_vol_p->ml+1);
}

static void constructor_Solver_Face_h_ref
	(struct Solver_Face*const s_face, const struct Adaptive_Solver_Face*const a_s_face_p,
	 const struct Simulation*const sim)
{
	constructor_derived_Solver_Face((struct Face*)s_face,sim);

	const struct Solver_Face*const s_face_p = (struct Solver_Face*) a_s_face_p;
	const_cast_i(&s_face->p_ref,s_face_p->p_ref);
	const_cast_i(&s_face->ml,s_face_p->ml+1);
	const_cast_c(&s_face->cub_type,s_face_p->cub_type);
	s_face->constructor_Boundary_Value_fcl = s_face_p->constructor_Boundary_Value_fcl;
}

static bool internal_face_has_boundary_edge (const int ind_vh, const int ind_lf, const struct Volume*const vol_p)
{
	const struct const_Element*const e = vol_p->element;
	if (e->d < 3 || !vol_p->boundary)
		return false;

	EXIT_ADD_SUPPORT; UNUSED(ind_vh); UNUSED(ind_lf);
}

static bool internal_face_is_curved
	(const int ind_vh, const int ind_lf, const struct Volume*const vol_p, const struct Simulation*const sim)
{
	if (sim->domain_type == DOM_PARAMETRIC)
		return true;

	const struct const_Element*const e = vol_p->element;
	if (e->d < 3 || !vol_p->curved)
		return false;

	EXIT_ADD_SUPPORT; UNUSED(ind_vh); UNUSED(ind_lf);
}

static void set_Neigh_Info_Face_i_new
	(struct Face*const face, const int ind_vh, const int ind_lf,
	 const struct Adaptive_Solver_Volume*const a_s_vol_p)
{
	const struct Info_I* info_i = NULL;

// make external
	const struct Volume*const vol_p = (struct Volume*) a_s_vol_p;
	switch (vol_p->element->type) {
	case LINE: {
		assert(ind_vh == 0); // Should already have found internal faces for other volumes.
		assert(ind_lf == 1); // Only internal face for volume 0.

		static const struct Info_I nii_01 = { .ind_vh_1 = 1, .ind_lf_1 = 0, .ind_ord[0] = 0, .ind_ord[1] = 0, };
		info_i = &nii_01;
		break;
	} default:
		EXIT_ERROR("Unsupported: %d.",vol_p->element->type);
		break;
	}

	int side_index = 0;
	face->neigh_info[side_index].ind_lf   = ind_lf;
	face->neigh_info[side_index].ind_href = 0;
	face->neigh_info[side_index].ind_sref = -1;
	face->neigh_info[side_index].ind_ord  = info_i->ind_ord[0];
	face->neigh_info[side_index].volume   = (struct Volume*) advance_Link(ind_vh,a_s_vol_p->child_0);
	const_cast_Face(&face->neigh_info[side_index].volume->faces[face->neigh_info[side_index].ind_lf][0],face);

	side_index = 1;
	face->neigh_info[side_index].ind_lf   = info_i->ind_lf_1;
	face->neigh_info[side_index].ind_href = 0;
	face->neigh_info[side_index].ind_sref = -1;
	face->neigh_info[side_index].ind_ord  = info_i->ind_ord[0];
	face->neigh_info[side_index].volume   = (struct Volume*) advance_Link(info_i->ind_vh_1,a_s_vol_p->child_0);
	const_cast_Face(&face->neigh_info[side_index].volume->faces[face->neigh_info[side_index].ind_lf][0],face);
}

// Level 6 ********************************************************************************************************** //

static const struct const_Vector_i* constructor_bc_faces_h_ref
	(const int ind_h, const struct Volume*const vol_p, const struct Volume*const vol_c)
{
	const struct const_Element*const element = vol_c->element;
	struct Vector_i* bc_faces = constructor_empty_Vector_i(element->n_f); // returned
	set_to_value_Vector_i(bc_faces,BC_INVALID);

	switch (vol_p->element->type) {
	case LINE:
		switch (ind_h) {
		case H_LINE2_V0:
			bc_faces->data[0] = vol_p->bc_faces->data[0];
			break;
		case H_LINE2_V1:
			bc_faces->data[1] = vol_p->bc_faces->data[1];
			break;
		default:
			EXIT_ERROR("Unsupported: %d",ind_h);
			break;
		}
		break;
	default:
		EXIT_ERROR("Unsupported: %d",vol_p->element->type);
		break;
	}

	return (struct const_Vector_i*) bc_faces;
}

static const struct const_Vector_i* constructor_bc_edges_h_ref
	(const int ind_h, const struct Volume*const vol_p, const struct Volume*const vol_c)
{
	const struct const_Element*const element = vol_c->element;
	struct Vector_i* bc_edges = constructor_empty_Vector_i(element->n_e); // returned
	set_to_value_Vector_i(bc_edges,BC_INVALID);

	if (DIM != DMAX)
		return (struct const_Vector_i*) bc_edges;

EXIT_ADD_SUPPORT; UNUSED(ind_h);
	switch (vol_p->element->type) {
	default:
		EXIT_ERROR("Unsupported: %d",vol_p->element->type);
		break;
	}

	return (struct const_Vector_i*) bc_edges;
}

static bool is_child_on_boundary (const struct Volume*const vol)
{
	const struct const_Vector_i* bc_fe = vol->bc_faces;
	for (int i = 0; i < bc_fe->ext_0; ++i) {
		if (bc_fe->data[i] != BC_INVALID)
			return true;
	}

	if (DIM == DMAX) {
		bc_fe = vol->bc_edges;
		for (int i = 0; i < bc_fe->ext_0; ++i) {
			if (bc_fe->data[i] != BC_INVALID)
				return true;
		}
	}
	return false;
}

static bool is_child_curved (const struct Volume*const vol, const struct Simulation*const sim)
{
	if (sim->domain_type == DOM_PARAMETRIC)
		return true;

	const struct const_Vector_i* bc_fe = vol->bc_faces;
	for (int i = 0; i < bc_fe->ext_0; ++i) {
		if (bc_fe->data[i] >= BC_CURVED_START)
			return true;
	}

	if (DIM == DMAX) {
		bc_fe = vol->bc_edges;
		for (int i = 0; i < bc_fe->ext_0; ++i) {
			if (bc_fe->data[i] >= BC_CURVED_START)
				return true;
		}
	}
	return false;
}

static int get_ind_child (const int ind_h, const int side_index, const struct Face*const face_p)
{
	const int ind_compound = get_ind_compound(ind_h,side_index,face_p);
	const struct Volume*const vol_p = face_p->neigh_info[side_index].volume;

	int ind_child = -1;
	switch (vol_p->element->type) {
	case LINE:
		switch (ind_compound) {
		case H_POINT1_V0+1:
			ind_child = H_LINE2_V0; break;
		case NFREFMAX+H_POINT1_V0+1:
			ind_child = H_LINE2_V1; break;
		default:
			EXIT_ERROR("Unsupported: %d",ind_compound);
			break;
		}
		break;
#if DIM > 1
	case TRI:
		switch (ind_compound) {
		case 1*NFREFMAX+H_LINE2_V0: // fallthrough
		case 2*NFREFMAX+H_LINE2_V0:
			ind_child = H_TRI4_V0; break;
		case 0*NFREFMAX+H_LINE2_V0: // fallthrough
		case 2*NFREFMAX+H_LINE2_V1:
			ind_child = H_TRI4_V1; break;
		case 0*NFREFMAX+H_LINE2_V1: // fallthrough
		case 1*NFREFMAX+H_LINE2_V1:
			ind_child = H_TRI4_V2; break;
		default:
			EXIT_ERROR("Unsupported: %d",ind_compound);
			break;
		}
		break;
	case QUAD:
		EXIT_ADD_SUPPORT;
#endif
	default:
		EXIT_ERROR("Unsupported: %d",vol_p->element->type);
		break;
	}
	return ind_child-1;
}

static int get_ind_lf_h_ref (const int side_index, const struct Face*const face_p)
{
	const struct Volume*const vol_p = face_p->neigh_info[side_index].volume;
	switch (vol_p->element->type) {
	case LINE: // fallthrough
	case TRI:  // fallthrough
	case QUAD:
		return face_p->neigh_info[side_index].ind_lf;
		break;
	default:
		EXIT_ERROR("Unsupported: %d",vol_p->element->type);
		break;
	}
}

static int get_ind_ord_h_ref (const int side_index, const struct Face*const face_p)
{
	const struct Volume*const vol_p = face_p->neigh_info[side_index].volume;
	switch (vol_p->element->type) {
	case LINE: // fallthrough
	case TRI:  // fallthrough
	case QUAD:
		return face_p->neigh_info[side_index].ind_ord;
		break;
	default:
		EXIT_ERROR("Unsupported: %d",vol_p->element->type);
		break;
	}
}

static int get_ind_compound (const int ind_h, const int side_index, const struct Face*const face_p)
{
	const int ind_lf      = face_p->neigh_info[side_index].ind_lf,
	          ind_ord     = get_ind_ord_h_ref(side_index,face_p),
		  face_e_type = face_p->element->type;

	const int base = ind_lf*NFREFMAX;
	if (ind_ord == 0 || ind_ord == -1)
		return base + ind_h;

	int ind_f_ref = -1;
	switch (face_e_type) {
	case POINT:
		EXIT_ERROR("Should not have made it here. ind_ord: %d",ind_ord);
		break;
	case LINE:
		assert(ind_ord == 1);
		switch (ind_h) {
			case H_LINE2_V0: ind_f_ref = H_LINE2_V1; break;
			case H_LINE2_V1: ind_f_ref = H_LINE2_V0; break;
			default: EXIT_ERROR("Unsupported: %d",ind_h); break;
		}
		break;
	case TRI:
		EXIT_ADD_SUPPORT;
	case QUAD:
		EXIT_ADD_SUPPORT;
	default:
		EXIT_ERROR("Unsupported: %d",face_e_type);
		break;
	}
	return base + ind_f_ref;
}
