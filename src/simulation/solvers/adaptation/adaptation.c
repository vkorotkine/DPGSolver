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
#include "definitions_simulation.h"
#include "definitions_test_case.h"
#include "definitions_tol.h"

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
#include "math_functions.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "solve.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

//MSB: Set this parameter to false to try and see if adaptation can be used with the NURBS
// parametric mesh case.

/** Flag for whether \ref correct_non_conforming_geometry should be executed.
 *  **This should only be disabled for comparative testing purposes.** See non-conforming and geometry related tests for
 *  additional relevant comments. */
#define CORRECT_NON_CONFORMING_GEOMETRY true

/** \brief Return the maximum number of global adaptation calls which must be made in order to achieve the adaptation
 *         strategy.
 *  \return See brief. */
static int compute_max_n_adapt
	(const int adapt_strategy,                     ///< Defined for \ref adapt_hp.
	 const struct Adaptation_Data*const adapt_data ///< Defined for \ref adapt_hp.
	);

/** \brief Mark volumes for adaptation based on the input strategy.
 *  \return `true` if at least one volume was marked for adaptation; `false` otherwise. */
static bool mark_volumes_to_adapt
	(const struct Simulation*const sim,            ///< Defined for \ref adapt_hp.
	 const int adapt_strategy,                     ///< Defined for \ref adapt_hp.
	 const struct Adaptation_Data*const adapt_data ///< Defined for \ref adapt_hp.
	);

/// \brief Perform hp adaptation of \ref Simulation::volumes.
static void adapt_hp_volumes
	(struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Perform hp adaptation of \ref Simulation::faces.
static void adapt_hp_faces
	(struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destruct any unused computational elements which have been replaced with an h-adapted parent or children.
static void destruct_unused_computational_elements
	(struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Get the pointer to the appropriate \ref Geometry_Element::cv0_vgc_fgc operator.
 *  \return See brief. */
static const struct Operator* get_operator__cv0_vgc_fgc
	(const int side_index,                  ///< The index of the side of the face under consideration.
	 const struct Solver_Face*const s_face, ///< The current \ref Solver_Face_T.
	 const bool use_pg_face,                ///< Defined for \ref constructor_geom_fg.
	 const bool use_full_face               ///< Defined for \ref constructor_geom_fg.
	);

/** \brief Get the pointer to the appropriate \ref Geometry_Element::nc_fg \ref const_Vector_T\*.
 *  \return See brief. */
static const struct const_Vector_i* get_operator__nc_fg
	(const int side_index_dest,             ///< Defined for \ref permute_Multiarray_T_fc.
	 const struct Solver_Face*const s_face, ///< Defined for \ref permute_Multiarray_T_fc.
	 const bool use_pg_face                 ///< Defined for \ref constructor_geom_fg.
	);

// Interface functions ********************************************************************************************** //

void adapt_hp (struct Simulation* sim, const int adapt_strategy, const struct Adaptation_Data*const adapt_data)
{
	assert(sim->volumes->name == IL_VOLUME_SOLVER);
	assert(sim->faces->name   == IL_FACE_SOLVER);
	assert(list_is_derived_from("solver",'e',sim));

	constructor_derived_computational_elements(sim,IL_SOLVER_ADAPTIVE); // destructed

	// MSB: On the call to compute_max_n_adapt from the orders of convergence test, we get 1
	const int n_adapt = compute_max_n_adapt(adapt_strategy,adapt_data);
	for (int i = 0; i < n_adapt; ++i) {
		if(!mark_volumes_to_adapt(sim,adapt_strategy,adapt_data))
			break;

		adapt_hp_volumes(sim);
		adapt_hp_faces(sim);
		destruct_unused_computational_elements(sim);
	}

	destructor_derived_computational_elements(sim,IL_SOLVER);

	update_ind_dof(sim);
}

const struct const_Multiarray_d* constructor_geom_fg
	(const int side_index, const int side_index_dest, const struct Solver_Face*const s_face, const bool use_pg_face,
	 const bool use_full_face)
{
	const struct Operator*const cv0_vgc_fgc = get_operator__cv0_vgc_fgc(side_index,s_face,use_pg_face,use_full_face);

	const struct Face*const face           = (struct Face*) s_face;
	const struct Solver_Volume*const s_vol = (struct Solver_Volume*) face->neigh_info[side_index].volume;

	const struct const_Multiarray_d*const geom_fg =
		constructor_mm_NN1_Operator_const_Multiarray_d(cv0_vgc_fgc,s_vol->geom_coef,'C','d',2,NULL); // returned
	if (side_index != side_index_dest) {
		const struct const_Vector_i*const nc_fg = get_operator__nc_fg(side_index_dest,s_face,use_pg_face);
		permute_Multiarray_d_V((struct Multiarray_d*)geom_fg,nc_fg,'R');
	}
	return geom_fg;
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

/// \brief Destruct the input number of \ref Adaptive_Solver_Volume\*s including all bases.
static void destruct_fully_Adaptive_Solver_Volumes
	(const int n_v,                  ///< The 'n'umber of 'v'olumes to destruct.
	 struct Intrusive_Link* first,   ///< Pointer to the first volume to be destructed.
	 struct Intrusive_Link** child_0 ///< Pointer to the first child of the parent or NULL if not h-refined.
	);

/// \brief Destruct the input number of \ref Adaptive_Solver_Face\*s including all bases.
static void destruct_fully_Adaptive_Solver_Faces
	(const int n_f,                  ///< The 'n'umber of 'f'aces to destruct.
	 struct Intrusive_Link* first,   ///< Pointer to the first face to be destructed.
	 struct Intrusive_Link** child_0 ///< Pointer to the first child of the parent or NULL if not h-refined.
	);

/// \brief Update \ref Volume::index members.
static void update_index_volumes
	(struct Simulation*const sim ///< \ref Simulation.
	);

/// \brief Update \ref Face::index members.
static void update_index_faces
	(struct Simulation*const sim ///< \ref Simulation.
	);

/** \brief Return whether the vertex in the specified row of the input multiarray is in \ref Volume::xyz_ve.
 *  \return See brief. */
static bool volume_has_specified_xyz_ve
	(const struct Volume*const vol,                ///< \ref Volume.
	 const struct const_Multiarray_d*const xyz_ve, ///< Multiarray of vertex coordinates.
	 const int row                                 ///< The index of the row of xyz_ve to check.
	);

/** \brief Ensure that the space will be 1-irregular after adaptation is performed by increasing the the mesh level or
 *         the order in all volumes which would cause this constraint to be violated.
 *
 *  The term 1-irregular is used to denote:
 *  - h: A mesh for which the \ref Solver_Volume_T::ml of neighbouring volumes differs by at most 1;
 *  - p: A mesh for which the \ref Solver_Volume_T::p_ref of neighbouring volumes differs by at most 1.
 */
static void ensure_1_irregular
	(const struct Simulation*const sim ///< \ref Simulation.
	);

/** \brief Exit if the finite element space will not be 1-irregular in h (mesh) and p (order) after adaptation based on
 *         the \ref Adaptive_Solver_Volume::adapt_type values for each element.
 *  \return `true` if satisfied.
 *
 *  See comments for \ref ensure_1_irregular for the definition of "1-irregular".
 */
static bool space_will_be_1_irregular
	(const struct Simulation*const sim ///< \ref Simulation.
	);

static int compute_max_n_adapt (const int adapt_strategy, const struct Adaptation_Data*const adapt_data)
{
	switch (adapt_strategy) {
	case ADAPT_S_P_REFINE: // fallthrough
	case ADAPT_S_P_COARSE: // fallthrough
	case ADAPT_S_H_REFINE: // fallthrough
	case ADAPT_S_H_COARSE:
		return 1; // MSB: 1 is returned for the orders convergence test
	case ADAPT_S_XYZ_VE: {
		const struct const_Vector_i*const xyz_ve_ml = adapt_data->xyz_ve_ml,
		                           *const xyz_ve_p  = adapt_data->xyz_ve_p;
		for (int i = 0; xyz_ve_ml && i < xyz_ve_ml->ext_0; ++i)
			assert(xyz_ve_ml->data[i] <= ML_MAX);
		for (int i = 0; xyz_ve_p && i < xyz_ve_p->ext_0; ++i)
			assert(xyz_ve_p->data[i] <= P_MAX);

		return ML_MAX+P_MAX;
	} default:
		EXIT_ERROR("Unsupported: %d\n",adapt_strategy);
		break;
	}
}

static bool mark_volumes_to_adapt
	(const struct Simulation*const sim, const int adapt_strategy, const struct Adaptation_Data*const adapt_data)
{
	bool adapting = false;
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Adaptive_Solver_Volume*const a_s_vol = (struct Adaptive_Solver_Volume*) curr;
		a_s_vol->adapt_type = ADAPT_NONE;
	}

	switch (adapt_strategy) {
	case ADAPT_S_XYZ_VE: {
		const struct const_Multiarray_d*const xyz_ve_ref = adapt_data->xyz_ve_refine;
		const struct const_Vector_i*const xyz_ve_ml      = adapt_data->xyz_ve_ml,
		                           *const xyz_ve_p       = adapt_data->xyz_ve_p;
		const ptrdiff_t n_ve = xyz_ve_ref->extents[0];
		assert(!xyz_ve_ml || (n_ve == xyz_ve_ml->ext_0));
		assert(!xyz_ve_p  || (n_ve == xyz_ve_p->ext_0));

		for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
			const struct Volume*const vol               = (struct Volume*) curr;
			const struct Solver_Volume*const s_vol      = (struct Solver_Volume*) curr;
			struct Adaptive_Solver_Volume*const a_s_vol = (struct Adaptive_Solver_Volume*) curr;

			for (int i = 0; i < n_ve; ++i) {
				if (xyz_ve_ml &&
				    (s_vol->ml < xyz_ve_ml->data[i] && volume_has_specified_xyz_ve(vol,xyz_ve_ref,i))) {
					if (xyz_ve_ml->data[i] <= ML_MAX) {
						adapting = true;
						a_s_vol->adapt_type = ADAPT_H_REFINE;
					}
				} else if (xyz_ve_p &&
				           (s_vol->p_ref < xyz_ve_p->data[i] && volume_has_specified_xyz_ve(vol,xyz_ve_ref,i))) {
					if (xyz_ve_p->data[i] <= P_MAX) {
						adapting = true;
						a_s_vol->adapt_type = ADAPT_P_REFINE;
					}
				} else if (xyz_ve_p &&
				           (s_vol->p_ref > xyz_ve_p->data[i] && volume_has_specified_xyz_ve(vol,xyz_ve_ref,i))) {
					if (xyz_ve_p->data[i] >= P_MIN) {
						adapting = true;
						a_s_vol->adapt_type = ADAPT_P_COARSE;
					}
				}
			}
		}
		ensure_1_irregular(sim);
		break;
	}
	case ADAPT_S_P_REFINE:
		adapting = true;
		for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
			struct Adaptive_Solver_Volume* a_s_vol = (struct Adaptive_Solver_Volume*) curr;
			a_s_vol->adapt_type = ADAPT_P_REFINE;
		}
		break;
	case ADAPT_S_P_COARSE:
		adapting = true;
		for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
			struct Adaptive_Solver_Volume* a_s_vol = (struct Adaptive_Solver_Volume*) curr;
			a_s_vol->adapt_type = ADAPT_P_COARSE;
		}
		break;
	case ADAPT_S_H_REFINE:
		adapting = true;
		for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
			struct Adaptive_Solver_Volume* a_s_vol = (struct Adaptive_Solver_Volume*) curr;
			a_s_vol->adapt_type = ADAPT_H_REFINE;
		}
		break;
	case ADAPT_S_H_COARSE:
		adapting = true;
		for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
			struct Adaptive_Solver_Volume* a_s_vol = (struct Adaptive_Solver_Volume*) curr;
			a_s_vol->adapt_type = ADAPT_H_COARSE;
		}
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",adapt_strategy);
		break;
	}
	assert(space_will_be_1_irregular(sim));

	return adapting;
}

static void adapt_hp_volumes (struct Simulation* sim)
{
	// MSB: To_Show_Philip

	// MSB: Perform the adaptation
	// - Here, we will get the xyz_ve_p2 values and use this to find the 
	// 	xyz_ve values of the refined volumes.
	update_hp_members_volumes(sim);

	// MSB: Update the list by adding in the volume structures for the new volumes.
	// This means new volumes will be added to the linked list structure holding them
	update_list_volumes(sim);

	// MSB: Update the geometry information for the given volume (metric terms). 
	// Here is where we will potentially add in the NURBS metric computation
	// option to be able to perform refinement using the NURBS mapping approach
	update_geometry_volumes(sim);

	// MSB: Use the operator cc0_vs_vs (coefficient to coefficient of type 
	// volume solution to volume solution) to get the coefficients of the solution
	// basis functions for the updated volume (split volume)
	project_solution_volumes(sim);

	// MSB: Now that there are new volume elements in the volume list, 
	// modify the index of the volumes in the list.
	update_index_volumes(sim);
}

static void adapt_hp_faces (struct Simulation* sim)
{
	// MSB: To_Show_Philip

	// MSB: Create the child faces from the parent ones that were adapted. 
	// Create the intrusive linked list first.
	update_hp_members_faces(sim);

	// MSB: Update the intrusive list of faces by adding in face data structures
	// for the child faces from the refined elements
	update_list_existing_faces(sim);

	// MSB: Get the coefficients for the solution for the refined face by projecting
	// the solution.
	project_solution_faces(sim);

	// MSB: For the volumes in the mesh, after the adaptation, set the pointers to the
	// faces (newly created) in the adaptation process
	update_volume_face_pointers(sim);

	// MSB: Add face data structures to the face linked list
	update_list_new_faces(sim);

	// MSB: Compute the geometry/metric terms for the new faces
	update_geometry_faces(sim);

	// MSB: With new faces now in the face linked list, update the indeces of the faces
	// in the list
	update_index_faces(sim);
}

static void destruct_unused_computational_elements (struct Simulation* sim)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Adaptive_Solver_Volume* a_s_vol = (struct Adaptive_Solver_Volume*) curr;

		switch (a_s_vol->adapt_type) {
		case ADAPT_NONE:     // fallthrough
		case ADAPT_P_REFINE: // fallthrough
		case ADAPT_P_COARSE:
			break; // Do nothing.
		case ADAPT_H_COARSE:
			EXIT_ADD_SUPPORT; // Remove children (pass NULL instead of curr in destruct_fully_...)
			break; // Do nothing.
		case ADAPT_H_REFINE:
			destruct_fully_Adaptive_Solver_Volumes(1,a_s_vol->parent,&curr);
			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",a_s_vol->adapt_type);
			break;
		}
	}
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next)
		initialize_Adaptive_Solver_Volume((struct Adaptive_Solver_Volume*)curr);

	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Adaptive_Solver_Face* a_s_face = (struct Adaptive_Solver_Face*) curr;

		switch (a_s_face->adapt_type) {
		case ADAPT_NONE:     // fallthrough
		case ADAPT_P_REFINE: // fallthrough
		case ADAPT_P_COARSE: // fallthrough
		case ADAPT_H_CREATE:
		case ADAPT_GEOM:
			break; // Do nothing.
		case ADAPT_H_COARSE:
			EXIT_ADD_SUPPORT; // Remove children (pass NULL instead of curr in destruct_fully_...)
			break; // Do nothing.
		case ADAPT_H_REFINE:
			destruct_fully_Adaptive_Solver_Faces(1,a_s_face->parent,&curr);
			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",a_s_face->adapt_type);
			break;
		}
	}
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next)
		initialize_Adaptive_Solver_Face((struct Adaptive_Solver_Face*)curr);
}

static const struct Operator* get_operator__cv0_vgc_fgc
	(const int side_index, const struct Solver_Face*const s_face, const bool use_pg_face, const bool use_full_face)
{
	const struct Face*const face            = (struct Face*) s_face;
	const struct Volume*const vol           = face->neigh_info[side_index].volume;
	const struct Solver_Volume*const s_vol  = (struct Solver_Volume*) vol;
	const struct Geometry_Element*const g_e = &((struct Solver_Element*)vol->element)->g_e;

	const int ind_lf   = face->neigh_info[side_index].ind_lf,
	          ind_href = ( use_full_face ? 0 : face->neigh_info[side_index].ind_href );
	const int p_v = s_vol->p_ref,
	          p_f = ( !use_pg_face ? s_face->p_ref : compute_face_geometry_order(s_face) );

	const int curved = ( (s_face->cub_type == 's') ? 0 : 1 );
	assert(curved);

	return get_Multiarray_Operator(g_e->cv0_vgc_fgc,(ptrdiff_t[]){ind_lf,ind_href,0,p_f,p_v});
}

static const struct const_Vector_i* get_operator__nc_fg
	(const int side_index_dest, const struct Solver_Face*const s_face, const bool use_pg_face)
{
	const struct Neigh_Info*const neigh_info = &((struct Face*)s_face)->neigh_info[side_index_dest];
	const struct Volume*const vol            = neigh_info->volume;
	const struct Geometry_Element*const g_e  = &((struct Solver_Element*)vol->element)->g_e;

	const int ind_ord = neigh_info->ind_ord,
	          ind_e   = get_face_element_index((struct Face*)s_face),
	          p_f     = ( !use_pg_face ? s_face->p_ref : compute_face_geometry_order(s_face) );
	const int curved = ( (s_face->cub_type == 's') ? 0 : 1 );

	return get_const_Multiarray_Vector_i(g_e->nc_fg[curved],(ptrdiff_t[]){ind_ord,ind_e,ind_e,0,0,p_f,p_f});
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
	 struct Flux_Input*const flux_i,                    ///< \ref Flux_Input_T.
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

/// \brief Advance the \ref Intrusive_Link\* volume pointer until the next link has a different parent.
static void advance_to_end_of_parent_volume
	(struct Intrusive_Link** curr ///< Pointer to the current \ref Adaptive_Solver_Volume.
	);

/// \brief Advance the \ref Intrusive_Link\* face pointer until the next link has a different parent.
static void advance_to_end_of_parent_face
	(struct Intrusive_Link** curr ///< Pointer to the current \ref Adaptive_Solver_Face.
	);

/** \brief Update the \ref Adaptive_Solver_Volume::adapt_type parameter in all neighbouring elements which would cause
 *         the 1-irregularity constraint to be violated after adaptation.
 *  \return `true` if a neighbour was updated; `false` otherwise. */
static bool update_neighbours_for_1_irregular
	(struct Adaptive_Solver_Volume*const a_s_vol, ///< The \ref Adaptive_Solver_Volume.
	 const char hp_type                           /**< The type of parameter. Options: mes'h' level,
	                                               *   'p'olynomial order. */
	);

/** \brief Exit if the input volume is not 1-irregular with respect to its neighbours for the given hp_type.
 *  \return `true` if satisfied. */
static bool volume_will_be_1_irregular
	(const struct Adaptive_Solver_Volume*const a_s_vol, ///< The \ref Adaptive_Solver_Volume.
	 const char hp_type                                 /**< The type of parameter. Options: mes'h' level,
	                                                     *   'p'olynomial order. */
	);

/** \brief Correct non-conforming geometry.
 *
 *  In the case of non-conforming meshes (either in mesh level (h) or in order (p)), the geometry representation of
 *  internally curved faces is corrected such that the same representation is obtained when interpolating the geometry
 *  values from either of the neighbouring volumes. The face geometry nodes are set such that:
 *  - p: the lowest of the two orders of the neighbouring faces is used.
 *  - h: each of the h-refined sub-faces have the geometry corresponding to the coarse face geometry interpolated to the
 *       sub-faces.
 *  - hp: As for 'h', but using the lowest order of both neighbouring volumes.
 */
static void correct_non_conforming_geometry
	(const struct Simulation*const sim ///< \ref Simulation.
	);

static void update_hp_members_volumes (const struct Simulation* sim)
{	
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Adaptive_Solver_Volume* a_s_vol = (struct Adaptive_Solver_Volume*) curr;

		const int adapt_type = a_s_vol->adapt_type;
		if (adapt_type == ADAPT_NONE)
			continue;

		a_s_vol->updated = true;

		// MSB: Loop through the volumes in the mesh. If we need to adapt the given volume,
		// we will be at this point, where we will either p refine or h refine it
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
	struct Flux_Input*const flux_i = constructor_Flux_Input_e(sim); // destructed
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Adaptive_Solver_Volume* a_s_vol = (struct Adaptive_Solver_Volume*) curr;

		switch (a_s_vol->adapt_type) {
		case ADAPT_NONE:     // fallthrough
		case ADAPT_P_REFINE: // fallthrough
		case ADAPT_P_COARSE: // fallthrough
		case ADAPT_H_COARSE:
			break; // Do nothing.
		case ADAPT_H_REFINE:
			insert_new_faces(a_s_vol,flux_i,sim);
			advance_to_end_of_parent_volume(&curr);
			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",a_s_vol->adapt_type);
			break;
		}
	}
	destructor_Flux_Input(flux_i);
}

static void update_geometry_volumes (struct Simulation* sim)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Adaptive_Solver_Volume* a_s_vol = (struct Adaptive_Solver_Volume*) curr;
		if (a_s_vol->adapt_type == ADAPT_NONE)
			continue;

		struct Solver_Volume* s_vol = (struct Solver_Volume*) curr;

		// MSB: Test this out here. We want to compute, for the adapted volume, the 
		// metrics using the NURBS mapping (for the NURBS enhanced case)
		if(NURBS_geometry){
			compute_NURBS_geometry_volume(true,s_vol,sim);
		}else{
			compute_geometry_volume(true,s_vol,sim);
		}
	}
}

static void update_geometry_faces (struct Simulation* sim)
{
	correct_non_conforming_geometry(sim);

	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Adaptive_Solver_Face* a_s_face = (struct Adaptive_Solver_Face*) curr;
		struct Solver_Face* s_face            = (struct Solver_Face*) curr;

		const int adapt_type = a_s_face->adapt_type;
		if (adapt_type == ADAPT_NONE)
			continue;

		// MSB: Test this out here. We want to compute, for the adapted face, the 
		// metrics using the NURBS mapping (for the NURBS enhanced case)
		if(NURBS_geometry){
			compute_NURBS_geometry_face(s_face,sim);
		} else{
			compute_geometry_face(s_face,sim);
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
			set_volume_face_ptr(face,1);
	}
}

static void destruct_fully_Adaptive_Solver_Volumes
	(const int n_v, struct Intrusive_Link* first, struct Intrusive_Link** child_0)
{
	advance_to_end_of_parent_volume(child_0);

	struct Intrusive_Link* curr = first;
	for (int n = 0; n < n_v; ++n) {
		struct Intrusive_Link* curr_n = curr->next;

		struct Volume* vol = (struct Volume*)curr;
		destructor_derived_Adaptive_Solver_Volume(vol);
		destructor_derived_Solver_Volume(vol);
		destructor_Volume(vol);

		curr = curr_n;
	}
}

static void destruct_fully_Adaptive_Solver_Faces
	(const int n_f, struct Intrusive_Link* first, struct Intrusive_Link** child_0)
{
	advance_to_end_of_parent_face(child_0);

	struct Intrusive_Link* curr = first;
	for (int n = 0; n < n_f; ++n) {
		struct Intrusive_Link* curr_n = curr->next;

		struct Face* face = (struct Face*)curr;
		destructor_derived_Adaptive_Solver_Face(face);
		destructor_derived_Solver_Face(face);
		destructor_Face(face);

		curr = curr_n;
	}
}

static void update_index_volumes (struct Simulation*const sim)
{
	int index = 0;
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Volume* vol = (struct Volume*) curr;
		const_cast_i(&vol->index,index);
		++index;
	}
}

static void update_index_faces (struct Simulation*const sim)
{
	int index = 0;
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Face* face = (struct Face*) curr;
		const_cast_i(&face->index,index);
		++index;
	}
}

static bool volume_has_specified_xyz_ve
	(const struct Volume*const vol, const struct const_Multiarray_d*const xyz_ve, const int row)
{
	bool contains_v = false;

	const double*const xyz_ve_data = get_row_const_Multiarray_d(row,xyz_ve);

	const bool requires_transpose = ( vol->xyz_ve->layout != 'R' ? true : false );
	if (requires_transpose)
		transpose_Multiarray_d((struct Multiarray_d*)vol->xyz_ve,true);

	const ptrdiff_t n_ve = vol->xyz_ve->extents[0];
	for (int i = 0; i < n_ve; ++i) {
		if (norm_diff_d(DIM,xyz_ve_data,get_row_const_Multiarray_d(i,vol->xyz_ve),"Inf") < NODETOL_MESH)
			contains_v = true;
	}
	if (requires_transpose)
		transpose_Multiarray_d((struct Multiarray_d*)vol->xyz_ve,true);

	return contains_v;
}

static void ensure_1_irregular (const struct Simulation*const sim)
{
	for (bool updated = 1; updated; ) {
		updated = 0;
		for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
			struct Adaptive_Solver_Volume*const a_s_vol = (struct Adaptive_Solver_Volume*) curr;
			if (update_neighbours_for_1_irregular(a_s_vol,'h') ||
			    update_neighbours_for_1_irregular(a_s_vol,'p'))
				updated = 1;
		}
	}
}

static bool space_will_be_1_irregular (const struct Simulation*const sim)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		const struct Adaptive_Solver_Volume*const a_s_vol = (struct Adaptive_Solver_Volume*) curr;
		assert(volume_will_be_1_irregular(a_s_vol,'h'));
		assert(volume_will_be_1_irregular(a_s_vol,'p'));
	}
	return true;
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
	 struct Flux_Input*const flux_i,                      ///< \ref Flux_Input_T.
	 const struct Simulation*const sim                    ///< \ref Simulation.
	);

/** \brief Get the index for \ref OP_IND_H for the given \ref Adaptive_Solver_Face::ind_h.
 *  \return See brief.
 *
 *  \note This function is only necessary to treat the special case of POINT faces which "refine"/"coarsen" into
 *        exacatly the same POINT computational element, forcing the index to be 0.
 */
static int get_ind_h_operator
	(const struct Adaptive_Solver_Face*const a_s_face ///< \ref Adaptive_Solver_Face.
	);

/// \brief Increment the \ref Adaptive_Solver_Volume::adapt_type for the current hp_type.
static void increment_adapt_type
	(struct Adaptive_Solver_Volume*const a_s_vol, ///< The \ref Adaptive_Solver_Volume.
	 const char hp_type                           /**< The type of parameter. Options: mes'h' level,
	                                               *   'p'olynomial order. */
	);

/** \brief Return the value of the mes'h' level or 'p'olynomial order in the volume after adaptation is performed.
 *  \return See brief. */
static int compute_updated_hp
	(const struct Adaptive_Solver_Volume*const a_s_vol, ///< The \ref Adaptive_Solver_Volume.
	 const char hp_type                                 /**< The type of parameter. Options: mes'h' level,
	                                                     *   'p'olynomial order. */
	);

/** \brief Check whether the input face requires a geometry correction.
 *  \return `true` if yes; `false` otherwise. */
static bool face_requires_geom_update
	(const struct Solver_Face*const s_face, ///< The current \ref Solver_Face_T.
	 const int ml,                          ///< The current 'm'esh 'l'evel under consideration.
	 const char ml_vf                       /**< Indicates whether the 'm'esh 'l'evel of the 'v'olumes or 'f'ace
	                                         *   should be considered to judge the necessity for the update. */
	);

/** \brief Return `true` if one of the neighbours of the face was updated, `false` otherwise.
 *  \return See brief. */
static bool face_neigh_was_updated
	(const struct Adaptive_Solver_Face*const a_s_face ///< \ref Adaptive_Solver_Face.
	);

/** \brief Correct the geometry coefficients for the two volumes adjacent to the input face if required due to
 *         non-conforming geometry representations on either side. */
static void correct_non_conforming_face_geometry
	(const struct Solver_Face*const s_face, ///< The current \ref Solver_Face_T.
	 const int ml                           ///< The 'm'esh 'l'evel of volumes to be corrected.
	);

/** \brief Correct additional faces which are only indirectly affected by non-conforming face geometry corrections
 * performed in \ref correct_non_conforming_face_geometry.
 *
 * This function handles cases where vertices (and edges in 3D) of certain volumes which are on conforming geometry
 * faces must also be corrected due to the movement of vertex/edge geometry nodes during previous corrections.
 */
static void correct_additional_non_conforming_face_geometry
	(const struct Simulation*const sim ///< \ref Simulation.
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
	a_s_vol->p_ref_prev = p_ref;

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
	const struct Solver_Face* s_face = (struct Solver_Face*) a_s_face;
	const_cast_i(&s_face->p_ref,compute_face_reference_order(s_face));

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
//	else if (s_vol_d->p_ref == s_face->p_ref && max_ml_d == s_face->ml) // Should not reach this point.
	else
		EXIT_ERROR("Did not find the adaptation type (face: %d %d, vol: %d %d).",
		           s_face->p_ref,s_face->ml,s_vol_d->p_ref,s_vol_d->ml);
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

	if (s_face->s_coef && s_face->s_coef->extents[0] > 0) {
		EXIT_ADD_SUPPORT;
	}
}

static void compute_projection_p_volume (struct Adaptive_Solver_Volume* a_s_vol, const struct Simulation* sim)
{
	const struct Volume* vol    = (struct Volume*) a_s_vol;
	struct Solver_Volume* s_vol = (struct Solver_Volume*) a_s_vol;
	const struct Adaptation_Element* a_e = &((struct Solver_Element*)vol->element)->a_e;

	const int p_i = a_s_vol->p_ref_prev,
	          p_o = s_vol->p_ref;
	assert(abs(p_i-p_o) <= 1);
	const struct Operator* cc0_vs_vs = get_Multiarray_Operator(a_e->cc0_vs_vs,(ptrdiff_t[]){0,0,p_o,p_i});

	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	UNUSED(sim);
	const char op_format = 'd';

	const struct Multiarray_d* s_coef_p = s_vol->sol_coef;
	struct Multiarray_d* s_coef =
		constructor_mm_NN1_Operator_Multiarray_d(cc0_vs_vs,s_coef_p,'C',op_format,s_coef_p->order,NULL); // moved

	destructor_Multiarray_d((struct Multiarray_d*)s_coef_p);
	s_vol->sol_coef = s_coef;

	struct Multiarray_d*const g_coef_p = s_vol->grad_coef;
	switch (sim->method) {
	case METHOD_DG: {
		const struct Operator* cc0_vr_vr = get_Multiarray_Operator(a_e->cc0_vr_vr,(ptrdiff_t[]){0,0,p_o,p_i});
		const ptrdiff_t ext_0 = cc0_vr_vr->op_std->ext_0;

		assert(g_coef_p->order == 3);
		const ptrdiff_t extents[3] = { ext_0, g_coef_p->extents[1], g_coef_p->extents[2], };
		resize_Multiarray_d(g_coef_p,g_coef_p->order,extents);
		break;
	} case METHOD_DPG: {
		if (compute_size(g_coef_p->order,g_coef_p->extents) != 0)
			EXIT_ADD_SUPPORT; // Project as for the solution.
		break;
	} default:
		EXIT_ERROR("Unsupported: %d.",sim->method);
		break;
	}
}

static void compute_projection_h_volume (struct Adaptive_Solver_Volume* a_s_vol, const struct Simulation* sim)
{
	const struct Test_Case*const test_case = (struct Test_Case*) sim->test_case_rc->tc;

	struct Solver_Volume*const s_vol         = (struct Solver_Volume*) a_s_vol;
	const struct Solver_Volume*const s_vol_p = (struct Solver_Volume*) a_s_vol->parent;
	const int p = s_vol->p_ref;
	if (a_s_vol->parent) {
		const struct Volume*const vol_p = (struct Volume*) a_s_vol->parent;

		const struct Solver_Element* s_e     = (struct Solver_Element*) vol_p->element;
		const struct Adaptation_Element* a_e = &s_e->a_e;

		assert(compute_size(s_vol_p->sol_coef->order,s_vol_p->sol_coef->extents) != 0);
		const int ind_h = a_s_vol->ind_h;
		const struct Operator*const cc0_vs_vs =
			get_Multiarray_Operator(a_e->cc0_vs_vs,(ptrdiff_t[]){ind_h,0,p,p});

		destructor_Multiarray_d(s_vol->sol_coef);
		s_vol->sol_coef = constructor_mm_NN1_Operator_Multiarray_d(
			cc0_vs_vs,s_vol_p->sol_coef,'C','d',2,NULL); // keep

		if (test_case->required_unknowns[2]) {
			assert(compute_size(s_vol_p->grad_coef->order,s_vol_p->grad_coef->extents) != 0);
			const struct Operator*const cc0_vr_vr =
				get_Multiarray_Operator(a_e->cc0_vr_vr,(ptrdiff_t[]){ind_h,0,p,p});

			destructor_Multiarray_d(s_vol->grad_coef);
			s_vol->grad_coef = constructor_mm_NN1_Operator_Multiarray_d(
				cc0_vr_vr,s_vol_p->grad_coef,'C','d',3,NULL); // keep
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

	const struct Multiarray_d*const s_coef_p = s_face->s_coef;
	switch (sim->method) {
	case METHOD_DG: {
		break; // Do nothing.
	} case METHOD_DPG: {
		if (s_coef_p && compute_size(s_coef_p->order,s_coef_p->extents) != 0)
			EXIT_ADD_SUPPORT; // Project as for nf_coef.
		break;
	} default:
		EXIT_ERROR("Unsupported: %d.",sim->method);
		break;
	}
}

static void compute_projection_h_refine_face (struct Adaptive_Solver_Face* a_s_face, const struct Simulation* sim)
{
	const struct Test_Case*const test_case = (struct Test_Case*) sim->test_case_rc->tc;

	const struct Face*const face_p          = (struct Face*) a_s_face->parent;
	const struct Solver_Face*const s_face_p = (struct Solver_Face*) face_p;
	const struct Volume*const vol_p         = face_p->neigh_info[0].volume;
	const struct Adaptation_Element* a_e = &((struct Solver_Element*)vol_p->element)->a_e;

	struct Solver_Face* s_face = (struct Solver_Face*) a_s_face;

	// Note: nf_coef not present on boundaries.
	const struct Multiarray_d*const nf_coef_p = s_face_p->nf_coef;

	if (test_case->required_unknowns[1] && (compute_size(nf_coef_p->order,nf_coef_p->extents) == 0))
		assert(face_p->boundary);

	if (test_case->required_unknowns[1] &&
	    (compute_size(nf_coef_p->order,nf_coef_p->extents) != 0)) {
		const int ind_e = get_face_element_index(face_p),
		          p_i   = a_s_face->p_ref_prev,
		          p_o   = s_face->p_ref,
		          ind_h = get_ind_h_operator(a_s_face);
		const struct Operator* cc0_ff_ff =
			get_Multiarray_Operator(a_e->cc0_ff_ff,(ptrdiff_t[]){ind_e,ind_e,ind_h,0,p_o,p_i});

		// sim may be used to store a parameter establishing which type of operator to use for the computation.
		UNUSED(sim);
		const char op_format = 'd';

		destructor_Multiarray_d(s_face->nf_coef);
		s_face->nf_coef = constructor_mm_NN1_Operator_Multiarray_d(
			cc0_ff_ff,nf_coef_p,'C',op_format,nf_coef_p->order,NULL); // moved
	}

	assert(test_case->required_unknowns[3] == false); // Add support;
}

static void insert_new_faces
	(const struct Adaptive_Solver_Volume*const a_s_vol, struct Flux_Input*const flux_i,
	 const struct Simulation*const sim)
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
				constructor_Adaptive_Solver_Face_i_new(ind_vh,ind_lf,a_s_vol_p,flux_i,sim));
		}
		curr = curr->next;
	}
	push_back_IL_List(sim->faces,faces_new);
	destructor_IL(faces_new,false);
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

static void advance_to_end_of_parent_volume (struct Intrusive_Link** curr)
{
	struct Adaptive_Solver_Volume* a_s_vol = (struct Adaptive_Solver_Volume*) (*curr);
	while (*curr) {
		struct Adaptive_Solver_Volume* a_s_vol_n = (struct Adaptive_Solver_Volume*) (*curr)->next;
		if (!a_s_vol_n || (a_s_vol_n->parent != a_s_vol->parent))
			return;
		*curr = (*curr)->next;
	}
}

static void advance_to_end_of_parent_face (struct Intrusive_Link** curr)
{
	struct Adaptive_Solver_Face* a_s_face = (struct Adaptive_Solver_Face*) (*curr);
	while (*curr) {
		struct Adaptive_Solver_Face* a_s_face_n = (struct Adaptive_Solver_Face*) (*curr)->next;
		if (!a_s_face_n || (a_s_face_n->parent != a_s_face->parent))
			return;
		*curr = (*curr)->next;
	}
}

static bool update_neighbours_for_1_irregular (struct Adaptive_Solver_Volume*const a_s_vol, const char hp_type)
{
	bool updated = false;

	const struct Volume*const vol = (struct Volume*) a_s_vol;
	const int mp_up = compute_updated_hp(a_s_vol,hp_type);
	for (int i = 0; i < NFMAX;    ++i) {
	for (int j = 0; j < NSUBFMAX; ++j) {
		const struct Volume*const vol_n = get_volume_neighbour(vol,vol->faces[i][j]);
		if (!vol_n)
			continue;

		struct Adaptive_Solver_Volume*const a_s_vol_n = (struct Adaptive_Solver_Volume*) vol_n;
		const int mp_up_n = compute_updated_hp(a_s_vol_n,hp_type);
		if (mp_up-mp_up_n > 1) {
			updated = true;
			increment_adapt_type(a_s_vol_n,hp_type);
		}
	}}
	return updated;
}

static bool volume_will_be_1_irregular (const struct Adaptive_Solver_Volume*const a_s_vol, const char hp_type)
{
	const struct Volume*const vol = (struct Volume*) a_s_vol;
	const int mp_up = compute_updated_hp(a_s_vol,hp_type);
	for (int i = 0; i < NFMAX;    ++i) {
	for (int j = 0; j < NSUBFMAX; ++j) {
		const struct Volume*const vol_n = get_volume_neighbour(vol,vol->faces[i][j]);
		if (!vol_n)
			continue;

		const int mp_up_n = compute_updated_hp((struct Adaptive_Solver_Volume*)vol_n,hp_type);
		if (abs(mp_up-mp_up_n) > 1)
			EXIT_ERROR("More than 1-irregular mesh levels (vol: %d, neigh: %d).\n",vol->index,vol_n->index);
	}}
	return true;
}

static void correct_non_conforming_geometry (const struct Simulation*const sim)
{
	if (!CORRECT_NON_CONFORMING_GEOMETRY || NURBS_geometry)
		return;
	if (DIM == 1)
		return;

	/** Correct sequentially all volumes for each mesh level such that it can always be assumed that the geometry
	 *  from coarser (dominant for geometry) volumes is correct. */
	for (int ml = ML_MIN; ml < ML_MAX; ++ml) {
		for (int ml_type = 0; ml_type < 2; ++ml_type) {
			const char ml_vf = ( ml_type == 0 ? 'f' : 'v' );
			for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
				const struct Solver_Face*const s_face            = (struct Solver_Face*) curr;
				const struct Adaptive_Solver_Face*const a_s_face = (struct Adaptive_Solver_Face*) curr;

				if (!face_requires_geom_update(s_face,ml,ml_vf) || !face_neigh_was_updated(a_s_face))
					continue;

				correct_non_conforming_face_geometry(s_face,ml);
			}
			correct_additional_non_conforming_face_geometry(sim);
		}
	}

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		const struct Adaptive_Solver_Volume*const a_s_vol = (struct Adaptive_Solver_Volume*) curr;
		if (!a_s_vol->updated_geom_nc)
			continue;

		struct Solver_Volume*const s_vol = (struct Solver_Volume*) curr;
		correct_internal_xyz_blended(s_vol,sim);
		compute_geometry_volume(false,s_vol,sim);
	}

	// Mark conforming faces, requiring geometry updating if currently marked with ADAPT_NONE.
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		const struct Face*const face = (struct Face*) curr;
		const struct Adaptive_Solver_Volume*const a_s_vol[2] =
			{ (struct Adaptive_Solver_Volume*) face->neigh_info[0].volume,
                    (struct Adaptive_Solver_Volume*) face->neigh_info[1].volume, };
		struct Adaptive_Solver_Face*const a_s_face = (struct Adaptive_Solver_Face*) curr;

		if (a_s_face->adapt_type != ADAPT_NONE)
			continue;

		if (a_s_vol[0]->updated_geom_nc || (!face->boundary && a_s_vol[1]->updated_geom_nc))
			a_s_face->adapt_type = ADAPT_GEOM;
	}
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

/** \brief Constructor for the \ref Solver_Face_T base of the new internal \ref Adaptive_Solver_Face.
 *  \return See brief. */
static void constructor_Solver_Face_i_new
	(struct Solver_Face*const s_face,                     ///< Defined for \ref constructor_Face_i_new.
	 const struct Adaptive_Solver_Volume*const a_s_vol_p, ///< Defined for \ref constructor_Face_i_new.
	 struct Flux_Input*const flux_i,                      ///< \ref Flux_Input_T.
	 const struct Simulation*const sim                    ///< Defined for \ref constructor_Face_i_new.
	);

/** \brief Constructor for the 'f'ace 'g'eometry coefficients of the geometry order of the volume having input
 *         side_index.
 *  \return See brief. */
static const struct const_Multiarray_d* constructor_coef_fg_vol
	(const int side_index,                 ///< Side index of the volume whose geometry is being corrected.
	 const struct Solver_Face*const s_face ///< \ref Solver_Face_T.
	);

/** \brief Constructor for the 'f'ace 'g'eometry coefficients of the geometry order of the volume having input
 *         side_index using the opposite volume as dominant.
 *  \return See brief. */
static const struct const_Multiarray_d* constructor_coef_fg_vol_opp_dom
	(const int side_index,                 ///< Side index of the volume whose geometry is being corrected.
	 const struct Solver_Face*const s_face ///< \ref Solver_Face_T.
	);

/** \brief Constructor for the indices of the columns of the \ref Geometry_Element::cc0_vgc_fgc operator which
 *         correspond to the face coefficients.
 *  \return See brief.
 *
 *  Note that the geometry basis in the volume forms a complete basis on each of the faces corresponding exactly to the
 *  volume basis of the geometry of the face element. Each of the rows of the cc0_vgc_fgc operator thus contains a
 *  single '1' and is otherwise '0'.
 */
static const struct const_Vector_i* constructor_cc0_vgc_fgc_indices
	(const int side_index,                 ///< Side index of the volume whose geometry is being corrected.
	 const struct Solver_Face*const s_face ///< \ref Solver_Face_T.
	);

static void constructor_volumes_h_refine
	(struct Adaptive_Solver_Volume*const a_s_vol_p, const struct Simulation*const sim)
{
	const struct Volume*const vol_p = (struct Volume*) a_s_vol_p;

	const struct Solver_Element* s_e     = (struct Solver_Element*) vol_p->element;
	const struct Adaptation_Element* a_e = &s_e->a_e;
	const struct Operator*const vv0_vv_vv = get_Multiarray_Operator(a_e->vv0_vv_vv,(ptrdiff_t[]){0,0,2,1});

	const struct const_Multiarray_d*const xyz_ve_p2 =
		constructor_mm_NN1_Operator_const_Multiarray_d(vv0_vv_vv,vol_p->xyz_ve,'R','d',2,NULL); // destructed

	struct Intrusive_List* volumes_c = constructor_empty_IL(IL_VOLUME_SOLVER_ADAPTIVE,NULL); // destructed

	const int n_children = get_n_children(vol_p->element);

	// MSB: Here, we are getting the vertices of the parent volume (P1). Use these
	// to find the P2 vertices (there should be (p+1)^2 in 2D). Then, these new
	// points (equidistant for vertices) will be used to split up this orginal
	// parent volume into 4 volumes

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
	 struct Flux_Input*const flux_i, const struct Simulation*const sim)
{
	struct Adaptive_Solver_Face*const a_s_face = calloc(1,sizeof *a_s_face); // returned

	constructor_Face_i_new((struct Face*)a_s_face,ind_vh,ind_lf,a_s_vol_p,sim);
	constructor_Solver_Face_i_new((struct Solver_Face*)a_s_face,a_s_vol_p,flux_i,sim);

	a_s_face->adapt_type = ADAPT_H_CREATE;
	a_s_face->p_ref_prev = -1;
	a_s_face->ind_h      = 0;

	a_s_face->updated    = true;
	a_s_face->child_0    = NULL;
	a_s_face->parent     = NULL;

	return a_s_face;
}

static int get_ind_h_operator (const struct Adaptive_Solver_Face*const a_s_face)
{
	const struct Face*const face = (struct Face*) a_s_face;
	if (face->element->type == POINT) {
		assert(a_s_face->ind_h == 1);
		return a_s_face->ind_h-1;
	}
	return a_s_face->ind_h;
}


static void increment_adapt_type (struct Adaptive_Solver_Volume*const a_s_vol, const char hp_type)
{
	const int adapt_type = a_s_vol->adapt_type;
	switch (hp_type) {
	case 'h':
		switch (adapt_type) {
			case ADAPT_H_COARSE: a_s_vol->adapt_type = ADAPT_NONE;             break;
			case ADAPT_NONE:     a_s_vol->adapt_type = ADAPT_H_REFINE;         break;
			case ADAPT_H_REFINE: EXIT_ERROR("Should not be entering here.\n"); break;
			case ADAPT_P_COARSE: a_s_vol->adapt_type = ADAPT_HP_RC;            break;
			case ADAPT_P_REFINE: a_s_vol->adapt_type = ADAPT_HP_RR;            break;
			default:             EXIT_ERROR("Unsupported: %d\n",adapt_type);   break;
		}
		break;
	case 'p':
		switch (adapt_type) {
			case ADAPT_P_COARSE: a_s_vol->adapt_type = ADAPT_NONE;             break;
			case ADAPT_NONE:     a_s_vol->adapt_type = ADAPT_P_REFINE;         break;
			case ADAPT_P_REFINE: EXIT_ERROR("Should not be entering here.\n"); break;
			case ADAPT_H_COARSE: a_s_vol->adapt_type = ADAPT_HP_CR;            break;
			case ADAPT_H_REFINE: a_s_vol->adapt_type = ADAPT_HP_RR;            break;
			default:             EXIT_ERROR("Unsupported: %d\n",adapt_type);   break;
		}
		break;
	default:
		EXIT_ERROR("Unsupported: %c\n",hp_type);
		break;
	}
}

static int compute_updated_hp (const struct Adaptive_Solver_Volume*const a_s_vol, const char hp_type)
{
	const int adapt_type = a_s_vol->adapt_type;
	const struct Solver_Volume*const s_vol = (struct Solver_Volume*) a_s_vol;
	switch (hp_type) {
	case 'h': {
		const int ml = s_vol->ml;
		switch (adapt_type) {
		case ADAPT_H_COARSE:
			return ml-1;
			break;
		case ADAPT_H_REFINE:
			return ml+1;
			break;
		case ADAPT_P_COARSE: // fallthrough
		case ADAPT_P_REFINE: // fallthrough
		case ADAPT_NONE:
			return ml;
		default:
			EXIT_ERROR("Unsupported: %d\n",adapt_type);
			break;
		}
		break;
	} case 'p': {
		const int p = s_vol->p_ref;
		switch (adapt_type) {
		case ADAPT_P_COARSE:
			return p-1;
			break;
		case ADAPT_P_REFINE:
			return p+1;
			break;
		case ADAPT_H_COARSE: // fallthrough
		case ADAPT_H_REFINE: // fallthrough
		case ADAPT_NONE:
			return p;
		default:
			EXIT_ERROR("Unsupported: %d\n",adapt_type);
			break;
		}
		break;
	} default:
		EXIT_ERROR("Unsupported: %c\n",hp_type);
		break;
	}
}

static bool face_requires_geom_update (const struct Solver_Face*const s_face, const int ml, const char ml_vf)
{
	const struct Face*const face = (struct Face*) s_face;
	switch (ml_vf) {
	case 'v': {
		const struct Solver_Volume*const s_vol[2] = { (struct Solver_Volume*) face->neigh_info[0].volume,
		                                              (struct Solver_Volume*) face->neigh_info[1].volume, };

		if (face_is_conforming(s_face) || !(s_vol[0]->ml == ml || s_vol[1]->ml == ml) || !face->curved)
			return false;
		break;
	} case 'f':
		if (face_is_conforming(s_face) || !(s_face->ml == ml) || !face->curved)
			return false;
		break;
	default:
		EXIT_ERROR("Unsupported: %c\n",ml_vf);
		break;
	}
	return true;
}

static bool face_neigh_was_updated (const struct Adaptive_Solver_Face*const a_s_face)
{
	const struct Face*const face                         = (struct Face*) a_s_face;
	for (int i = 0; i < 2; ++i) {
		const struct Adaptive_Solver_Volume*const a_s_vol =
			(struct Adaptive_Solver_Volume*) face->neigh_info[i].volume;
		if (a_s_vol->updated)
			return true;
	}
	return false;
}

static void correct_non_conforming_face_geometry (const struct Solver_Face*const s_face, const int ml)
{
	const struct Face*const face = (struct Face*) s_face;
	const struct Solver_Volume*const s_vol[2] = { (struct Solver_Volume*) face->neigh_info[0].volume,
	                                              (struct Solver_Volume*) face->neigh_info[1].volume, };

	for (int si = 0; si < 2; ++si) {
		if (s_vol[si]->ml != ml)
			continue;

		const int si_n = (si+1) % 2;
		if (!((s_vol[si_n]->ml < s_vol[si]->ml) || (s_vol[si_n]->p_ref < s_vol[si]->p_ref)))
			continue;

		struct Adaptive_Solver_Volume*const a_s_vol = (struct Adaptive_Solver_Volume*) s_vol[si];
		a_s_vol->updated_geom_nc = true;

		const struct const_Multiarray_d*const geom_coef_fg = constructor_coef_fg_vol(si,s_face); // destructed
		const struct const_Vector_i*const coef_inds = constructor_cc0_vgc_fgc_indices(si,s_face); // destructed

		update_rows_Multiarray_d((struct Multiarray_d*)s_vol[si]->geom_coef,geom_coef_fg,coef_inds);

		destructor_const_Multiarray_d(geom_coef_fg);
		destructor_const_Vector_i(coef_inds);
	}
}

static void correct_additional_non_conforming_face_geometry (const struct Simulation*const sim)
{
	// The loop until no corrections were made is required for the treatment of cases where several layers of
	// volumes all share a vertex/edge on a non-conforming face.
	for (bool updated = true; updated; ) {
		updated = false;
		for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
			const struct Face*const face          = (struct Face*) curr;
			const struct Solver_Face*const s_face = (struct Solver_Face*) curr;
			if (face->boundary || !face->curved || !face_is_conforming(s_face))
				continue;

			struct Adaptive_Solver_Volume*const a_s_vol[2] =
				{ (struct Adaptive_Solver_Volume*) face->neigh_info[0].volume,
	                          (struct Adaptive_Solver_Volume*) face->neigh_info[1].volume, };

			if (a_s_vol[0]->updated_geom_nc == a_s_vol[1]->updated_geom_nc)
				continue;

			const int si_dom   = ( a_s_vol[0]->updated_geom_nc ? 0 : 1 ),
			          si_neigh = (si_dom+1) % 2;
			updated = true;

			// Do not mark 'updated_geom_nc' as `true` such that multiple faces may be updated if needed.
			a_s_vol[si_neigh]->updated_geom_nc_face = true;

			const struct const_Multiarray_d*const geom_coef_fg =
				constructor_coef_fg_vol_opp_dom(si_neigh,s_face); // destructed
			const struct const_Vector_i*const coef_inds =
				constructor_cc0_vgc_fgc_indices(si_neigh,s_face); // destructed

			const struct Solver_Volume*const s_vol = (struct Solver_Volume*) a_s_vol[si_neigh];
			update_rows_Multiarray_d((struct Multiarray_d*)s_vol->geom_coef,geom_coef_fg,coef_inds);

			destructor_const_Multiarray_d(geom_coef_fg);
			destructor_const_Vector_i(coef_inds);
		}

		for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
			struct Adaptive_Solver_Volume*const a_s_vol = (struct Adaptive_Solver_Volume*) curr;
			if (a_s_vol->updated_geom_nc_face)
				a_s_vol->updated_geom_nc = true;
		}
	}
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

/// \brief Set the members of \ref Face::Neigh_Info for the newly created internal face.
static void set_Neigh_Info_Face_i_new
	(struct Face*const face,                             ///< The new face.
	 const int ind_vh,                                   ///< The 'ind'ex of the child 'h'-refined 'v'olume.
	 const int ind_lf,                                   ///< Local face index.
	 const struct Adaptive_Solver_Volume*const a_s_vol_p ///< The parent volume.
	);

/** \brief Get the pointer to the appropriate \ref Geometry_Element::cc0_vgc_fgc operator.
 *  \return See brief. */
static const struct Operator* get_operator__cc0_vgc_fgc
	(const int side_index,                 ///< The index of the side of the face under consideration.
	 const struct Solver_Face*const s_face ///< The current \ref Solver_Face_T.
	);

/** \brief Constructor for the 'f'ace 'g'eometry coefficients of the geometry order of the volume having input
 *         side_index with specified additional parameters.
 *  \return See brief. */
static const struct const_Multiarray_d* constructor_coef_fg_vol_specified
	(const int side_index,                 ///< Side index of the volume for which the coefficients are computed.
	 const int side_index_d,               /**< Side index of the volume from which to interpolate the coefs to the
	                                        *   face. */
	 const bool use_full_face,             ///< Flag for whether the full face should be used.
	 const struct Solver_Face*const s_face ///< \ref Solver_Face_T.
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

static void constructor_Solver_Face_i_new
	(struct Solver_Face*const s_face, const struct Adaptive_Solver_Volume*const a_s_vol_p,
	 struct Flux_Input*const flux_i, const struct Simulation*const sim)
{
	constructor_derived_Solver_Face((struct Face*)s_face,sim);

	const struct Solver_Volume*const s_vol_p = (struct Solver_Volume*) a_s_vol_p;
	const_cast_i(&s_face->p_ref,compute_face_reference_order(s_face));
	const_cast_i(&s_face->ml,s_vol_p->ml+1);
	const_cast_c(&s_face->cub_type,(check_for_curved_neigh((struct Face*)s_face) ? 'c' : 's'));
	set_function_pointers_face_num_flux(s_face,sim);

	switch (sim->method) {
	case METHOD_DG:
		break; // do nothing.
	case METHOD_DPG:
		constructor_Solver_Face__nf_coef(s_face,flux_i,sim);

		const struct Test_Case*const test_case = (struct Test_Case*) sim->test_case_rc->tc;
		assert(!test_case->has_2nd_order); // add constructor s_coef.
		break;
	default:
		EXIT_ERROR("Unsupported: %d.",sim->method);
		break;
	}
}

static const struct const_Multiarray_d* constructor_coef_fg_vol
	(const int side_index, const struct Solver_Face*const s_face)
{
	const int side_index_d = get_dominant_geom_vol_side_index(s_face);
	const struct Face*const face           = (struct Face*) s_face;
	const struct Solver_Volume*const s_vol = (struct Solver_Volume*) face->neigh_info[side_index].volume;
	const bool use_full_face = ( s_vol->ml != s_face->ml ? true : false );

	return constructor_coef_fg_vol_specified(side_index,side_index_d,use_full_face,s_face);
}

static const struct const_Multiarray_d* constructor_coef_fg_vol_opp_dom
	(const int side_index, const struct Solver_Face*const s_face)
{
	const int side_index_d = (side_index+1) % 2;
	const bool use_full_face = true;

	return constructor_coef_fg_vol_specified(side_index,side_index_d,use_full_face,s_face);
}

static const struct const_Vector_i* constructor_cc0_vgc_fgc_indices
	(const int side_index, const struct Solver_Face*const s_face)
{
	const struct Operator*const cc0_vgc_fgc_op = get_operator__cc0_vgc_fgc(side_index,s_face);
	return constructor_const_Vector_i_inds_eq_1_const_Matrix_d(cc0_vgc_fgc_op->op_std);
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
 *  \return See brief.
 *
 *  \todo Add reference with visualization of the returned values.
 */
static int get_ind_child
	(const int ind_h,               ///< The index of the h-refinement of the face.
	 const int side_index,          ///< The index of the side of the face under consideration.
	 const struct Face*const face_p ///< The parent \ref Face.
	);

/** \brief Return the value of \ref Face::Neigh_Info::ind_lf for the desired \ref Face::neigh_info for an h-refined
 *         face.
 *  \return See brief. */
static int get_ind_lf_h_ref
	(const int side_index,          ///< The index of the side of the face under consideration.
	 const struct Face*const face_p ///< The parent \ref Face.
	);

/** \brief Return the value of \ref Face::Neigh_Info::ind_href for the non-dominant \ref Face::neigh_info for an
 *         h-refined face.
 *  \return See brief.
 *
 *  \note The non-dominant \ref Face::neigh_info is that which has index [1].
 */
static int get_ind_href_h_ref_non_dominant
	(const int ind_h,                ///< The index of the h-refinement of the face.
	 const struct Face*const face_p, ///< The parent \ref Face.
	 const struct Face*const face    ///< The current \ref Face.
	);

/** \brief Return the value of \ref Face::Neigh_Info::ind_ord for the desired \ref Face::neigh_info for an h-refined
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

/** \brief Return a pointer to a statically allocated \ref Info_I container corresponding to the current internal face.
 *  \return See brief. */
const struct Info_I* get_Info_I
	(const int ind_vh,                                   ///< The 'ind'ex of the child 'h'-refined 'v'olume.
	 const int ind_lf,                                   ///< Local face index.
	 const struct Adaptive_Solver_Volume*const a_s_vol_p ///< The parent volume.
	);

/** \brief Get the pointer to the appropriate \ref Geometry_Element::vv0_fgc_fgc operator.
 *  \return See brief. */
static const struct Operator* get_operator__vv0_fgc_fgc
	(const int side_index,                 ///< The index of the side of the face under consideration.
	 const struct Solver_Face*const s_face ///< The current \ref Solver_Face_T.
	);

/** \brief Get the pointer to the appropriate \ref Geometry_Element::vc0_fgc_fgc operator.
 *  \return See brief. */
static const struct Operator* get_operator__vc0_fgc_fgc
	(const int side_index,                 ///< The index of the side of the face under consideration.
	 const struct Solver_Face*const s_face ///< The current \ref Solver_Face_T.
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
		constructor_mm_NN1_Operator_const_Multiarray_d(vv0_vv_vv,xyz_ve_p2,'R','d',2,NULL)); // keep
	vol->h = compute_h_volume(vol);

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
	face->neigh_info[0].ind_sref = 0;
	face->neigh_info[0].ind_ord  = get_ind_ord_h_ref(0,face_p);
	face->neigh_info[0].volume   = (struct Volume*) advance_Link(ind_child_0,a_s_vol_p_0->child_0);
	face->h = compute_h_face(face);

	if (!face->boundary) {
		const struct Volume*const vol_p_1                     = face_p->neigh_info[1].volume;
		const struct Adaptive_Solver_Volume*const a_s_vol_p_1 = (struct Adaptive_Solver_Volume*) vol_p_1;

		if (a_s_vol_p_1->adapt_type == ADAPT_H_REFINE) { // neighbour is also being h-refined
			const int ind_child_1 = get_ind_child(ind_h,1,face_p);
			face->neigh_info[1].ind_lf   = get_ind_lf_h_ref(1,face_p);
			face->neigh_info[1].ind_ord  = get_ind_ord_h_ref(1,face_p);
			face->neigh_info[1].ind_href = get_ind_href_h_ref_non_dominant(ind_h,face_p,face);
			face->neigh_info[1].ind_sref = 0;
			face->neigh_info[1].volume   = (struct Volume*) advance_Link(ind_child_1,a_s_vol_p_1->child_0);
		} else { // neighbour is not being refined
			const int ind_compound_1 = get_ind_compound(ind_h,1,face_p);

			face->neigh_info[1].ind_lf   = face_p->neigh_info[1].ind_lf,
			face->neigh_info[1].ind_href = ind_compound_1 % NFREFMAX;
			face->neigh_info[1].ind_sref = 0;
			face->neigh_info[1].ind_ord  = get_ind_ord_h_ref(1,face_p);
			face->neigh_info[1].volume   = face_p->neigh_info[1].volume;
		}
	} else {
		face->neigh_info[1].ind_lf   = -1;
		face->neigh_info[1].ind_href = -1;
		face->neigh_info[1].ind_sref = -1;
		face->neigh_info[1].ind_ord  = -1;
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
	const_cast_i(&s_face->p_ref,compute_face_reference_order(s_face));
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
	const struct Info_I*const info_i = get_Info_I(ind_vh,ind_lf,a_s_vol_p);

	int side_index = 0;
	face->neigh_info[side_index].ind_lf   = ind_lf;
	face->neigh_info[side_index].ind_href = 0;
	face->neigh_info[side_index].ind_sref = 0;
	face->neigh_info[side_index].ind_ord  = info_i->ind_ord[0];
	face->neigh_info[side_index].volume   = (struct Volume*) advance_Link(ind_vh,a_s_vol_p->child_0);
	const_cast_Face(&face->neigh_info[side_index].volume->faces[face->neigh_info[side_index].ind_lf][0],face);

	side_index = 1;
	face->neigh_info[side_index].ind_lf   = info_i->ind_lf_1;
	face->neigh_info[side_index].ind_href = 0;
	face->neigh_info[side_index].ind_sref = 0;
	face->neigh_info[side_index].ind_ord  = info_i->ind_ord[0];
	face->neigh_info[side_index].volume   = (struct Volume*) advance_Link(info_i->ind_vh_1,a_s_vol_p->child_0);
	const_cast_Face(&face->neigh_info[side_index].volume->faces[face->neigh_info[side_index].ind_lf][0],face);
}

static const struct const_Multiarray_d* constructor_coef_fg_vol_specified
	(const int side_index, const int side_index_d, const bool use_full_face, const struct Solver_Face*const s_face)
{
	const struct Face*const face           = (struct Face*) s_face;
	const struct Solver_Volume*const s_vol = (struct Solver_Volume*) face->neigh_info[side_index].volume;

	const struct const_Multiarray_d*const geom_fg_min =
		constructor_geom_fg(side_index_d,side_index,s_face,true,use_full_face); // destructed/moved

	const int p_fg = compute_face_geometry_order(s_face);
	const struct const_Multiarray_d* geom_fg = NULL;
	if (s_vol->p_ref == p_fg) {
		geom_fg = geom_fg_min; // destructed
	} else {
		const struct Operator*const vv0_fgc_fgc = get_operator__vv0_fgc_fgc(side_index,s_face);
		geom_fg = constructor_mm_NN1_Operator_const_Multiarray_d(vv0_fgc_fgc,geom_fg_min,'C','d',2,NULL); // d.
		destructor_const_Multiarray_d(geom_fg_min);
	}

	const struct Operator*const vc0_fgc_fgc = get_operator__vc0_fgc_fgc(side_index,s_face);
	const struct const_Multiarray_d*const geom_coef_fg =
		constructor_mm_NN1_Operator_const_Multiarray_d(vc0_fgc_fgc,geom_fg,'R','d',2,NULL); // returned
	destructor_const_Multiarray_d(geom_fg);

	return geom_coef_fg;
}

static const struct Operator* get_operator__cc0_vgc_fgc (const int side_index, const struct Solver_Face*const s_face)
{
	const struct Face*const face            = (struct Face*) s_face;
	const struct Volume*const vol           = face->neigh_info[side_index].volume;
	const struct Solver_Volume*const s_vol  = (struct Solver_Volume*) vol;
	const struct Geometry_Element*const g_e = &((struct Solver_Element*)vol->element)->g_e;

	const bool use_full_face = ( s_vol->ml != s_face->ml ? true : false );

	const int ind_lf   = face->neigh_info[side_index].ind_lf,
	          ind_href = ( use_full_face ? 0 : face->neigh_info[side_index].ind_href );
	const int p_v = s_vol->p_ref;

	const int curved = ( (s_face->cub_type == 's') ? 0 : 1 );
	assert(curved);

	return get_Multiarray_Operator(g_e->cc0_vgc_fgc,(ptrdiff_t[]){ind_lf,ind_href,0,p_v,p_v});
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
	case TRI:
		switch (ind_h) {
		case H_TRI4_V0:
			bc_faces->data[1] = vol_p->bc_faces->data[1];
			bc_faces->data[2] = vol_p->bc_faces->data[2];
			break;
		case H_TRI4_V1:
			bc_faces->data[0] = vol_p->bc_faces->data[0];
			bc_faces->data[2] = vol_p->bc_faces->data[2];
			break;
		case H_TRI4_V2:
			bc_faces->data[0] = vol_p->bc_faces->data[0];
			bc_faces->data[1] = vol_p->bc_faces->data[1];
			break;
		case H_TRI4_V3:
			break; // do nothing.
		default:
			EXIT_ERROR("Unsupported: %d",ind_h);
			break;
		}
		break;
	case QUAD:
		switch (ind_h) {
		case H_QUAD4_V0:
			bc_faces->data[0] = vol_p->bc_faces->data[0];
			bc_faces->data[2] = vol_p->bc_faces->data[2];
			break;
		case H_QUAD4_V1:
			bc_faces->data[1] = vol_p->bc_faces->data[1];
			bc_faces->data[2] = vol_p->bc_faces->data[2];
			break;
		case H_QUAD4_V2:
			bc_faces->data[0] = vol_p->bc_faces->data[0];
			bc_faces->data[3] = vol_p->bc_faces->data[3];
			break;
		case H_QUAD4_V3:
			bc_faces->data[1] = vol_p->bc_faces->data[1];
			bc_faces->data[3] = vol_p->bc_faces->data[3];
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

	const int ind_href = face_p->neigh_info[side_index].ind_href;

	int ind_child = -1;
	switch (vol_p->element->type) {
	case LINE:
		assert(ind_href == H_POINT1_V0);
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
		switch (ind_href) {
		case H_LINE1_V0:
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
		case H_LINE2_V0:
			switch (ind_compound) {
			case 1*NFREFMAX+H_LINE2_V0: // fallthrough
			case 1*NFREFMAX+H_LINE2_V1: // fallthrough
			case 2*NFREFMAX+H_LINE2_V0: // fallthrough
			case 2*NFREFMAX+H_LINE2_V1:
				ind_child = H_TRI4_V0; break;
			case 0*NFREFMAX+H_LINE2_V0: // fallthrough
			case 0*NFREFMAX+H_LINE2_V1:
				ind_child = H_TRI4_V1; break;
			default:
				EXIT_ERROR("Unsupported: %d",ind_compound);
				break;
			}
			break;
		case H_LINE2_V1:
			switch (ind_compound) {
			case 2*NFREFMAX+H_LINE2_V0: // fallthrough
			case 2*NFREFMAX+H_LINE2_V1:
				ind_child = H_TRI4_V1; break;
			case 0*NFREFMAX+H_LINE2_V0: // fallthrough
			case 0*NFREFMAX+H_LINE2_V1: // fallthrough
			case 1*NFREFMAX+H_LINE2_V0: // fallthrough
			case 1*NFREFMAX+H_LINE2_V1:
				ind_child = H_TRI4_V2; break;
			default:
				EXIT_ERROR("Unsupported: %d",ind_compound);
				break;
			}
			break;
		default:
			EXIT_ERROR("Unsupported: %d",ind_href);
			break;
		}
		break;
	case QUAD:
		switch (ind_href) {
		case H_LINE1_V0:
			switch (ind_compound) {
			case 0*NFREFMAX+H_LINE2_V0: // fallthrough
			case 2*NFREFMAX+H_LINE2_V0:
				ind_child = H_QUAD4_V0; break;
			case 1*NFREFMAX+H_LINE2_V0: // fallthrough
			case 2*NFREFMAX+H_LINE2_V1:
				ind_child = H_QUAD4_V1; break;
			case 0*NFREFMAX+H_LINE2_V1: // fallthrough
			case 3*NFREFMAX+H_LINE2_V0:
				ind_child = H_QUAD4_V2; break;
			case 1*NFREFMAX+H_LINE2_V1: // fallthrough
			case 3*NFREFMAX+H_LINE2_V1:
				ind_child = H_QUAD4_V3; break;
			default:
				EXIT_ERROR("Unsupported: %d",ind_compound);
				break;
			}
			break;
		case H_LINE2_V0:
			switch (ind_compound) {
			case 0*NFREFMAX+H_LINE2_V0: // fallthrough
			case 0*NFREFMAX+H_LINE2_V1: // fallthrough
			case 2*NFREFMAX+H_LINE2_V0: // fallthrough
			case 2*NFREFMAX+H_LINE2_V1:
				ind_child = H_QUAD4_V0; break;
			case 1*NFREFMAX+H_LINE2_V0: // fallthrough
			case 1*NFREFMAX+H_LINE2_V1:
				ind_child = H_QUAD4_V1; break;
			case 3*NFREFMAX+H_LINE2_V0: // fallthrough
			case 3*NFREFMAX+H_LINE2_V1:
				ind_child = H_QUAD4_V2; break;
			default:
				EXIT_ERROR("Unsupported: %d",ind_compound);
				break;
			}
			break;
		case H_LINE2_V1:
			switch (ind_compound) {
			case 2*NFREFMAX+H_LINE2_V0: // fallthrough
			case 2*NFREFMAX+H_LINE2_V1:
				ind_child = H_QUAD4_V1; break;
			case 0*NFREFMAX+H_LINE2_V0: // fallthrough
			case 0*NFREFMAX+H_LINE2_V1:
				ind_child = H_QUAD4_V2; break;
			case 1*NFREFMAX+H_LINE2_V0: // fallthrough
			case 1*NFREFMAX+H_LINE2_V1: // fallthrough
			case 3*NFREFMAX+H_LINE2_V0: // fallthrough
			case 3*NFREFMAX+H_LINE2_V1:
				ind_child = H_QUAD4_V3; break;
			default:
				EXIT_ERROR("Unsupported: %d",ind_compound);
				break;
			}
			break;
		default:
			EXIT_ERROR("Unsupported: %d",ind_href);
			break;
		}
		break;
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

static int get_ind_href_h_ref_non_dominant
	(const int ind_h, const struct Face*const face_p, const struct Face*const face)
{
	const int side_index = 1;
	const int ind_href_p = face_p->neigh_info[side_index].ind_href,
	          ind_ord    = face->neigh_info[side_index].ind_ord;

	const struct Volume*const vol_p = face_p->neigh_info[side_index].volume;
	switch (vol_p->element->type) {
	case LINE:
		/** +1 as the refined POINT is the same as the original point and there is consequently no increment for
		 *  the refined POINT in \ref definitions_h_ref.h. */
		assert(ind_h == H_POINT1_V0+1);
		return H_POINT1_V0;
		break;
	case TRI:  // fallthrough
	case QUAD:
		assert(ind_ord == 0 || ind_ord == 1);
		assert(ind_h == H_LINE2_V0 || ind_h == H_LINE2_V1);
		if (ind_ord == 0) {
			switch (ind_href_p) {
				case H_LINE1_V0:
					return H_LINE1_V0;
					break;
				case H_LINE2_V0: // fallthrough
				case H_LINE2_V1:
					return ind_h;
					break;
				default: EXIT_ERROR("Unsupported: %d",ind_href_p); break;
			}
		} else {
			switch (ind_href_p) {
				case H_LINE1_V0:
					return H_LINE1_V0;
					break;
				case H_LINE2_V0: // fallthrough
				case H_LINE2_V1:
					return ( ind_h == H_LINE2_V0 ? H_LINE2_V1 : H_LINE2_V0 );
					break;
				default: EXIT_ERROR("Unsupported: %d",ind_href_p); break;
			}
		}
		break;
	default:
		EXIT_ERROR("Unsupported: %d",vol_p->element->type);
		break;
	}
	EXIT_ERROR("Should not have reached this point.\n");
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

	if (side_index == 0 || ind_ord == 0 || ind_ord == -1)
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

const struct Info_I* get_Info_I
	(const int ind_vh, const int ind_lf, const struct Adaptive_Solver_Volume*const a_s_vol_p)
{
	const struct Info_I* info_i = NULL;
	const struct Volume*const vol_p = (struct Volume*) a_s_vol_p;

	/**
	 *  \note The values of the members can be generated for each of the supported element/h-refiment combinations
	 *        using the [h_refinement_info.py] script.
	 *
	 *  <!-- References: -->
	 *  [h_refinement_info.py]: h_refinement/h_refinement_info.py
	 */
	switch (vol_p->element->type) {
	case LINE: {
		assert(ind_vh == 0); // Should already have found internal faces for other volumes.
		assert(ind_lf == 1); // Only internal face for volume 0.

		static const struct Info_I i01 = { .ind_vh_1 = 1, .ind_lf_1 = 0, .ind_ord[0] = 0, .ind_ord[1] = 0, };
		info_i = &i01;
		break;
	} case TRI: {
		static const struct Info_I i00 = { .ind_vh_1 = 3, .ind_lf_1 = 0, .ind_ord[0] = 1, .ind_ord[1] = 1, };
		static const struct Info_I i11 = { .ind_vh_1 = 3, .ind_lf_1 = 1, .ind_ord[0] = 1, .ind_ord[1] = 1, };
		static const struct Info_I i22 = { .ind_vh_1 = 3, .ind_lf_1 = 2, .ind_ord[0] = 1, .ind_ord[1] = 1, };
		switch (ind_vh) {
		case 0:
			assert(ind_lf == 0); // Only internal face for volume 0.
			info_i = &i00;
			break;
		case 1:
			assert(ind_lf == 1); // Only internal face for volume 1.
			info_i = &i11;
			break;
		case 2:
			assert(ind_lf == 2); // Only internal face for volume 2.
			info_i = &i22;
			break;
		case 3: // fallthrough
		default:
			EXIT_ERROR("Should already have found all internal faces (%d, %d).",ind_vh,ind_lf);
			break;
		}
		break;
	} case QUAD: {
		static const struct Info_I i01 = { .ind_vh_1 = 1, .ind_lf_1 = 0, .ind_ord[0] = 0, .ind_ord[1] = 0, };
		static const struct Info_I i03 = { .ind_vh_1 = 2, .ind_lf_1 = 2, .ind_ord[0] = 0, .ind_ord[1] = 0, };
		static const struct Info_I i13 = { .ind_vh_1 = 3, .ind_lf_1 = 2, .ind_ord[0] = 0, .ind_ord[1] = 0, };
		static const struct Info_I i21 = { .ind_vh_1 = 3, .ind_lf_1 = 0, .ind_ord[0] = 0, .ind_ord[1] = 0, };
		switch (ind_vh) {
		case 0:
			switch (ind_lf) {
				case 1: info_i = &i01; break;
				case 3: info_i = &i03; break;
				default: EXIT_ERROR("Unsupported: %d.\n",ind_lf); break;
			}
			break;
		case 1:
			assert(ind_lf == 3); // Only internal face for volume 1.
			info_i = &i13;
			break;
		case 2:
			assert(ind_lf == 1); // Only internal face for volume 2.
			info_i = &i21;
			break;
		case 3: // fallthrough
		default:
			EXIT_ERROR("Should already have found all internal faces (%d, %d).",ind_vh,ind_lf);
			break;
		}
		break;
	} default:
		EXIT_ERROR("Unsupported: %d.",vol_p->element->type);
		break;
	}
	return info_i;
}

static const struct Operator* get_operator__vv0_fgc_fgc (const int side_index, const struct Solver_Face*const s_face)
{
	const struct Face*const face            = (struct Face*) s_face;
	const struct Volume* vol                = face->neigh_info[side_index].volume;
	const struct Solver_Volume*const s_vol  = (struct Solver_Volume*) vol;
	const struct Geometry_Element*const g_e = &((struct Solver_Element*)vol->element)->g_e;

	const int ind_e = get_face_element_index(face),
	          p_fg  = compute_face_geometry_order(s_face),
	          p_f    = s_vol->p_ref;
	const int curved = ( (s_face->cub_type == 's') ? 0 : 1 );
	assert(curved == 1);

	return get_Multiarray_Operator(g_e->vv0_fgc_fgc,(ptrdiff_t[]){ind_e,ind_e,0,0,p_f,p_fg});
}

static const struct Operator* get_operator__vc0_fgc_fgc (const int side_index, const struct Solver_Face*const s_face)
{
	const struct Face*const face            = (struct Face*) s_face;
	const struct Volume*const vol           = face->neigh_info[side_index].volume;
	const struct Solver_Volume*const s_vol  = (struct Solver_Volume*) vol;
	const struct Geometry_Element*const g_e = &((struct Solver_Element*)vol->element)->g_e;

	const int ind_e  = get_face_element_index(face),
	          p_f    = s_vol->p_ref;
	const int curved = ( (s_face->cub_type == 's') ? 0 : 1 );
	assert(curved == 1);

	return get_Multiarray_Operator(g_e->vc0_fgc_fgc,(ptrdiff_t[]){ind_e,ind_e,0,0,p_f,p_f});
}
