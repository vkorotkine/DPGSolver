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
#include "definitions_intrusive.h"

#include "face_solver_adaptive.h"
#include "element_adaptation.h"
#include "volume_solver_adaptive.h"

#include "matrix.h"
#include "multiarray.h"

#include "computational_elements.h"
#include "const_cast.h"
#include "geometry.h"
#include "intrusive.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "solve.h"

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

// Interface functions ********************************************************************************************** //

void adapt_hp (struct Simulation* sim, const int adapt_strategy)
{
	constructor_derived_computational_elements(sim,IL_SOLVER_ADAPTIVE); // destructed
/** \todo Remove the constructor/destructor of derived elements within the solver functions.
 *  Likely directly derive the adaptive solver element at the start of the solve and include all required adaptation
 *  operators in a single element. Or form an adaptive solver element from multiple other elements and access each
 *  through the appropriate pointer.
 */

	mark_volumes_to_adapt(sim,adapt_strategy);
	adapt_hp_volumes(sim);
	adapt_hp_faces(sim);
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

/// \brief Update geometry for \ref Solver_Volume\*s.
static void update_geometry_volumes
	(struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Update geometry for \ref Solver_Face\*s.
static void update_geometry_faces
	(struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Project the solution for \ref Solver_Volume\*s.
static void project_solution_volumes
	(struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Project the solution for \ref Solver_Face\*s.
static void project_solution_faces
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
	default:
		EXIT_ERROR("Unsupported: %d\n",adapt_strategy);
		break;
	}
}

static void adapt_hp_volumes (struct Simulation* sim)
{
	update_hp_members_volumes(sim);
	update_geometry_volumes(sim);
	project_solution_volumes(sim);
}

static void adapt_hp_faces (struct Simulation* sim)
{
	update_hp_members_faces(sim);
	update_geometry_faces(sim);
	project_solution_faces(sim);
}

// Level 1 ********************************************************************************************************** //

/// \brief Update \ref Solver_Volume::p_ref.
static void update_volume_p_ref
	(struct Adaptive_Solver_Volume* a_s_vol, ///< \ref Adaptive_Solver_Volume.
	 const struct Simulation* sim            ///< \ref Simulation.
	);

/// \brief Update \ref Solver_Face::p_ref.
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

/// \brief Project all \ref Solver_Volume solution-related data to the p-adaptive volume.
static void compute_projection_p_volume
	(struct Adaptive_Solver_Volume* a_s_vol, ///< \ref Adaptive_Solver_Volume.
	 const struct Simulation* sim            ///< \ref Simulation.
	);

/// \brief Project all \ref Solver_Face solution-related data to the p-adaptive face.
static void compute_projection_p_face
	(struct Adaptive_Solver_Face* a_s_face, ///< \ref Adaptive_Solver_Face.
	 const struct Simulation* sim           ///< \ref Simulation.
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

		switch (adapt_type) {
		case ADAPT_P_REFINE: // fallthrough
		case ADAPT_P_COARSE:
			update_face_p_ref(a_s_face,sim);
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

static void update_geometry_volumes (struct Simulation* sim)
{
	constructor_derived_Elements(sim,IL_ELEMENT_GEOMETRY);
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Adaptive_Solver_Volume* a_s_vol = (struct Adaptive_Solver_Volume*) curr;
		struct Solver_Volume* s_vol            = (struct Solver_Volume*) curr;

		const int adapt_type = a_s_vol->adapt_type;
		switch (adapt_type) {
		case ADAPT_NONE:
			break; // Do nothing.
		case ADAPT_P_REFINE: // fallthrough
		case ADAPT_P_COARSE:
			compute_geometry_volume(s_vol,sim);
			break;
		case ADAPT_H_REFINE:
			EXIT_ADD_SUPPORT; // Loop over child volumes and set geometry.
		default:
			EXIT_ERROR("Unsupported: %d\n",adapt_type);
			break;
		}
	}
	destructor_derived_Elements(sim,IL_ELEMENT);
}

static void update_geometry_faces (struct Simulation* sim)
{
	constructor_derived_Elements(sim,IL_ELEMENT_GEOMETRY);
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
	destructor_derived_Elements(sim,IL_ELEMENT);
}

static void project_solution_volumes (struct Simulation* sim)
{
	constructor_derived_Elements(sim,IL_ELEMENT_ADAPTATION);
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
		default:
			EXIT_ERROR("Unsupported: %d\n",adapt_type);
			break;
		}

	}
	destructor_derived_Elements(sim,IL_ELEMENT);
}

static void project_solution_faces (struct Simulation* sim)
{
	constructor_derived_Elements(sim,IL_ELEMENT_ADAPTATION);
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
		default:
			EXIT_ERROR("Unsupported: %d\n",adapt_type);
			break;
		}

	}
	destructor_derived_Elements(sim,IL_ELEMENT);
}

// Level 2 ********************************************************************************************************** //

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
		const struct Solver_Volume* s_vol_l = (struct Solver_Volume*) face->neigh_info[0].volume;
		const struct Solver_Volume* s_vol_r = (struct Solver_Volume*) face->neigh_info[1].volume;
		if (s_vol_l->p_ref != s_vol_r->p_ref)
			return (s_vol_l->p_ref > s_vol_r->p_ref ? 0 : 1);
		else if (s_vol_l->ml != s_vol_r->ml)
			return (s_vol_l->ml > s_vol_r->ml ? 0 : 1);
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
	const struct Face* face          = (struct Face*) a_s_face;

	const int ind_dom_vol = get_dominant_volume_index(a_s_face);
	const struct Adaptive_Solver_Volume* a_s_vol_d =
		(struct Adaptive_Solver_Volume*) face->neigh_info[ind_dom_vol].volume;

	if (!a_s_vol_d->updated)
		a_s_face->adapt_type = ADAPT_NONE;

	const struct Solver_Face* s_face    = (struct Solver_Face*) a_s_face;
	const struct Solver_Volume* s_vol_d = (struct Solver_Volume*) a_s_vol_d;
	if (s_vol_d->p_ref > s_face->p_ref)
		a_s_face->adapt_type = ADAPT_P_REFINE;
	else if (s_vol_d->p_ref < s_face->p_ref)
		a_s_face->adapt_type = ADAPT_P_COARSE;
	else if (s_vol_d->ml > s_face->ml)
		a_s_face->adapt_type = ADAPT_H_REFINE;
	else if (s_vol_d->ml < s_face->ml)
		a_s_face->adapt_type = ADAPT_H_COARSE;
	else
		EXIT_ERROR("Did not find the adaptation type.");
}

static void compute_projection_p_volume (struct Adaptive_Solver_Volume* a_s_vol, const struct Simulation* sim)
{
	const struct Volume* vol    = (struct Volume*) a_s_vol;
	struct Solver_Volume* s_vol = (struct Solver_Volume*) a_s_vol;
	const struct Adaptation_Element* a_e = (struct Adaptation_Element*) vol->element;

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

#include "test_case.h"
static void compute_projection_p_face (struct Adaptive_Solver_Face* a_s_face, const struct Simulation* sim)
{
	const struct Face* face    = (struct Face*) a_s_face;
	struct Solver_Face* s_face = (struct Solver_Face*) a_s_face;
	const struct Volume* vol   = face->neigh_info[0].volume;
	const struct Adaptation_Element* a_e = (struct Adaptation_Element*) vol->element;

	const struct Multiarray_d* nf_coef_p = s_face->nf_coef;
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

	assert(sim->test_case->has_2nd_order == false); // Add support for 2nd order. Remove include above.
}
