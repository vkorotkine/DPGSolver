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
	assert(sim->volumes->name == IL_VOLUME_SOLVER);
	assert(sim->faces->name   == IL_FACE_SOLVER);
	assert(list_is_derived_from("solver",'e',sim));

	constructor_derived_computational_elements(sim,IL_SOLVER_ADAPTIVE); // destructed

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
	(struct Adaptive_Solver_Volume** a_s_vol, ///< Pointer to pointer of current \ref Adaptive_Solver_Volume.
	 const struct Simulation*const sim        ///< \ref Simulation.
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

/// \brief Project all \ref Solver_Volume_T solution-related data to the p-adaptive volume.
static void compute_projection_p_volume
	(struct Adaptive_Solver_Volume* a_s_vol, ///< \ref Adaptive_Solver_Volume.
	 const struct Simulation* sim            ///< \ref Simulation.
	);

/// \brief Project all \ref Solver_Face_T solution-related data to the p-adaptive face.
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
		case ADAPT_H_REFINE: // fallthrough
		case ADAPT_H_COARSE:
			update_volume_h(&a_s_vol,sim);
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
		default:
			EXIT_ERROR("Unsupported: %d\n",adapt_type);
			break;
		}

	}
}

// Level 2 ********************************************************************************************************** //

/// \brief Constructs the children of the current volume and updates the volume list.
static void constructor_volumes_h_refine
	(struct Adaptive_Solver_Volume** a_s_vol, ///< Pointer to pointer of current \ref Adaptive_Solver_Volume.
	 const struct Simulation*const sim        ///< \ref Simulation.
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
}

static void update_volume_h (struct Adaptive_Solver_Volume** a_s_vol, const struct Simulation*const sim)
{
	const struct Solver_Volume* s_vol = (struct Solver_Volume*) *a_s_vol;
	const int ml = s_vol->ml;

	const int adapt_type = (*a_s_vol)->adapt_type;
	switch (adapt_type) {
	case ADAPT_H_REFINE:
		constructor_volumes_h_refine(a_s_vol,sim);
		assert(ml < sim->ml[1]);
		break;
	case ADAPT_H_COARSE:
		assert(ml > sim->ml[0]);
		break;
	case ADAPT_P_REFINE: // fallthrough
	case ADAPT_P_COARSE: // fallthrough
	default:
		EXIT_ERROR("Unsupported: %d\n",adapt_type);
		break;
	}
EXIT_ADD_SUPPORT;
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

#include "test_case.h"
static void compute_projection_p_face (struct Adaptive_Solver_Face* a_s_face, const struct Simulation* sim)
{
	const struct Face* face    = (struct Face*) a_s_face;
	struct Solver_Face* s_face = (struct Solver_Face*) a_s_face;
	const struct Volume* vol   = face->neigh_info[0].volume;
	const struct Adaptation_Element* a_e = &((struct Solver_Element*)vol->element)->a_e;

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

	assert(((struct Test_Case*)sim->test_case_rc->tc)->has_2nd_order == false); // Add support. Remove include above.
}

// Level 3 ********************************************************************************************************** //

/** \brief Return the number of children associated with the isotropically h-refined volume of input type.
 *  \return See brief. */
static int get_n_children
	(const struct Adaptive_Solver_Volume*const a_s_vol ///< The current volume.
	);

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

static void constructor_volumes_h_refine
	(struct Adaptive_Solver_Volume** a_s_vol, const struct Simulation*const sim)
{
	struct Adaptive_Solver_Volume* a_s_vol_p = *a_s_vol;
	const struct Volume*const vol_p          = (struct Volume*) a_s_vol_p;

	const struct Solver_Element* s_e     = (struct Solver_Element*) vol_p->element;
	const struct Adaptation_Element* a_e = &s_e->a_e;
	const struct Operator*const vv0_vv_vv = get_Multiarray_Operator(a_e->vv0_vv_vv,(ptrdiff_t[]){0,0,2,1});

	const struct const_Multiarray_d*const xyz_ve_p2 =
		constructor_mm_NN1_Operator_const_Multiarray_d(vv0_vv_vv,vol_p->xyz_ve,'C','d',2,NULL); // destructed

	struct Intrusive_List* volumes_c = constructor_empty_IL(IL_VOLUME_SOLVER_ADAPTIVE,NULL); // destructed

print_const_Multiarray_d(xyz_ve_p2);
	const int n_children = get_n_children(*a_s_vol);
	for (int n = 1; n <= n_children; ++n) {

		push_back_IL(volumes_c,(struct Intrusive_Link*)
			constructor_Adaptive_Solver_Volume_h_ref(n,a_s_vol_p,xyz_ve_p2,sim));

		if (n == 1) {
			a_s_vol_p->child_0 = (struct Adaptive_Solver_Volume*) volumes_c->first;
		}

	}
	destructor_const_Multiarray_d(xyz_ve_p2);
	destructor_IL(volumes_c,false);
EXIT_UNSUPPORTED;
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

/** \brief Constructor for the \ref Solver_Volume_T base of the h-refined \ref Adaptive_Solver_Volume.
 *  \return See brief. */
static void constructor_Solver_Volume_h_ref
	(struct Solver_Volume*const s_vol,                    ///< The child volume.
	 const int ind_h,                                     ///< The index of the child h-refined element.
	 const struct Adaptive_Solver_Volume*const a_s_vol_p, ///< The parent volume.
	 const struct Simulation*const sim                    ///< \ref Simulation.
	);

static int get_n_children (const struct Adaptive_Solver_Volume*const a_s_vol)
{
	struct Volume* vol = (struct Volume*) a_s_vol;
	switch (vol->element->type) {
	case LINE:  return 2;  break;
	case TRI:   return 4;  break;
	case QUAD:  return 4;  break;
	case TET:   return 8;  break;
	case HEX:   return 8;  break;
	case WEDGE: return 8;  break;
	case PYR:   return 10; break;
	default:
		EXIT_ERROR("Unsupported: %d",vol->element->type);
		break;
	}
}

static struct Adaptive_Solver_Volume* constructor_Adaptive_Solver_Volume_h_ref
	(const int ind_h, const struct Adaptive_Solver_Volume*const a_s_vol_p,
	 const struct const_Multiarray_d*const xyz_ve_p2_i, const struct Simulation*const sim)
{
	struct Adaptive_Solver_Volume*const a_s_vol = malloc(sizeof *a_s_vol); // returned

	constructor_Volume_h_ref((struct Volume*)a_s_vol,ind_h,a_s_vol_p,xyz_ve_p2_i,sim);
	constructor_Solver_Volume_h_ref((struct Solver_Volume*)a_s_vol,ind_h,a_s_vol_p,sim);

	const struct Solver_Volume*const s_vol_p = (struct Solver_Volume*) a_s_vol_p;
	a_s_vol->adapt_type = ADAPT_H;
	a_s_vol->p_ref_prev = s_vol_p->p_ref;
	a_s_vol->updated    = true;

	return a_s_vol;
}

// Level 5 ********************************************************************************************************** //

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

static void constructor_Volume_h_ref
	(struct Volume*const vol, const int ind_h, const struct Adaptive_Solver_Volume*const a_s_vol_p,
	 const struct const_Multiarray_d*const xyz_ve_p2_i, const struct Simulation*const sim)
{
	const struct Volume*const vol_p          = (struct Volume*) a_s_vol_p;
	const struct Solver_Volume*const s_vol_p = (struct Solver_Volume*) a_s_vol_p;
	const_cast_i(&vol->index,-1);

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

static void constructor_Solver_Volume_h_ref
	(struct Solver_Volume*const s_vol, const int ind_h, const struct Adaptive_Solver_Volume*const a_s_vol_p,
	 const struct Simulation*const sim)
{
	struct Volume*const vol                  = (struct Volume*) s_vol;
	const struct Volume*const vol_p          = (struct Volume*) a_s_vol_p;
	const struct Solver_Volume*const s_vol_p = (struct Solver_Volume*) a_s_vol_p;

	constructor_derived_Solver_Volume(vol,sim);

	const_cast_i(&s_vol->p_ref,s_vol_p->p_ref);
	const_cast_i(&s_vol->ml,s_vol_p->ml+1);

	compute_geometry_volume(s_vol,sim);

	const struct Solver_Element* s_e     = (struct Solver_Element*) vol_p->element;
	const struct Adaptation_Element* a_e = &s_e->a_e;

	const int p = s_vol->p_ref;
	if (compute_size(s_vol_p->sol_coef->order,s_vol_p->sol_coef->extents) != 0) {
		const struct Operator*const cc0_vs_vs = get_Multiarray_Operator(a_e->cc0_vs_vs,(ptrdiff_t[]){ind_h,0,p,p});

		destructor_Multiarray_d(s_vol->sol_coef);
		s_vol->sol_coef = constructor_mm_NN1_Operator_Multiarray_d(
			cc0_vs_vs,s_vol_p->sol_coef,'C','d',2,NULL); // keep
	}

	if (compute_size(s_vol_p->grad_coef->order,s_vol_p->grad_coef->extents) != 0) {
		EXIT_ADD_SUPPORT;
	}

	if (test_case_explicitly_enforces_conservation(sim))
		set_Multiarray_d(s_vol->l_mult,(struct const_Multiarray_d*)s_vol_p->l_mult);
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
