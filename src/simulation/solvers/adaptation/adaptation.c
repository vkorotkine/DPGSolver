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

#include "macros.h"
#include "definitions_adaptation.h"
#include "definitions_intrusive.h"

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

// Interface functions ********************************************************************************************** //

void adapt_hp (struct Simulation* sim, const int adapt_strategy)
{
	constructor_derived_computational_elements(sim,IL_SOLVER_ADAPTIVE); // destructed

	mark_volumes_to_adapt(sim,adapt_strategy);
	adapt_hp_volumes(sim);

	destructor_derived_computational_elements(sim,IL_SOLVER);
EXIT_UNSUPPORTED;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Update \ref Solver_Volume::p_ref.
static void update_volume_p_ref
	(struct Adaptive_Solver_Volume* a_s_vol, ///< \ref Adaptive_Solver_Volume.
	 const struct Simulation* sim            ///< \ref Simulation.
	);

/// \brief Project all \ref Solver_Volume solution-related data to the p-adaptved volume.
static void compute_projection_p_volume
	(struct Adaptive_Solver_Volume* a_s_vol, ///< \ref Adaptive_Solver_Volume.
	 const struct Simulation* sim            ///< \ref Simulation.
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
/// \todo Separate functions when usage is set.
	// Update volume flags
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Adaptive_Solver_Volume* a_s_vol = (struct Adaptive_Solver_Volume*) curr;

		const int adapt_type = a_s_vol->adapt_type;
		switch (adapt_type) {
		case ADAPT_NONE:
			break; // Do nothing.
		case ADAPT_P_REFINE: // fallthrough
		case ADAPT_P_COARSE:
			update_volume_p_ref(a_s_vol,sim);
			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",adapt_type);
			break;
		}
	}

	// Update volume geometry
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

	// Project volume solution
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
EXIT_UNSUPPORTED;
}

// Level 1 ********************************************************************************************************** //

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

print_const_Matrix_d(cc0_vs_vs->op_std);
print_Multiarray_d(s_coef_p);
print_Multiarray_d(s_coef);
EXIT_UNSUPPORTED;
}
