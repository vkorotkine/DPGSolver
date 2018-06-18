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

#include "test_support_solve.h"
#include "test_support_computational_elements.h"
#include "test_support_multiarray.h"

#include <assert.h>

#include "macros.h"
#include "definitions_elements.h"
#include "definitions_intrusive.h"
#include "definitions_test_case.h"
#include "definitions_test_integration.h"

#include "face_solver_dg.h"
#include "face_solver_dpg.h"
#include "face_solver_opg.h"
#include "volume_solver_dg.h"
#include "volume_solver_dpg.h"
#include "volume_solver_opg.h"

#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void perturb_solution (const struct Simulation*const sim)
{
	switch (get_set_method(NULL)) {
	case METHOD_DG:
		break; // do nothing
	case METHOD_DPG: // fallthrough
	case METHOD_OPG:
		assert(get_set_has_1st_2nd_order(NULL)[1] == false); // Add support.
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",get_set_method(NULL));
		break;
	}

	// template this if maintenance becomes prohibitive.
	if (sim->test_case_rc->is_real) {
		for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
			struct Solver_Volume* s_vol = (struct Solver_Volume*) curr;
			perturb_Multiarray_d(s_vol->sol_coef,MAX_PERTURB);
		}

		for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
			struct Solver_Face* s_face = (struct Solver_Face*) curr;
			perturb_Multiarray_d(s_face->nf_coef,MAX_PERTURB);
		}

		if (test_case_explicitly_enforces_conservation(sim)) {
			for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
				struct Solver_Volume* s_vol = (struct Solver_Volume*) curr;
				perturb_Multiarray_d(s_vol->l_mult,1e3*MAX_PERTURB);
			}
		}
	} else {
		for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
			struct Solver_Volume_c* s_vol = (struct Solver_Volume_c*) curr;
			perturb_Multiarray_c(s_vol->sol_coef,MAX_PERTURB);
		}

		for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
			struct Solver_Face_c* s_face = (struct Solver_Face_c*) curr;
			perturb_Multiarray_c(s_face->nf_coef,MAX_PERTURB);
		}

		if (test_case_explicitly_enforces_conservation(sim)) {
			for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
				struct Solver_Volume_c* s_vol = (struct Solver_Volume_c*) curr;
				perturb_Multiarray_c(s_vol->l_mult,1e3*MAX_PERTURB);
			}
		}
	}
}

struct Intrusive_List* constructor_Volumes_local_centre_only (const struct Volume*const vol)
{
	struct Data { int list_name; size_t sizeof_base; } data;
	switch (get_set_method(NULL)) {
	case METHOD_DG:
		data.list_name   = IL_VOLUME_SOLVER_DG;
		data.sizeof_base = sizeof(struct DG_Solver_Volume_c);
		break;
	case METHOD_DPG:
		data.list_name   = IL_VOLUME_SOLVER_DPG;
		data.sizeof_base = sizeof(struct DPG_Solver_Volume_c);
		break;
	case METHOD_OPG:
		data.list_name   = IL_VOLUME_SOLVER_OPG;
		data.sizeof_base = sizeof(struct OPG_Solver_Volume_c);
		break;
	default:
		EXIT_ERROR("Unsupported: %d.\n",get_set_method(NULL));
		break;
	}

	struct Intrusive_List*const volumes = constructor_empty_IL(data.list_name,NULL); // returned

	struct Intrusive_Link* curr = (struct Intrusive_Link*) vol;

	// A copy is required such that the link in the global list is not modified.
	push_back_IL(volumes,constructor_copied_Intrusive_Link(curr,data.sizeof_base,data.sizeof_base));
	return volumes;
}

struct Intrusive_List* constructor_Faces_local_neigh_only
	(const struct Volume*const vol, const struct Simulation*const sim)
{
	struct Data { int list_name; size_t sizeof_base; } data;
	switch (get_set_method(NULL)) {
	case METHOD_DG:
		data.list_name   = IL_FACE_SOLVER_DG;
		data.sizeof_base = sizeof(struct DG_Solver_Face_c);
		break;
	case METHOD_DPG:
		data.list_name   = IL_FACE_SOLVER_DPG;
		data.sizeof_base = sizeof(struct DPG_Solver_Face_c);
		break;
	case METHOD_OPG:
		data.list_name   = IL_FACE_SOLVER_OPG;
		data.sizeof_base = sizeof(struct OPG_Solver_Face_c);
		break;
	default:
		EXIT_ERROR("Unsupported: %d.\n",get_set_method(NULL));
		break;
	}

	struct Intrusive_List*const faces = constructor_empty_IL(data.list_name,NULL);

	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		if (!is_face_neighbour((struct Face*)curr,vol))
			continue;

		// A copy is required such that the link in the global list is not modified.
		push_back_IL(faces,constructor_copied_Intrusive_Link(curr,data.sizeof_base,data.sizeof_base));
	}
	return faces;
}

bool is_face_neighbour (const struct Face*const face, const struct Volume*const vol_curr)
{
	for (int f = 0; f < NFMAX; ++f) {
	for (int sf = 0; sf < NSUBFMAX; ++sf) {
		const struct Face*const face_n = vol_curr->faces[f][sf];

		if (!face_n)
			continue;

		if (face->index == face_n->index)
			return true;
	}}
	return false;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
