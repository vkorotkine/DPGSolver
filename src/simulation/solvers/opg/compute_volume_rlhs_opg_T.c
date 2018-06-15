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

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"
#include "definitions_intrusive.h"


#include "def_templates_compute_volume_rlhs_opg.h"

#include "def_templates_multiarray.h"

#include "def_templates_volume_solver.h"
#include "def_templates_volume_solver_opg.h"

#include "def_templates_compute_volume_rlhs.h"
#include "def_templates_flux.h"
#include "def_templates_test_case.h"
#include "def_templates_operators.h"

// Static function declarations ************************************************************************************* //

/// \brief Container for solver-related parameters.
struct S_Params_T {
	struct S_Params_Volume_Structor_T spvs; ///< \ref S_Params_Volume_Structor_T.

	compute_rlhs_v_fptr_T compute_rlhs; ///< Pointer to the appropriate function.
};

/** \brief Set the parameters of \ref S_Params_T.
 *  \return A statically allocated \ref S_Params_T container. */
static struct S_Params_T set_s_params_T
	(const struct Simulation* sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void compute_volume_rlhs_opg_T
	(const struct Simulation* sim, struct Solver_Storage_Implicit* ssi, struct Intrusive_List* volumes)
{
	assert(sim->volumes->name == IL_VOLUME_SOLVER_OPG);
	assert(sim->elements->name == IL_ELEMENT_SOLVER_OPG);

	struct S_Params_T s_params = set_s_params_T(sim);
	struct Flux_Input_T* flux_i = constructor_Flux_Input_T(sim); // destructed

	for (struct Intrusive_Link* curr = volumes->first; curr; curr = curr->next) {
		struct Solver_Volume_T*const s_vol = (struct Solver_Volume_T*) curr;

		struct Flux_Ref_T* flux_r = constructor_Flux_Ref_vol_T(&s_params.spvs,flux_i,s_vol,sim);
		print_const_Multiarray_T(s_vol->geom_coef);

		// Compute the rhs and the lhs terms.
		s_params.compute_rlhs(flux_r,s_vol,ssi);
		destructor_Flux_Ref_T(flux_r);
	}
	destructor_Flux_Input_T(flux_i);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct S_Params_T set_s_params_T (const struct Simulation*const sim)
{
	struct S_Params_T s_params;

	set_S_Params_Volume_Structor_T(&s_params.spvs,sim);

	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	switch (test_case->solver_method_curr) {
#if TYPE_RC == TYPE_COMPLEX
	case 'e':
		s_params.compute_rlhs = compute_rhs_v_dg_like_T;
		break;
#elif TYPE_RC == TYPE_REAL
	case 'i':
		if (test_case->has_1st_order && !test_case->has_2nd_order)
			s_params.compute_rlhs = compute_rlhs_1;
		else if (!test_case->has_1st_order && test_case->has_2nd_order)
			EXIT_ADD_SUPPORT; //s_params.compute_rlhs = compute_rlhs_2;
		else if (test_case->has_1st_order && test_case->has_2nd_order)
			EXIT_ADD_SUPPORT; //s_params.compute_rlhs = compute_rlhs_12;
		else
			EXIT_ADD_SUPPORT;
		break;
#endif
	default:
		EXIT_ERROR("Unsupported: %c (type_rc: %d)\n",test_case->solver_method_curr,TYPE_RC);
		break;
	}

	return s_params;
}
