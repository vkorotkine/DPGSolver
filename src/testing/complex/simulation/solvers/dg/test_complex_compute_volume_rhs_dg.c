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

#include "test_complex_compute_volume_rhs_dg.h"

#include <assert.h>

#include "macros.h"
#include "definitions_core.h"
#include "definitions_intrusive.h"

#include "test_complex_flux.h"
#include "test_complex_operators.h"
#include "test_complex_compute_volume_rhs.h"

#include "volume_solver_dg_complex.h"

#include "complex_multiarray.h"
#include "multiarray.h"

#include "compute_volume_rlhs.h"
#include "intrusive.h"
#include "multiarray_operator.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

/// \brief `complex` version of \ref compute_rhs_v_dg.
static void compute_rhs_v_dg_c
	(const struct Flux_Ref_c* flux_r,           ///< See brief.
	 struct Complex_DG_Solver_Volume* dg_s_vol, ///< See brief.
	 const struct Simulation* sim               ///< See brief.
	);

// Interface functions ********************************************************************************************** //

void compute_volume_rhs_dg_c (const struct Simulation* sim, struct Intrusive_List* volumes)
{
	assert(sim->volumes->name == IL_VOLUME_SOLVER_DG_COMPLEX);
	assert(sim->elements->name == IL_ELEMENT_SOLVER_DG);

	struct S_Params_Volume_Structor_c spvs = { .constructor_sol_vc = constructor_sol_vc_dg_c, };
	struct Flux_Input_c* flux_i = constructor_Flux_Input_c(sim); // destructed

	for (struct Intrusive_Link* curr = volumes->first; curr; curr = curr->next) {
		struct Solver_Volume* s_vol               = (struct Solver_Volume*) curr;
		struct Complex_DG_Solver_Volume* dg_s_vol = (struct Complex_DG_Solver_Volume*) curr;

		struct Flux_Ref_c* flux_r = constructor_Flux_Ref_vol_c(&spvs,flux_i,s_vol,sim); // destructed

		compute_rhs_v_dg_c(flux_r,dg_s_vol,sim);
		destructor_Flux_Ref_c(flux_r);
	}
	destructor_Flux_Input_c(flux_i);
}

const struct const_Multiarray_c* constructor_sol_vc_dg_c
	(const struct Solver_Volume* s_vol, const struct Simulation* sim)
{
	UNUSED(sim);
	const struct Operator* cv0_vs_vc = get_operator__cv0_vs_vc(s_vol);

	struct Complex_DG_Solver_Volume* s_vol_c = (struct Complex_DG_Solver_Volume*) s_vol;
	const struct const_Multiarray_c* s_coef = (const struct const_Multiarray_c*) s_vol_c->sol_coef;

	return constructor_mm_NN1_Operator_const_Multiarray_c(cv0_vs_vc,s_coef,'C','d',s_coef->order,NULL);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void compute_rhs_v_dg_c
	(const struct Flux_Ref_c* flux_r, struct Complex_DG_Solver_Volume* dg_s_vol, const struct Simulation* sim)
{
	UNUSED(sim);
	struct Solver_Volume* s_vol = (struct Solver_Volume*) dg_s_vol;
	const struct Multiarray_Operator tw1_vt_vc = get_operator__tw1_vt_vc(s_vol);

	for (ptrdiff_t dim = 0; dim < DIM; ++dim)
		mm_NNC_Operator_Multiarray_c(1.0,1.0,tw1_vt_vc.data[dim],flux_r->fr,dg_s_vol->rhs,'d',2,&dim,NULL);
}
