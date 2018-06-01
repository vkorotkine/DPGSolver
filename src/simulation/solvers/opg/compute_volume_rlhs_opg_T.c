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

/** \brief Function pointer to the function used to evaluate the rhs and lhs terms.
 *
 *  \param flux_r    \ref Flux_Ref_T.
 *  \param opg_s_vol \ref OPG_Solver_Volume_T.
 *  \param ssi       \ref Solver_Storage_Implicit.
 *  \param sim       \ref Simulation.
 */
typedef void (*compute_rlhs_opg_fptr_T)
	(const struct Flux_Ref_T*const flux_r,
	 struct OPG_Solver_Volume_T*const opg_s_vol,
	 struct Solver_Storage_Implicit*const ssi,
	 const struct Simulation*const sim
	);

/// \brief Container for solver-related parameters.
struct S_Params_T {
	struct S_Params_Volume_Structor_T spvs; ///< \ref S_Params_Volume_Structor_T.

	compute_rlhs_opg_fptr_T compute_rlhs; ///< Pointer to the appropriate function.
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
		struct Solver_Volume_T* s_vol         = (struct Solver_Volume_T*) curr;
		struct OPG_Solver_Volume_T* opg_s_vol = (struct OPG_Solver_Volume_T*) curr;

		struct Flux_Ref_T* flux_r = constructor_Flux_Ref_vol_T(&s_params.spvs,flux_i,s_vol,sim);

		// Compute the rhs and the lhs terms.
		s_params.compute_rlhs(flux_r,opg_s_vol,ssi,sim);
		destructor_Flux_Ref_T(flux_r);
	}
	destructor_Flux_Input_T(flux_i);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Version of \ref compute_rlhs_opg_fptr_T computing only the rhs term.
static void compute_rhs_all
	(const struct Flux_Ref_T*const flux_r,       ///< See brief.
	 struct OPG_Solver_Volume_T*const opg_s_vol, ///< See brief.
	 struct Solver_Storage_Implicit*const ssi,   ///< See brief.
	 const struct Simulation*const sim           ///< See brief.
	);

static struct S_Params_T set_s_params_T (const struct Simulation*const sim)
{
	struct S_Params_T s_params;

	set_S_Params_Volume_Structor_T(&s_params.spvs,sim);

	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	switch (test_case->solver_method_curr) {
#if TYPE_RC == TYPE_COMPLEX
	case 'e':
		s_params.compute_rlhs = compute_rhs_all;
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

// Level 1 ********************************************************************************************************** //

static void compute_rhs_all
	(const struct Flux_Ref_T*const flux_r, struct OPG_Solver_Volume_T*const opg_s_vol,
	 struct Solver_Storage_Implicit*const ssi, const struct Simulation* sim)
{
	UNUSED(sim);

	struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) opg_s_vol;
	const struct Multiarray_Operator tw1_vt_vc = get_operator__tw1_vt_vc_T(s_vol);

	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	const char op_format = 'd';

	const ptrdiff_t extents[] = { tw1_vt_vc.data[0]->op_std->ext_0, flux_r->fr->extents[1], };
	struct Multiarray_T*const rhs = constructor_zero_Multiarray_T('C',2,extents); // destructed;
	for (ptrdiff_t dim = 0; dim < DIM; ++dim)
		mm_NNC_Operator_Multiarray_T(1.0,1.0,tw1_vt_vc.data[dim],flux_r->fr,rhs,op_format,2,&dim,NULL);

	UNUSED(ssi);
	// Try to use the same function as DG to avoid duplication (think about face/source terms).
	EXIT_UNSUPPORTED; // Add -ve rhs to ssi->b.
}
