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

#include "compute_all_rlhs_dpg.h"

#include <assert.h>
#include <stdlib.h>

#include "macros.h"
#include "definitions_dpg.h"
#include "definitions_intrusive.h"

#include "face_solver_dpg.h"
#include "volume_solver_dpg.h"
#include "element_solver_dpg.h"

#include "matrix.h"
#include "multiarray.h"

#include "compute_volume_rlhs.h"
#include "flux.h"
#include "intrusive.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Container for solver-related parameters.
struct S_Params {
	struct S_Params_Volume_Structor spvs; ///< \ref S_Params_Volume_Structor.

	// compute_test_norm_operator; ///< Pointer to the appropriate function.
//	compute_rlhs_fptr compute_rlhs; ///< Pointer to the appropriate function.
};

/** \brief Set the parameters of \ref S_Params.
 *  \return A statically allocated \ref S_Params container. */
static struct S_Params set_s_params
	(const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Constructor for the norm operator used to compute the optimal test functions.
 *  \return See brief. */
/// \todo Change doxygen comments.
static const struct const_Matrix_d* constructor_norm_op__h1_upwind
	(const struct DPG_Solver_Volume* dpg_s_vol, ///< \ref The current volume.
	 const struct Flux_Ref* flux_r,             ///< \ref Flux_Ref.
	 const struct Simulation* sim               ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void compute_all_rlhs_dpg (const struct Simulation* sim, struct Solver_Storage_Implicit* ssi)
{
	assert(sim->test_case->solver_method_curr == 'i');
	assert(sim->volumes->name == IL_VOLUME_SOLVER_DPG);
	assert(sim->faces->name == IL_FACE_SOLVER_DPG);
	assert(sim->elements->name == IL_ELEMENT_SOLVER_DPG);

	struct S_Params s_params = set_s_params(sim);
	struct Flux_Input* flux_i = constructor_Flux_Input(sim); // destructed

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
//struct Volume*        vol   = (struct Volume*) curr;
//printf("v_ind: %d\n",vol->index);
		struct Solver_Volume* s_vol         = (struct Solver_Volume*) curr;
		struct DPG_Solver_Volume* dpg_s_vol = (struct DPG_Solver_Volume*) curr;

		struct Flux_Ref* flux_r  = constructor_Flux_Ref_vol(&s_params.spvs,flux_i,s_vol,sim);

/// \todo Think whether it would be better to store the test_norm operator as part of the dpg solver volume.
		const struct const_Matrix_d* norm_op = constructor_norm_op__h1_upwind(dpg_s_vol,flux_r,sim); // destructed
EXIT_ERROR("Change this to function pointer.\n");

		destructor_const_Matrix_d(norm_op);
// Compute_op_norm
EXIT_UNSUPPORTED;

		// Compute the rhs and the lhs terms.
//		s_params.compute_rlhs(flux_r,vol,ssi,sim);
		destructor_Flux_Ref(flux_r);
UNUSED(ssi);
EXIT_UNSUPPORTED;
	}
	destructor_Flux_Input(flux_i);
//EXIT_UNSUPPORTED;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct S_Params set_s_params (const struct Simulation* sim)
{
	struct S_Params s_params;

	set_S_Params_Volume_Structor(&s_params.spvs,sim);

	struct Test_Case* test_case = sim->test_case;
	switch (test_case->solver_method_curr) {
	case 'i':
		if (test_case->has_1st_order && !test_case->has_2nd_order)
//			s_params.compute_rlhs = compute_rlhs_1;
;
		else if (!test_case->has_1st_order && test_case->has_2nd_order)
			EXIT_ADD_SUPPORT; //s_params.compute_rlhs = compute_rlhs_2;
		else if (test_case->has_1st_order && test_case->has_2nd_order)
			EXIT_ADD_SUPPORT; //s_params.compute_rlhs = compute_rlhs_12;
		else
			EXIT_ERROR("Unsupported: %d %d\n",test_case->has_1st_order,test_case->has_2nd_order);
		break;
	case 'e': // fallthrough
	default:
		EXIT_ERROR("Unsupported: %c\n",test_case->solver_method_curr);
		break;
	}

	return s_params;
}

static const struct const_Matrix_d* constructor_norm_op__h1_upwind
	(const struct DPG_Solver_Volume* dpg_s_vol, const struct Flux_Ref* flux_r, const struct Simulation* sim)
{
// function for this.
	const struct Multiarray_Operator cvt1_vt_vc;
	set_MO_from_MO(&cvt1_vt_vc,e->cvt1_vt_vc[curved],1,(ptrdiff_t[]){0,0,p,p});

	const int n_eq = sim->test_case->n_eq,
	          n_vr = sim->test_case->n_var;

// cv't'ransposed
	const ptrdiff_t ext_0 = cvt1_vt_vc.data[0]->op_std->ext_0,
	                ext_1 = cvt1_vt_vc.data[0]->op_std->ext_1;

EXIT_UNSUPPORTED;
	struct Matrix_d* cvt1_r[n_eq] = { NULL, };
	for (int i = 0; i < n_eq; ++i)
		cvt1_r[i] = constructor_empty_Matrix_d('R',ext_0,ext_1); // destructed

	struct Matrix_d* norm = constructor_empty_Matrix_d('R',n_eq*ext_0,n_eq*ext_0); // destructed

	const struct const_Multiarray_d* dfr_ds_Ma = flux_r->dfr_ds;
	struct Vector_d dfr_ds = { .ext_0 = dfr_ds_Ma->extents[0], .owns_data = false, .data = NULL, };

	for (int vr = 0; vr < n_vr; ++vr) {
	for (int eq = 0; eq < n_eq; ++eq) {
		set_to_value_Matrix_d(cvt1_r,0.0);
		for (int dim = 0; dim < d; ++dim) {
			const ptrdiff_t ind =
				compute_index_sub_container(dfr_ds_Ma->order,1,dfr_ds_Ma->extents,(ptrdiff_t[]){eq,vr,dim});
			dfr_ds.data = (double*)&dfr_ds_Ma->data[ind];
			mm_diag_d('R',1.0,1.0,cvt1_vt_vc.data[dim]->op_std,(struct const_Vector_d*)&dfr_ds,cvt1_r,false);
		}

		mm_diag_d('R',1.0,1.0,cvt1_r,(struct const_Vector_d*)&dfr_ds,cvt1_r,false);
		mm_d('N','N',1.0,0.0,(struct const_Matrix_d*)cvt1_r,cv0_vs_vc->op_std,lhs);
//print_Matrix_d(lhs);

		set_petsc_Mat_row_col(ssi,s_vol,eq,s_vol,vr);
		add_to_petsc_Mat(ssi,(struct const_Matrix_d*)lhs);
	}}
	destructor_Matrix_d(cvt1_r);
	destructor_Matrix_d(lhs);
}





//	const struct Volume* vol = (struct Volume*) dpg_s_vol;

	// compute cv1_vt_vc (dot) dfr_ds (Use analogue of dg function).
	// compute w/det_J
UNUSED(dpg_s_vol);
UNUSED(sim);
EXIT_ADD_SUPPORT;
return NULL;
}
