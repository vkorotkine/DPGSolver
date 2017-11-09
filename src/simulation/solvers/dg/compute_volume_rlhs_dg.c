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

#include "compute_volume_rlhs_dg.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"
#include "definitions_intrusive.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "volume_solver_dg.h"
#include "element_solver_dg.h"

#include "compute_volume_rlhs.h"
#include "flux.h"
#include "intrusive.h"
#include "math_functions.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "solve_dg.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Function pointer to the function used to evaluate the rhs (and optionally lhs) terms.
 *
 *  \param flux_r   \ref Flux_Ref.
 *  \param dg_s_vol \ref DG_Solver_Volume.
 *  \param ssi      \ref Solver_Storage_Implicit.
 *  \param sim      \ref Simulation.
 */
typedef void (*compute_rlhs_dg_fptr)
	(const struct Flux_Ref* flux_r,
	 struct DG_Solver_Volume* dg_s_vol,
	 struct Solver_Storage_Implicit* ssi,
	 const struct Simulation* sim
	);

/// \brief Container for solver-related parameters.
struct S_Params {
	struct S_Params_Volume_Structor spvs; ///< \ref S_Params_Volume_Structor.

	compute_rlhs_dg_fptr compute_rlhs; ///< Pointer to the appropriate function.
};

/** \brief Set the parameters of \ref S_Params.
 *  \return A statically allocated \ref S_Params container. */
static struct S_Params set_s_params
	(const struct Simulation* sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void compute_volume_rlhs_dg (const struct Simulation* sim, struct Solver_Storage_Implicit* ssi)
{
	assert(sim->volumes->name == IL_VOLUME_SOLVER_DG);
	assert(sim->elements->name == IL_ELEMENT_SOLVER_DG);

	struct S_Params s_params = set_s_params(sim);
	struct Flux_Input* flux_i = constructor_Flux_Input(sim); // destructed

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume* s_vol       = (struct Solver_Volume*) curr;
		struct DG_Solver_Volume* dg_s_vol = (struct DG_Solver_Volume*) curr;

		struct Flux_Ref* flux_r  = constructor_Flux_Ref_vol(&s_params.spvs,flux_i,s_vol,sim);

		// Compute the rhs (and optionally the lhs) terms.
		s_params.compute_rlhs(flux_r,dg_s_vol,ssi,sim);
		destructor_Flux_Ref(flux_r);
//EXIT_UNSUPPORTED;
	}
	destructor_Flux_Input(flux_i);
//EXIT_UNSUPPORTED;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Version of \ref compute_rlhs_dg_fptr computing only the rhs term.
static void compute_rhs_v_dg
	(const struct Flux_Ref* flux_r,       ///< See brief.
	 struct DG_Solver_Volume* dg_s_vol,   ///< See brief.
	 struct Solver_Storage_Implicit* ssi, ///< See brief.
	 const struct Simulation* sim         ///< See brief.
	);

/// \brief Version of \ref compute_rlhs_dg_fptr computing the rhs and lhs terms for 1st order equations only.
static void compute_rlhs_1
	(const struct Flux_Ref* flux_r,       ///< See brief.
	 struct DG_Solver_Volume* dg_s_vol,   ///< See brief.
	 struct Solver_Storage_Implicit* ssi, ///< See brief.
	 const struct Simulation* sim         ///< See brief.
	);

static struct S_Params set_s_params (const struct Simulation* sim)
{
	struct S_Params s_params;

	set_S_Params_Volume_Structor(&s_params.spvs,sim);

	struct Test_Case* test_case = sim->test_case;
	switch (test_case->solver_method_curr) {
	case 'e':
		s_params.compute_rlhs = compute_rhs_v_dg;
		break;
	case 'i':
		if (test_case->has_1st_order && !test_case->has_2nd_order)
			s_params.compute_rlhs = compute_rlhs_1;
		else if (!test_case->has_1st_order && test_case->has_2nd_order)
			EXIT_ADD_SUPPORT; //s_params.compute_rlhs = compute_rlhs_2;
		else if (test_case->has_1st_order && test_case->has_2nd_order)
			EXIT_ADD_SUPPORT; //s_params.compute_rlhs = compute_rlhs_12;
		else
			EXIT_ERROR("Unsupported: %d %d\n",test_case->has_1st_order,test_case->has_2nd_order);
		break;
	default:
		EXIT_ERROR("Unsupported: %c\n",test_case->solver_method_curr);
		break;
	}

	return s_params;
}

// Level 1 ********************************************************************************************************** //

static void compute_rhs_v_dg
	(const struct Flux_Ref* flux_r, struct DG_Solver_Volume* dg_s_vol, struct Solver_Storage_Implicit* ssi,
	 const struct Simulation* sim)
{
	UNUSED(ssi);

	struct Solver_Volume* s_vol = (struct Solver_Volume*) dg_s_vol;
	const struct Multiarray_Operator tw1_vt_vc = get_operator__tw1_vt_vc__rlhs(s_vol);

	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	const char op_format = 'd';

	const ptrdiff_t d = sim->d;
	for (ptrdiff_t dim = 0; dim < d; ++dim)
		mm_NNC_Operator_Multiarray_d(1.0,1.0,tw1_vt_vc.data[dim],flux_r->fr,dg_s_vol->rhs,op_format,2,&dim,NULL);
}

static void compute_rlhs_1
	(const struct Flux_Ref* flux_r, struct DG_Solver_Volume* dg_s_vol, struct Solver_Storage_Implicit* ssi,
	 const struct Simulation* sim)
{
/// \todo Add special case for collocated.
// If collocation is enabled, note that the diagonal weight scaling must be added back in to recover the symmetry of the
// residual Jacobian. Add it just before adding the contribution to the petsc mat. Also add for face terms and RHS
// terms (volume, face, source or simply the complete rhs).
assert(sim->collocated == false); // Add support in future.
	const ptrdiff_t d = sim->d;

	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	const char op_format = 'd';

	struct Volume* vol          = (struct Volume*) dg_s_vol;
	struct Solver_Volume* s_vol = (struct Solver_Volume*) vol;

	const struct Solver_Element* e = (const struct Solver_Element*) vol->element;

	const int p      = s_vol->p_ref,
	          curved = vol->curved;
	const struct Multiarray_Operator tw1_vt_vc;
	set_MO_from_MO(&tw1_vt_vc,e->tw1_vt_vc[curved],1,(ptrdiff_t[]){0,0,p,p});

	// rhs
	for (ptrdiff_t dim = 0; dim < d; ++dim)
		mm_NNC_Operator_Multiarray_d(1.0,1.0,tw1_vt_vc.data[dim],flux_r->fr,dg_s_vol->rhs,op_format,2,&dim,NULL);
//print_Multiarray_d(dg_s_vol->rhs);

	// lhs
	const struct Operator* cv0_vs_vc = get_Multiarray_Operator(e->cv0_vs_vc[curved],(ptrdiff_t[]){0,0,p,p});

	const ptrdiff_t ext_0 = tw1_vt_vc.data[0]->op_std->ext_0,
	                ext_1 = tw1_vt_vc.data[0]->op_std->ext_1;

	struct Matrix_d* tw1_r = constructor_empty_Matrix_d('R',ext_0,ext_1);                    // destructed
	struct Matrix_d* lhs   = constructor_empty_Matrix_d('R',ext_0,cv0_vs_vc->op_std->ext_1); // destructed

	const struct const_Multiarray_d* dfr_ds_Ma = flux_r->dfr_ds;
	struct Vector_d dfr_ds = { .ext_0 = dfr_ds_Ma->extents[0], .owns_data = false, .data = NULL, };

	const int n_eq = sim->test_case->n_eq,
	          n_vr = sim->test_case->n_var;
	for (int vr = 0; vr < n_vr; ++vr) {
	for (int eq = 0; eq < n_eq; ++eq) {
		set_to_value_Matrix_d(tw1_r,0.0);
		for (int dim = 0; dim < d; ++dim) {
			const ptrdiff_t ind =
				compute_index_sub_container(dfr_ds_Ma->order,1,dfr_ds_Ma->extents,(ptrdiff_t[]){eq,vr,dim});
			dfr_ds.data = (double*)&dfr_ds_Ma->data[ind];
			mm_diag_d('R',1.0,1.0,tw1_vt_vc.data[dim]->op_std,(struct const_Vector_d*)&dfr_ds,tw1_r,false);
		}

		mm_d('N','N',1.0,0.0,(struct const_Matrix_d*)tw1_r,cv0_vs_vc->op_std,lhs);
//print_Matrix_d(lhs);

		set_petsc_Mat_row_col(ssi,s_vol,eq,s_vol,vr);
		add_to_petsc_Mat(ssi,(struct const_Matrix_d*)lhs);
	}}
	destructor_Matrix_d(tw1_r);
	destructor_Matrix_d(lhs);
}
