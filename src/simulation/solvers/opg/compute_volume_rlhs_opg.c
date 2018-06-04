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

#include "compute_volume_rlhs_opg.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "face_solver_opg.h"
#include "volume_solver_opg.h"
#include "element_solver_opg.h"

#include "compute_rlhs.h"
#include "compute_volume_rlhs.h"
#include "flux.h"
#include "intrusive.h"
#include "math_functions.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "solve.h"
#include "solve_opg.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Version of \ref compute_rlhs_opg_fptr_T computing the rhs and lhs terms for 1st order equations only.
static void compute_rlhs_1
	(const struct Flux_Ref*const flux_r,      ///< See brief.
	 struct Solver_Volume*const s_vol,        ///< See brief.
	 struct Solver_Storage_Implicit*const ssi ///< See brief.
	);

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "compute_volume_rlhs_opg_T.c"

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Version of \ref compute_rlhs_opg_fptr_T computing the lhs term for linearization wrt the solution.
static void compute_lhs_1
	(const struct Flux_Ref*const flux_r,       ///< See brief.
	 struct OPG_Solver_Volume*const opg_s_vol, ///< See brief.
	 struct Solver_Storage_Implicit*const ssi  ///< See brief.
	);

static void compute_rlhs_1
	(const struct Flux_Ref*const flux_r, struct Solver_Volume*const s_vol, struct Solver_Storage_Implicit*const ssi)
{
	struct OPG_Solver_Volume*const opg_s_vol = (struct OPG_Solver_Volume*) s_vol;
	compute_rhs_v_dg_like(flux_r,s_vol,ssi);
	compute_lhs_1(flux_r,opg_s_vol,ssi);
}

// Level 1 ********************************************************************************************************** //

/** \brief Constructor for the volume contribution to the lhs term for 1st order equations.
 *  \return See brief. */
static const struct const_Matrix_d* constructor_lhs_v_1_opg
	(const struct Flux_Ref*const flux_r,            ///< Standard.
	 const struct OPG_Solver_Volume*const opg_s_vol ///< Standard.
	);

static void compute_lhs_1
	(const struct Flux_Ref*const flux_r, struct OPG_Solver_Volume*const opg_s_vol,
	 struct Solver_Storage_Implicit*const ssi)
{
	assert(get_set_collocated(NULL) == false); // Ensure that premultiplication by inv(w_vc) is not present if true.

	const struct const_Matrix_d*const lhs = constructor_lhs_v_1_opg(flux_r,opg_s_vol); // destructed
	set_petsc_Mat_row_col_opg(ssi,opg_s_vol,0,opg_s_vol,0);
	add_to_petsc_Mat(ssi,lhs);
	destructor_const_Matrix_d(lhs);
}

// Level 2 ********************************************************************************************************** //

static const struct const_Matrix_d* constructor_lhs_v_1_opg
	(const struct Flux_Ref*const flux_r, const struct OPG_Solver_Volume*const opg_s_vol)
{
	const int*const n_var_eq = get_set_n_var_eq(NULL);
	const int n_vr = n_var_eq[0],
	          n_eq = n_var_eq[1];

	struct Solver_Volume* s_vol = (struct Solver_Volume*) opg_s_vol;
	const struct Multiarray_Operator cv1_vt_vc = get_operator__cv1_vt_vc(s_vol);

	const ptrdiff_t ext_0 = cv1_vt_vc.data[0]->op_std->ext_0,
	                ext_1 = cv1_vt_vc.data[0]->op_std->ext_1;

	struct Matrix_d* cv1r = constructor_empty_Matrix_d('R',n_vr*ext_0,n_eq*ext_1); // destructed

	struct Matrix_d* cv1r_l = constructor_empty_Matrix_d('R',ext_0,ext_1); // destructed
	const struct const_Multiarray_d* dfr_ds_Ma = flux_r->dfr_ds;
	struct Vector_d dfr_ds = { .ext_0 = dfr_ds_Ma->extents[0], .owns_data = false, .data = NULL, };

	for (int vr = 0; vr < n_vr; ++vr) {
	for (int eq = 0; eq < n_eq; ++eq) {
		set_to_value_Matrix_d(cv1r_l,0.0);
		for (int dim = 0; dim < DIM; ++dim) {
			const ptrdiff_t ind =
				compute_index_sub_container(dfr_ds_Ma->order,1,dfr_ds_Ma->extents,(ptrdiff_t[]){eq,vr,dim});
			dfr_ds.data = (Type*)&dfr_ds_Ma->data[ind];
			mm_diag_d('L',1.0,1.0,cv1_vt_vc.data[dim]->op_std,(struct const_Vector_d*)&dfr_ds,cv1r_l,false);
		}
		set_block_Matrix_d(cv1r,eq*ext_0,vr*ext_1,
		                   (struct const_Matrix_d*)cv1r_l,0,0,cv1r_l->ext_0,cv1r_l->ext_1,'i');
	}}
	destructor_Matrix_d(cv1r_l);

	const struct const_Vector_d* w_vc = get_operator__w_vc__s_e(s_vol);
	const struct const_Vector_d J_vc  = interpret_const_Multiarray_as_Vector_d(s_vol->jacobian_det_vc);

	const struct const_Vector_d* J_inv_vc = constructor_inverse_const_Vector_d(&J_vc);                   // destructed
	const struct const_Vector_d* wJ_vc    = constructor_dot_mult_const_Vector_d(1.0,w_vc,J_inv_vc,n_vr); // destructed
	destructor_const_Vector_d(J_inv_vc);

	const struct const_Matrix_d* n1_lt =
		constructor_mm_diag_const_Matrix_d_d(1.0,(struct const_Matrix_d*)cv1r,wJ_vc,'L',false); // destructed
	destructor_const_Vector_d(wJ_vc);

	const struct const_Matrix_d*const lhs =
		constructor_mm_const_Matrix_d('T','N',1.0,n1_lt,(struct const_Matrix_d*)cv1r,'R'); // returned
	destructor_const_Matrix_d(n1_lt);
	destructor_Matrix_d(cv1r);

	return lhs;
}
