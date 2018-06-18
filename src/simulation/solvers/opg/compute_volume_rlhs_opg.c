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

/// \brief Version of \ref compute_rlhs_v_fptr_T computing the rhs and lhs terms for 1st order equations only.
static void compute_rlhs_1
	(const struct Flux_Ref*const flux_r,      ///< See brief.
	 struct Solver_Volume*const s_vol,        ///< See brief.
	 struct Solver_Storage_Implicit*const ssi ///< See brief.
	);

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "compute_volume_rlhs_opg_T.c"
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "compute_volume_rlhs_opg_T.c"
#include "undef_templates_type.h"

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Version of \ref compute_rlhs_v_fptr_T computing the lhs term for linearization wrt the solution.
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
	(const struct Flux_Ref*const flux_r,    ///< Standard.
	 const struct Solver_Volume*const s_vol ///< Standard.
	);

static void compute_lhs_1
	(const struct Flux_Ref*const flux_r, struct OPG_Solver_Volume*const opg_s_vol,
	 struct Solver_Storage_Implicit*const ssi)
{
	assert(get_set_collocated(NULL) == false); // Ensure that premultiplication by inv(w_vc) is not present if true.

	struct Solver_Volume*const s_vol = (struct Solver_Volume*) opg_s_vol;
	const struct const_Matrix_d*const lhs = constructor_lhs_v_1_opg(flux_r,s_vol); // destructed
	set_petsc_Mat_row_col_opg(ssi,opg_s_vol,0,opg_s_vol,0);
	add_to_petsc_Mat(ssi,lhs);
	destructor_const_Matrix_d(lhs);
}

// Level 2 ********************************************************************************************************** //

static const struct const_Matrix_d* constructor_lhs_v_1_opg
	(const struct Flux_Ref*const flux_r, const struct Solver_Volume*const s_vol)
{
	const struct const_Matrix_d*const cv1r = constructor_test_diff_op_1v_opg_d(flux_r,s_vol,false); // destructed

	const int n_vr = get_set_n_var_eq(NULL)[0];

	const struct const_Vector_d*const w_vc = get_operator__w_vc__s_e(s_vol);
	const struct const_Vector_d J_vc       = interpret_const_Multiarray_as_Vector_d(s_vol->jacobian_det_vc);

	const struct const_Vector_d*const wJ_vc =
		constructor_dot_mult_inverse_2nd_const_Vector_d(1.0,w_vc,&J_vc,n_vr); // destructed

	const struct const_Matrix_d*const n1_lt = constructor_mm_diag_const_Matrix_d_d(1.0,cv1r,wJ_vc,'L',false); // dest.
	destructor_const_Vector_d(wJ_vc);

	// -ve sign from the definition of u == - test_diff_op*test_s_coef.
	const struct const_Matrix_d*const lhs = constructor_mm_const_Matrix_d('T','N',-1.0,n1_lt,cv1r,'R'); // returned
	destructor_const_Matrix_d(n1_lt);
	destructor_const_Matrix_d(cv1r);

	return lhs;
}
