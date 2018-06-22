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
#include <math.h>

#include "macros.h"
#include "definitions_tol.h"

#include "def_templates_face_solver_opg.h"
#include "def_templates_volume_solver_opg.h"

#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"

#include "def_templates_compute_face_rlhs_opg.h"
#include "def_templates_compute_rlhs.h"
#include "def_templates_math_functions.h"
#include "def_templates_numerical_flux.h"
#include "def_templates_penalty_opg_advection.h"
#include "def_templates_face_solver.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void constructor_rlhs_f_test_penalty_advection_upwind_T
	(const struct Flux_Ref_T*const flux_r, const struct Numerical_Flux_T*const num_flux,
	 struct Solver_Face_T*const s_face, struct Solver_Storage_Implicit*const ssi)
{
	/* printf("returning without penalty\n"); */
	/* return; */
	UNUSED(flux_r);

	/** It is assumed that the initial solution for \ref Solver_Volume_T::test_s_coef is equal to zero on all
	 *  boundary faces which require the addition of the penalty term (outflow faces). It is also currently assumed
	 *  that the exact solution test functions are equal to zero (i.e. that the value of `g` described in
	 *  \ref constructor_rlhs_f_b_test_penalty_T is equal to zero).
	 *
	 *  While this is consistent with the test space used by Brunken et al. \cite Brunken2018, it has the
	 *  disadvantage of the exact test functions having jump discontinuities on faces which have a combination of
	 *  inflow and outflow boundary conditions. It would likely be desirable to set the value of the boundary
	 *  condition for the test functions based on the values of the solution test functions interpolated to the
	 *  inflow boundary nodes. \todo Investigate when applicable.
	 *
	 *  It is also possible that imposing the boundary conditions for the solution test functions at the face
	 *  cubature nodes does not set the physically correct number of boundary conditions when the number of cubature
	 *  nodes is not equal to the number of face basis functions. It may be required to set values for the numerical
	 *  solution in a face basis and subsequently interpolate to the face cubature nodes. Note that an identical
	 *  problem may be present for the treatment of the face terms for the DG scheme as well...
	 */
	const int*const n_vr_eq = get_set_n_var_eq(NULL);
	const int n_vr = n_vr_eq[0];
	const int n_eq = n_vr_eq[1];
	assert(n_vr == 1 && n_eq == 1); // Ensure that all is working properly when removed; need to add loop over eq/var.

	const struct Face*const face                     = (struct Face*) s_face;
	const struct Solver_Volume_T*const s_vol         = (struct Solver_Volume_T*) face->neigh_info[0].volume;
	const struct OPG_Solver_Face_T*const opg_s_face  = (struct OPG_Solver_Face_T*) s_face;

	const struct const_Multiarray_T*const n_dot_b = num_flux->neigh_info[0].dnnf_ds;
	assert(num_flux->nnf != NULL);
	assert(n_dot_b != NULL);

	const ptrdiff_t n_fc = n_dot_b->extents[0];
	struct Vector_R*const indicator = constructor_zero_Vector_R(n_fc); // destructed
	for (int n = 0; n < n_fc; ++n) {
		/* if (real_T(n_dot_b->data[n]) > -EPS) */
		if (real_T(n_dot_b->data[n]) > EPS)
		/* if (real_T(n_dot_b->data[n]) > 0.0) */
			indicator->data[n] = 1.0;
	}

#if 0
	static int count = 0;
	const int count_max = 2;
	/* const int count_max = (int) round(2.0/face->h); // number of faces along one boundary face */
	/* const int count_max = (int) round(4.0/face->h); // number of faces along two boundary faces */
	if (count == count_max)
		return;
	/* if (face->index == 7 || face->index == 10) */
	/* if (face->index == 11 || face->index == 7 || face->index == 10) */
	if (face->index == 6 || face->index == 7 || face->index == 10)
		return;
	if (indicator->data[0]) {
		printf("%d\n",face->index);
		print_const_Multiarray_T(n_dot_b);
		print_const_Multiarray_T(s_face->xyz_fc);
		print_const_Multiarray_T(s_face->normals_fc);
		++count;
	}
//	printf("%d %e\n",count_max,face->h);
#endif

	const struct Lhs_Operators_OPG_T*const ops = constructor_Lhs_Operators_OPG_T(opg_s_face); // destructed

	const int ind_sc = (s_face->cub_type == 'c');
	const int p_t_p = get_set_degree_poly(NULL,"tp")[ind_sc];
	/* const double scale = pow(face->h,s_vol->p_ref+p_t_p); */
	const double scale = 1.0+0.0*pow(face->h,s_vol->p_ref+p_t_p);

	// Note: -ve sign included here. Would need to be moved for g != 0.
	const struct const_Vector_T*const diag = constructor_dot_mult_const_Vector_T_RT
		(-PENALTY_SCALING_OPG/scale,(struct const_Vector_R*)indicator,ops->wJ_fc,1); // destructed
	destructor_Vector_R(indicator);

	const struct const_Matrix_T*const lhs_r =
		constructor_mm_diag_const_Matrix_R_T(1.0,ops->cv0_vt_fc[0],diag,'L',false); // destructed
	destructor_const_Vector_T(diag);

	const struct const_Matrix_T*const lhs =
		constructor_mm_RT_const_Matrix_T('T','N',1.0,ops->cv0_vt_fc[0],lhs_r,'R'); // destructed
	destructor_Lhs_Operators_OPG_T(s_face,ops);

#if TYPE_RC == TYPE_REAL
	// Currently assumed that w == 0 (i.e. always initializing to 0) and using a Lagrange basis.
//	* \f$ \eps^{-1} <v,g-w>_{\Gamma^\text{characteristic out}} \forall v \f$
	struct Vector_T*const g_V = constructor_empty_Vector_T(n_fc); // destructed
	double g = 0.1;
	if (face->index == 6)
		g += 0.2;
	set_to_value_Vector_T(g_V,g);
	struct const_Matrix_T g_M = { .ext_0 = g_V->ext_0, .ext_1 = 1, .layout = 'C', .data = g_V->data, };
	struct Matrix_T rhs_M = interpret_Multiarray_as_Matrix_T(s_vol->rhs);
	mm_T('T','N',-1.0,1.0,lhs_r,&g_M,&rhs_M);

	destructor_Vector_T(g_V);

	if (ssi != NULL) {
		const struct OPG_Solver_Volume_T*const opg_s_vol = (struct OPG_Solver_Volume_T*) face->neigh_info[0].volume;
		for (int vr = 0; vr < n_vr; ++vr) {
		for (int eq = 0; eq < n_eq; ++eq) {
			if (eq != vr)
				continue;
			set_petsc_Mat_row_col_opg(ssi,opg_s_vol,eq,opg_s_vol,vr);
			add_to_petsc_Mat(ssi,lhs);
		}}
	}
#elif TYPE_RC == TYPE_COMPLEX
	assert(ssi == NULL);
	mm_NNC_Multiarray_TTT(1.0,1.0,lhs,(struct const_Multiarray_T*)s_vol->test_s_coef,s_vol->rhs);
#endif
	destructor_const_Matrix_T(lhs_r);
	destructor_const_Matrix_T(lhs);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

#include "undef_templates_face_solver_opg.h"
#include "undef_templates_volume_solver_opg.h"

#include "undef_templates_matrix.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"

#include "undef_templates_compute_face_rlhs_opg.h"
#include "undef_templates_compute_rlhs.h"
#include "undef_templates_math_functions.h"
#include "undef_templates_numerical_flux.h"
#include "undef_templates_penalty_opg_advection.h"
#include "undef_templates_face_solver.h"
