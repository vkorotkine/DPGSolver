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

#include "test_complex_compute_face_rhs.h"

#include "test_complex_boundary.h"
#include "test_complex_numerical_flux.h"

#include <assert.h>

#include "macros.h"

#include "complex_matrix.h"
#include "complex_multiarray.h"
#include "complex_vector.h"
#include "matrix.h"
#include "multiarray.h"

#include "compute_face_rlhs.h"
#include "operator.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void destructor_Numerical_Flux_Input_c_data (struct Numerical_Flux_Input_c* num_flux_i)
{
	destructor_Boundary_Value_Input_c(&num_flux_i->bv_l);
	destructor_Boundary_Value_c(&num_flux_i->bv_r);
}

struct Matrix_c* constructor_lhs_f_1_c
	(const int side_index[2], const struct Numerical_Flux_c* num_flux, const struct Solver_Face* s_face)
{
	const struct Operator* tw0_vt_fc    = get_operator__tw0_vt_fc(side_index[0],s_face);
	const struct Operator* cv0_vs_fc_op = get_operator__cv0_vs_fc(side_index[1],s_face);

	const struct const_Matrix_d* cv0_vs_fc = cv0_vs_fc_op->op_std;
	bool need_free_cv0 = false;
	if (side_index[0] != side_index[1]) {
		need_free_cv0 = true;
		cv0_vs_fc = constructor_copy_const_Matrix_d(cv0_vs_fc); // destructed
		permute_Matrix_d_fc((struct Matrix_d*)cv0_vs_fc,'R',side_index[0],s_face);
	}

	const struct const_Multiarray_c* dnnf_ds_Ma = num_flux->neigh_info[side_index[1]].dnnf_ds;

	const ptrdiff_t ext_0 = tw0_vt_fc->op_std->ext_0,
	                ext_1 = tw0_vt_fc->op_std->ext_1;
	const int n_eq = dnnf_ds_Ma->extents[1],
	          n_vr = dnnf_ds_Ma->extents[2];

	struct Matrix_c* tw0_nf = constructor_empty_Matrix_c('R',ext_0,ext_1);                         // destructed
	struct Matrix_c* lhs_l  = constructor_empty_Matrix_c('R',ext_0,cv0_vs_fc->ext_1);              // destructed
	struct Matrix_c* lhs    = constructor_empty_Matrix_c('R',n_eq*lhs_l->ext_0,n_vr*lhs_l->ext_1); // returned
	set_to_value_Matrix_c(tw0_nf,0.0);

	struct Vector_c dnnf_ds = { .ext_0 = dnnf_ds_Ma->extents[0], .owns_data = false, .data = NULL, };

	for (int vr = 0; vr < n_vr; ++vr) {
	for (int eq = 0; eq < n_eq; ++eq) {
		const ptrdiff_t ind =
			compute_index_sub_container(dnnf_ds_Ma->order,1,dnnf_ds_Ma->extents,(ptrdiff_t[]){eq,vr});
		dnnf_ds.data = (double complex*)&dnnf_ds_Ma->data[ind];
		mm_diag_c('R',1.0,0.0,tw0_vt_fc->op_std,(struct const_Vector_c*)&dnnf_ds,tw0_nf,false);

		mm_cdc('N','N',-1.0,0.0,(struct const_Matrix_c*)tw0_nf,cv0_vs_fc,lhs_l);

		set_block_Matrix_c(lhs,(struct const_Matrix_c*)lhs_l,eq*lhs_l->ext_0,vr*lhs_l->ext_1,'i');
	}}
	destructor_Matrix_c(tw0_nf);
	destructor_Matrix_c(lhs_l);

	if (need_free_cv0)
		destructor_const_Matrix_d(cv0_vs_fc);

	return lhs;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
