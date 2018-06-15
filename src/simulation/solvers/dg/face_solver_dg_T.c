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

#include "macros.h"

#include "def_templates_face_solver.h"
#include "def_templates_face_solver_dg.h"
#include "def_templates_volume_solver.h"

#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"

#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void constructor_derived_DG_Solver_Face_T (struct Face* face_ptr, const struct Simulation* sim)
{
	const struct Face*const face            = face_ptr;
	struct DG_Solver_Face_T*const dg_s_face = (struct DG_Solver_Face_T*) face_ptr;

	const struct Test_Case_T*const test_case = (struct Test_Case_T*) sim->test_case_rc->tc;
	if (test_case->has_2nd_order) {
		const int n_neigh = ( face->boundary ? 1 : 2 );
		for (int i = 0; i < n_neigh; ++i) {
			const struct Solver_Volume_T*const s_vol = (struct Solver_Volume_T*) face->neigh_info[i].volume;
			const int order = s_vol->grad_coef->order;
			ptrdiff_t* extents = s_vol->grad_coef->extents;

			struct Neigh_Info_DG*const ni = &dg_s_face->neigh_info[i];
			ni->grad_coef_f = constructor_zero_Multiarray_T('C',order,extents); // destructed

			if (test_case->solver_method_curr == 'i') {
				for (int i = 0; i < DIM; ++i) {
					ni->d_g_coef_f__d_s_coef[0][i] = constructor_empty_const_Matrix_T('R',0,0); // destructed
					ni->d_g_coef_f__d_s_coef[1][i] = constructor_empty_const_Matrix_T('R',0,0); // destructed
				}
			} else {
				assert(test_case->solver_method_curr == 'e');
			}
		}
	}
}

void destructor_derived_DG_Solver_Face_T (struct Face* face_ptr)
{
	const struct Face*const face            = face_ptr;
	struct DG_Solver_Face_T*const dg_s_face = (struct DG_Solver_Face_T*) face_ptr;

	const int n_neigh = ( face->boundary ? 1 : 2 );
	for (int i = 0; i < n_neigh; ++i) {
		struct Neigh_Info_DG*const ni = &dg_s_face->neigh_info[i];

		destructor_conditional_Multiarray_T(ni->grad_coef_f);
		for (int i = 0; i < DIM; ++i) {
			destructor_conditional_const_Matrix_T(ni->d_g_coef_f__d_s_coef[0][i]);
			destructor_conditional_const_Matrix_T(ni->d_g_coef_f__d_s_coef[1][i]);
		}
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

#include "undef_templates_face_solver.h"
#include "undef_templates_face_solver_dg.h"
#include "undef_templates_volume_solver.h"

#include "undef_templates_matrix.h"
#include "undef_templates_multiarray.h"

#include "undef_templates_test_case.h"
