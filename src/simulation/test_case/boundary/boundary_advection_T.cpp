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
 *  \brief Provides the templated linear advection boundary condition functions.
 */

#include <assert.h>
#include <math.h>

#include "macros.h"
#include "definitions_mesh.h"


#include "def_templates_boundary.h"
#include "boundary_pde_T.cpp"

#include "def_templates_multiarray.h"

#include "def_templates_face_solver.h"

#include "def_templates_math_functions.h"
#include "def_templates_solution_advection.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void constructor_Boundary_Value_T_advection_upwind
	(struct Boundary_Value_T* bv, const struct Boundary_Value_Input_T* bv_i, const struct Solver_Face_T* face,
	 const struct Simulation* sim)
{
	UNUSED(face);

	static bool need_input = true;
	static struct Sol_Data__Advection_T sol_data;
	if (need_input) {
		need_input = false;
		read_data_advection_T(&sol_data);
	}

	const bool* c_m = bv_i->compute_member;

	assert(c_m[0] == true);
	struct Multiarray_T*const s = (struct Multiarray_T*) constructor_sol_bv(bv_i->xyz,sim); // moved
	const struct const_Multiarray_T*const s_i = constructor_copy_const_Multiarray_T(bv_i->s); // destructed

	const Type*const xyz[DIM] = ARRAY_DIM( get_col_const_Multiarray_T(0,bv_i->xyz),
	                                       get_col_const_Multiarray_T(1,bv_i->xyz),
	                                       get_col_const_Multiarray_T(2,bv_i->xyz) );
	const struct const_Multiarray_T* normals = bv_i->normals;
	assert(normals->layout == 'R');

	const ptrdiff_t n_n = s->extents[0];
	for (int n = 0; n < n_n; n++) {
		const Type xyz_n[DIM] = ARRAY_DIM(xyz[0][n],xyz[1][n],xyz[2][n]);
		const Type*const data_n = get_row_const_Multiarray_T(n,normals);
		const Real*const b_adv  = sol_data.cppompute_b_adv(xyz_n);
		const Real b_dot_n = dot_R_from_RT(DIM,b_adv,data_n);

		if (b_dot_n <= 0.0)
			; // Inflow (do nothing)
		else
			s->data[n] = s_i->data[n];
	}
	destructor_const_Multiarray_T(s_i);
	bv->s = (struct const_Multiarray_T*)s; // keep

	if (c_m[1] == true) {
		const ptrdiff_t n_n  = bv->s->extents[0],
		                n_vr = bv->s->extents[1];
		struct Multiarray_T* ds_ds = constructor_zero_Multiarray_T('C',3,(ptrdiff_t[]){n_n,n_vr,n_vr}); // moved

		for (int n = 0; n < n_n; n++) {
			const Type xyz_n[DIM] = ARRAY_DIM(xyz[0][n],xyz[1][n],xyz[2][n]);
			const Type*const data_n = get_row_const_Multiarray_T(n,normals);
			const Real*const b_adv  = sol_data.cppompute_b_adv(xyz_n);
			const Real b_dot_n = dot_R_from_RT(DIM,b_adv,data_n);

			if (b_dot_n <= 0.0)
				; // Inflow (do nothing)
			else
				ds_ds->data[n] = 1.0;
		}
		bv->ds_ds = (const struct const_Multiarray_T*) ds_ds; // keep
	}
	assert(c_m[2] == false);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

#include "undef_templates_boundary.h"

#include "undef_templates_multiarray.h"

#include "undef_templates_face_solver.h"

#include "undef_templates_math_functions.h"
#include "undef_templates_solution_advection.h"
