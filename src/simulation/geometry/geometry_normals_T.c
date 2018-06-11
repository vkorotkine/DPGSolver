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
#include <math.h>
#include <string.h>

#include "macros.h"
//#include "definitions_alloc.h"
#include "definitions_bc.h"
#include "definitions_core.h"
#include "definitions_math.h"
#include "definitions_tol.h"


#include "def_templates_geometry_normals.h"

#include "def_templates_face_solver.h"

#include "def_templates_multiarray.h"

#include "def_templates_math_functions.h"

// Static function declarations ************************************************************************************* //

/// \brief Set \ref Solver_Face_T::normals_fc to \ref Solver_Face_T::normals_fc_exact if necessary.
static void correct_normals_fc
	(struct Solver_Face_T*const s_face ///< Standard.
	);

// Interface functions ********************************************************************************************** //

correct_for_exact_normal_fptr_T set_correct_for_exact_normal_fptr_T (const struct Simulation*const sim)
{
	if (strstr(sim->pde_spec,"supersonic_vortex") ||
	    strstr(sim->pde_spec,"steady/vortex")) {
		return correct_for_exact_normal_cylinder_T;
	} else if (strstr(sim->pde_spec,"gaussian_bump")) {
		return correct_for_exact_normal_gaussian_bump_T;
	} else {
		EXIT_ERROR("Unsupported: %s\n",sim->pde_spec);
	}
}

void correct_for_exact_normal_cylinder_T (struct Solver_Face_T*const s_face)
{
	assert(DIM >= 2);
	enum { d_max = 2, };

	const struct Face*const face = (struct Face*) s_face;
	if (!face->boundary || face->bc < BC_CURVED_START)
		return;

	const struct const_Multiarray_T*const xyz_fc     = s_face->xyz_fc;
	const struct const_Multiarray_T*const normals_fc = s_face->normals_fc;

	destructor_const_Multiarray_T(s_face->normals_fc_exact);
	struct Multiarray_R*const normals_fc_exact = constructor_copy_Multiarray_R((struct Multiarray_R*)normals_fc); // k
	const_constructor_move_Multiarray_R(&s_face->normals_fc_exact,normals_fc_exact);

	const Real*const xyz[] = ARRAY_DIM( get_col_const_Multiarray_R(0,xyz_fc),
	                                    get_col_const_Multiarray_R(1,xyz_fc),
	                                    get_col_const_Multiarray_R(2,xyz_fc) );

	const ptrdiff_t n_n = xyz_fc->extents[0];
	for (int n = 0; n < n_n; ++n) {
		const Real xyz_n[] = ARRAY_DIM( xyz[0][n], xyz[1][n], xyz[2][n] );
		const Real th = atan2(xyz_n[1],xyz_n[0]);
		const Real n_ex[] = { -cos(th), -sin(th), };
		Real*const n_fc = get_row_Multiarray_R(n,normals_fc_exact);

		const Real n_dot_n_ex = dot_R(d_max,n_ex,n_fc);

		for (int d = 0; d < d_max; ++d)
			n_fc[d] = n_ex[d];
		if (n_dot_n_ex < 0.0) {
			for (int d = 0; d < d_max; ++d)
				n_fc[d] *= -1.0;
		}

		if (DIM > 2) {
			for (int d = 2; d < DIM; ++d)
				assert(equal_R(n_fc[d],0.0,EPS));
		}
	}
	correct_normals_fc(s_face);
}

void correct_for_exact_normal_gaussian_bump_T (struct Solver_Face_T*const s_face)
{
	assert(DIM == 2); // Add support if needed.

	const struct Face*const face = (struct Face*) s_face;
	if (!face->boundary || face->bc < BC_CURVED_START)
		return;

	const struct const_Multiarray_T*const xyz_fc     = s_face->xyz_fc;
	const struct const_Multiarray_T*const normals_fc = s_face->normals_fc;

	destructor_const_Multiarray_T(s_face->normals_fc_exact);
	struct Multiarray_R*const normals_fc_exact = constructor_copy_Multiarray_R((struct Multiarray_R*)normals_fc); // k
	const_constructor_move_Multiarray_R(&s_face->normals_fc_exact,normals_fc_exact);

	const Real*const xyz[] = ARRAY_DIM( get_col_const_Multiarray_R(0,xyz_fc),
	                                    get_col_const_Multiarray_R(1,xyz_fc),
	                                    get_col_const_Multiarray_R(2,xyz_fc) );

	const ptrdiff_t n_n = xyz_fc->extents[0];
	for (int n = 0; n < n_n; ++n) {
		const Real xyz_n[] = ARRAY_DIM( xyz[0][n], xyz[1][n], xyz[2][n] );

		struct Function_Data_GP f_data = { .scale = 1.0, };
		const Real df_dx = f_gaussian_bump(xyz_n[0],1,&f_data);

		Real n_ex[] = { df_dx, -1.0, 0.0, };
		const Real norm_n = sqrt(n_ex[0]*n_ex[0]+n_ex[1]*n_ex[1]);
		for (int i = 0; i < 2; ++i)
			n_ex[i] /= norm_n;

		Real*const n_fc = get_row_Multiarray_R(n,normals_fc_exact);

		const Real n_dot_n_ex = dot_R(DIM,n_ex,n_fc);
		assert(n_dot_n_ex > 0.0);

		for (int d = 0; d < DIM; ++d)
			n_fc[d] = n_ex[d];
	}
	correct_normals_fc(s_face);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void correct_normals_fc (struct Solver_Face_T*const s_face)
{
	if (!using_exact_normals())
		return;
	set_Multiarray_R((struct Multiarray_R*)s_face->normals_fc,s_face->normals_fc_exact);
}
