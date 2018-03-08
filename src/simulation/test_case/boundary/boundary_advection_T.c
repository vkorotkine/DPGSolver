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

#include "macros.h"
#include "definitions_mesh.h"


#include "def_templates_boundary.h"
#include "boundary_pde_T.c"

#include "def_templates_multiarray.h"

#include "def_templates_face_solver.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void constructor_Boundary_Value_T_advection_inflow
	(struct Boundary_Value_T* bv, const struct Boundary_Value_Input_T* bv_i, const struct Solver_Face_T* face,
	 const struct Simulation* sim)
{
	UNUSED(face);
	const struct const_Multiarray_d*const xyz = bv_i->xyz;
	const bool* c_m = bv_i->compute_member;

	assert(c_m[0] == true);
	bv->s = constructor_sol_bv(xyz,sim); // keep

	if (c_m[1] == true) {
		const ptrdiff_t n_n  = bv->s->extents[0],
		                n_vr = bv->s->extents[1];
		struct Multiarray_T* ds_ds = constructor_empty_Multiarray_T('C',3,(ptrdiff_t[]){n_n,n_vr,n_vr}); // moved
		set_to_value_Multiarray_T(ds_ds,0.0);
		bv->ds_ds = (const struct const_Multiarray_T*) ds_ds; // keep
	}
	assert(c_m[2] == false);
}

void constructor_Boundary_Value_T_advection_outflow
	(struct Boundary_Value_T* bv, const struct Boundary_Value_Input_T* bv_i, const struct Solver_Face_T* face,
	 const struct Simulation* sim)
{
	UNUSED(face);
	UNUSED(sim);

	const bool* c_m = bv_i->compute_member;

	assert(c_m[0] == true);
	bv->s = constructor_copy_const_Multiarray_T(bv_i->s); // keep

	if (c_m[1] == true) {
		const ptrdiff_t n_n  = bv->s->extents[0],
		                n_vr = bv->s->extents[1];
		struct Multiarray_T* ds_ds = constructor_empty_Multiarray_T('C',3,(ptrdiff_t[]){n_n,n_vr,n_vr}); // moved
		set_to_value_Multiarray_T(ds_ds,1.0);
		bv->ds_ds = (const struct const_Multiarray_T*) ds_ds;
	}
	assert(c_m[2] == false);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
