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

#include "boundary_advection.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"
#include "definitions_mesh.h"

#include "multiarray.h"

#include "boundary.h"
#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void constructor_Boundary_Value_advection_inflow
	(struct Boundary_Value* bv, const struct Boundary_Value_Input* bv_i, const struct Solver_Face* face,
	 const struct Simulation* sim)
{
	UNUSED(face);

	const struct const_Multiarray_d* xyz = NULL;
	if (sim->domain_type == DOM_STRAIGHT) {
		xyz = bv_i->xyz;
	} else {
		/* Boundary values computed by evaluating the exact solution must be applied on the exact boundary which
		 * requires a correction of the boundary node coordinates. Otherwise, the approximate domain is being
		 * implicitly treated as the exact domain.
		 *
		 * Can correct the coordinates in the direction of the normal vector.
		 */
		EXIT_ADD_SUPPORT;
	}

	const bool* c_m = bv_i->compute_member;

	assert(c_m[0] == true);
	bv->s = sim->test_case->constructor_sol(xyz,sim);

	if (c_m[1] == true) {
		const ptrdiff_t n_n  = bv->s->extents[0],
		                n_vr = bv->s->extents[1];
		struct Multiarray_d* ds_ds = constructor_empty_Multiarray_d('C',3,(ptrdiff_t[]){n_n,n_vr,n_vr}); // moved
		set_to_value_Multiarray_d(ds_ds,0.0);
		bv->ds_ds = (const struct const_Multiarray_d*) ds_ds;
	}

	assert(c_m[2] == false);
}

void constructor_Boundary_Value_advection_outflow
	(struct Boundary_Value* bv, const struct Boundary_Value_Input* bv_i, const struct Solver_Face* face,
	 const struct Simulation* sim)
{
	UNUSED(face);
	UNUSED(sim);

	const bool* c_m = bv_i->compute_member;

	assert(c_m[0] == true);
	bv->s = constructor_copy_const_Multiarray_d(bv_i->s);

	if (c_m[1] == true) {
		const ptrdiff_t n_n  = bv->s->extents[0],
		                n_vr = bv->s->extents[1];
		struct Multiarray_d* ds_ds = constructor_empty_Multiarray_d('C',3,(ptrdiff_t[]){n_n,n_vr,n_vr}); // moved
		set_to_value_Multiarray_d(ds_ds,1.0);
		bv->ds_ds = (const struct const_Multiarray_d*) ds_ds;
	}

	assert(c_m[2] == false);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
