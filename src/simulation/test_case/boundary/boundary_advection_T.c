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

#include "boundary.h"

#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

#include "boundary_T.c"

// Interface functions ********************************************************************************************** //

/// \brief Version of \ref constructor_Boundary_Value_fptr computing members using the inflow (boundary) values.
void constructor_Boundary_Value_T_advection_inflow
	(struct Boundary_Value_T* bv,               ///< See brief.
	 const struct Boundary_Value_Input_T* bv_i, ///< See brief.
	 const struct Solver_Face* face,            ///< See brief.
	 const struct Simulation* sim               ///< See brief.
	)
{
	UNUSED(face);
	struct Boundary_Value_Input_R* bv_i_r = (struct Boundary_Value_Input_R*) bv_i;

	const struct const_Multiarray_d* xyz = NULL;
	if (sim->domain_type == DOM_STRAIGHT) {
		xyz = bv_i_r->xyz;
	} else {
		/* Boundary values computed by evaluating the exact solution must be applied on the exact boundary which
		 * requires a correction of the boundary node coordinates. Otherwise, the approximate domain is being
		 * implicitly treated as the exact domain.
		 *
		 * Can correct the coordinates in the direction of the normal vector.
		 */
		EXIT_ADD_SUPPORT;
	}
	const bool* c_m = get_compute_member(bv_i);

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

/// \brief Version of \ref constructor_Boundary_Value_fptr computing members using the outflow (extrapolated) values.
void constructor_Boundary_Value_T_advection_outflow
	(struct Boundary_Value_T* bv,               ///< See brief.
	 const struct Boundary_Value_Input_T* bv_i, ///< See brief.
	 const struct Solver_Face* face,            ///< See brief.
	 const struct Simulation* sim               ///< See brief.
	)
{
	UNUSED(face);
	UNUSED(sim);

	const bool* c_m = get_compute_member(bv_i);

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
