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
 *  \brief Provides the templated diffusion boundary condition functions.
 */

#include <assert.h>

#include "macros.h"
#include "definitions_mesh.h"


#include "def_templates_boundary.h"
#include "boundary_pde_T.c"

#include "def_templates_multiarray.h"

#include "def_templates_face_solver.h"

// Static function declarations ************************************************************************************* //

/** \brief Constructor for the Multiarray_\* holding the solution gradient at the given xyz coordinates.
 *  \return See brief. */
static const struct const_Multiarray_T* constructor_grad_bv
	(const struct const_Multiarray_d* xyz, ///< xyz coordinates.
	 const struct Simulation* sim          ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void constructor_Boundary_Value_T_diffusion_dirichlet
	(struct Boundary_Value_T* bv, const struct Boundary_Value_Input_T* bv_i, const struct Solver_Face_T* face,
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
		 *
		 * Move this correction out of the boundary condition functions and into the bv_i->xyz specification.
		 */
		EXIT_ADD_SUPPORT;
	}
	const bool* c_m = bv_i->compute_member;

	assert(c_m[0] == true);
	struct Multiarray_T*const s_ex = (struct Multiarray_T*)constructor_sol_bv(xyz,sim); // keep

	const ptrdiff_t size_s = compute_size(s_ex->order,s_ex->extents);
	for (int i = 0; i < size_s; ++i) {
		s_ex->data[i] *= 2.0;
		s_ex->data[i] -= bv_i->s->data[i];
	}
	bv->s = (struct const_Multiarray_T*) s_ex;

	if (c_m[1] == true) {
		const ptrdiff_t n_n  = bv->s->extents[0],
		                n_vr = bv->s->extents[1];
		struct Multiarray_T* ds_ds = constructor_empty_Multiarray_T('C',3,(ptrdiff_t[]){n_n,n_vr,n_vr}); // moved
		set_to_value_Multiarray_T(ds_ds,-1.0);
		bv->ds_ds = (const struct const_Multiarray_T*) ds_ds; // keep
	}

	assert(c_m[2] == true);
	bv->g = constructor_copy_const_Multiarray_T(bv_i->g); // keep

	assert(c_m[3] == false);
	assert(c_m[4] == false);
	assert(c_m[5] == false);
}

void constructor_Boundary_Value_T_diffusion_neumann
	(struct Boundary_Value_T* bv, const struct Boundary_Value_Input_T* bv_i, const struct Solver_Face_T* face,
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
		 *
		 * Move this correction out of the boundary condition functions and into the bv_i->xyz specification.
		 */
		EXIT_ADD_SUPPORT;
	}
	const bool* c_m = bv_i->compute_member;

	assert(c_m[0] == true);
	bv->s = constructor_copy_const_Multiarray_T(bv_i->s); // keep

	if (c_m[1] == true) {
		const ptrdiff_t n_n  = bv->s->extents[0],
		                n_vr = bv->s->extents[1];
		struct Multiarray_T* ds_ds = constructor_empty_Multiarray_T('C',3,(ptrdiff_t[]){n_n,n_vr,n_vr}); // moved
		set_to_value_Multiarray_T(ds_ds,1.0);
		bv->ds_ds = (const struct const_Multiarray_T*) ds_ds; // keep
	}

	assert(c_m[2] == true);
	struct Multiarray_T*const g_ex = (struct Multiarray_T*)constructor_grad_bv(xyz,sim); // keep

	const ptrdiff_t size_g = compute_size(g_ex->order,g_ex->extents);
	for (int i = 0; i < size_g; ++i) {
		g_ex->data[i] *= 2.0;
		g_ex->data[i] -= bv_i->g->data[i];
	}
	bv->g = (struct const_Multiarray_T*) g_ex;

	assert(c_m[3] == false);
	assert(c_m[4] == false);
	assert(c_m[5] == false);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static const struct const_Multiarray_T* constructor_grad_bv
	(const struct const_Multiarray_d* xyz, const struct Simulation* sim)
{
	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	return test_case->constructor_grad(xyz,sim); // returned
}
