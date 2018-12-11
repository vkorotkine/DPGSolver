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

#include "test_integration_support_face.h"

#include <assert.h>

#include "macros.h"

#include "multiarray.h"

#include "face_solver.h"

#include "boundary.h"
#include "numerical_flux.h"
#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void constructor_Numerical_Flux_Input_data_with_gradients
	(struct Numerical_Flux_Input*const num_flux_i, const struct Solver_Face*const s_face,
	 const struct Simulation*const sim)
{
	const struct Test_Case*const test_case = (struct Test_Case*)sim->test_case_rc->tc;

	constructor_Boundary_Value_Input_face_fptr
	constructor_Boundary_Value_Input_face_fcl = test_case->constructor_Boundary_Value_Input_face_fcl;
	if (test_case->has_2nd_order &&
	    constructor_Boundary_Value_Input_face_fcl == constructor_Boundary_Value_Input_face_s_fcl_interp)
		constructor_Boundary_Value_Input_face_fcl = constructor_Boundary_Value_Input_face_sg_fcl_interp;

	constructor_Boundary_Value_Input_face_fcl(&num_flux_i->bv_l,s_face,sim);                // keep
	s_face->constructor_Boundary_Value_fcl(&num_flux_i->bv_r,&num_flux_i->bv_l,s_face,sim); // keep
}

void constructor_Numerical_Flux_Input_c_data_members
	(struct Numerical_Flux_Input_c*const num_flux_c_i, struct Numerical_Flux_Input*const num_flux_i, const char side)
{
	assert(side == 'l' || side == 'b'); // Add support.

	num_flux_c_i->bv_l.normals = constructor_copy_const_Multiarray_c_Multiarray_d(num_flux_i->bv_l.normals);
	num_flux_c_i->bv_l.xyz     = constructor_copy_const_Multiarray_c_Multiarray_d(num_flux_i->bv_l.xyz);
	num_flux_c_i->bv_l.s       = constructor_copy_const_Multiarray_c_Multiarray_d(num_flux_i->bv_l.s);
	if (num_flux_i->bv_l.g)
		num_flux_c_i->bv_l.g = constructor_copy_const_Multiarray_c_Multiarray_d(num_flux_i->bv_l.g);

	if (side == 'b') {
		num_flux_c_i->bv_r.s       = constructor_copy_const_Multiarray_c_Multiarray_d(num_flux_i->bv_r.s);
		if (num_flux_i->bv_r.g)
			num_flux_c_i->bv_r.g = constructor_copy_const_Multiarray_c_Multiarray_d(num_flux_i->bv_r.g);
	}
}

void destructor_Numerical_Flux_Input_c_data_members (struct Numerical_Flux_Input_c*const num_flux_c_i, const char side)
{
	assert(side == 'l' || side == 'b'); // Add support.

	destructor_const_Multiarray_c(num_flux_c_i->bv_l.normals);
	destructor_const_Multiarray_c(num_flux_c_i->bv_l.xyz);
	destructor_const_Multiarray_c(num_flux_c_i->bv_l.s);
	destructor_conditional_const_Multiarray_c(num_flux_c_i->bv_l.g);

	if (side == 'b') {
		destructor_conditional_const_Multiarray_c(num_flux_c_i->bv_r.s);
		destructor_conditional_const_Multiarray_c(num_flux_c_i->bv_r.g);
	}
}

void constructor_Boundary_Value_c_data
	(struct Numerical_Flux_Input_c* num_flux_i, const struct Solver_Face_c* s_face, const struct Simulation* sim)
{
	s_face->constructor_Boundary_Value_fcl(&num_flux_i->bv_r,&num_flux_i->bv_l,s_face,sim); // destructed
}

void destructor_Boundary_Value_c_data (struct Numerical_Flux_Input_c* num_flux_i)
{
	destructor_Boundary_Value_c(&num_flux_i->bv_r);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
