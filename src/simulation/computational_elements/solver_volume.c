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

#include "solver_volume.h"

#include <string.h>

#include "macros.h"
#include "definitions_intrusive.h"

#include "multiarray.h"
#include "vector.h"

#include "const_cast.h"
#include "geometry.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void constructor_derived_Solver_Volume (struct Volume* volume_ptr, const struct Simulation* sim)
{
	struct Solver_Volume* volume = (struct Solver_Volume*) volume_ptr;

	const_cast_i(&volume->p_ref,sim->p_s_v[0]);
	const_constructor_move_Multiarray_d(&volume->geom_coef,constructor_default_Multiarray_d());

	volume->sol_coef  = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){0,0});   // destructed
	volume->grad_coef = constructor_empty_Multiarray_d('C',3,(ptrdiff_t[]){0,0,0}); // destructed

	const_constructor_move_Multiarray_d(
		&volume->metrics_vm,constructor_empty_Multiarray_d('C',3,(ptrdiff_t[]){0,0,0}));  // destructed
	const_constructor_move_Multiarray_d(
		&volume->metrics_vc,constructor_empty_Multiarray_d('C',3,(ptrdiff_t[]){0,0,0}));  // destructed
	const_constructor_move_Multiarray_d(
		&volume->jacobian_det_vc,constructor_empty_Multiarray_d('C',1,(ptrdiff_t[]){0})); // destructed
}

void destructor_derived_Solver_Volume (struct Volume* volume_ptr)
{
	struct Solver_Volume* volume = (struct Solver_Volume*) volume_ptr;

	destructor_const_Multiarray_d(volume->geom_coef);
	destructor_Multiarray_d(volume->sol_coef);
	destructor_Multiarray_d(volume->grad_coef);
	destructor_const_Multiarray_d(volume->metrics_vm);
	destructor_const_Multiarray_d(volume->metrics_vc);
	destructor_const_Multiarray_d(volume->jacobian_det_vc);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
