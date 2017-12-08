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

#include <string.h>

#include "macros.h"
#include "definitions_intrusive.h"


#include "def_templates_volume_solver.h"

#include "def_templates_multiarray.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void constructor_derived_Solver_Volume_T (struct Volume* volume_ptr, const struct Simulation* sim)
{
	struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) volume_ptr;

	const_cast_ptrdiff(&s_vol->ind_dof,-1);
	const_cast_i(&s_vol->p_ref,sim->p_ref[0]);
	const_cast_i(&s_vol->ml,0);
	const_constructor_move_Multiarray_R(&s_vol->geom_coef,constructor_default_Multiarray_R());

	s_vol->sol_coef  = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){0,0});   // destructed
	s_vol->grad_coef = constructor_empty_Multiarray_T('C',3,(ptrdiff_t[]){0,0,0}); // destructed

	const_constructor_move_Multiarray_R(
		&s_vol->metrics_vm,constructor_empty_Multiarray_R('C',3,(ptrdiff_t[]){0,0,0}));  // destructed
	const_constructor_move_Multiarray_R(
		&s_vol->metrics_vc,constructor_empty_Multiarray_R('C',3,(ptrdiff_t[]){0,0,0}));  // destructed
	const_constructor_move_Multiarray_R(
		&s_vol->jacobian_det_vc,constructor_empty_Multiarray_R('C',1,(ptrdiff_t[]){0})); // destructed
}

void destructor_derived_Solver_Volume_T (struct Volume* volume_ptr)
{
	struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) volume_ptr;

	destructor_const_Multiarray_R(s_vol->geom_coef);
	destructor_Multiarray_T(s_vol->sol_coef);
	destructor_Multiarray_T(s_vol->grad_coef);
	destructor_const_Multiarray_R(s_vol->metrics_vm);
	destructor_const_Multiarray_R(s_vol->metrics_vc);
	destructor_const_Multiarray_R(s_vol->jacobian_det_vc);
}

const struct const_Vector_d* get_operator__w_vc__s_e_T (const struct Solver_Volume_T* s_vol)
{
	struct Volume* vol               = (struct Volume*) s_vol;
	const struct Solver_Element* s_e = (struct Solver_Element*) vol->element;

	const int p      = s_vol->p_ref,
	          curved = vol->curved;
	return get_const_Multiarray_Vector_d(s_e->w_vc[curved],(ptrdiff_t[]){0,0,p,p});
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
