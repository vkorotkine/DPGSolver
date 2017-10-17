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

#include "solution.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"
#include "definitions_intrusive.h"

#include "multiarray.h"
#include "matrix.h"

#include "computational_elements.h"
#include "element.h"
#include "element_solution.h"
#include "volume.h"
#include "volume_solver.h"

#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void set_initial_solution (struct Simulation* sim)
{
	const struct Test_Case* test_case = sim->test_case;
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume* volume = (struct Solver_Volume*) curr;

		test_case->set_sol_coef_v(sim,volume);
		test_case->set_grad_coef_v(sim,volume);
	}
}

void set_grad_coef_v_do_nothing (const struct Simulation* sim, struct Solver_Volume* volume)
{
	UNUSED(sim);
	UNUSED(volume);
	return;
}

void set_sol_coef_f_do_nothing (const struct Simulation* sim, struct Solver_Face* face)
{
	UNUSED(sim);
	UNUSED(face);
	return;
}

void set_grad_coef_f_do_nothing (const struct Simulation* sim, struct Solver_Face* face)
{
	UNUSED(sim);
	UNUSED(face);
	return;
}

const struct const_Multiarray_d* constructor_xyz_vs (const struct Simulation* sim, struct Solver_Volume* volume)
{
	assert(sim->elements->name == IL_SOLUTION_ELEMENT);

	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	const char op_format = 'd';

	struct Volume* base_volume = (struct Volume*) volume;
	struct const_Solution_Element* element = (struct const_Solution_Element*) base_volume->element;

	const int d = ((struct const_Element*)element)->d;

	const int p = volume->p_ref;

	const struct Operator* cv0_vg_vs =
		(!base_volume->curved ? get_Multiarray_Operator(element->cv0_vgs_vs,(ptrdiff_t[]){0,0,p,1})
		                      : get_Multiarray_Operator(element->cv0_vgc_vs,(ptrdiff_t[]){0,0,p,p}) );

	const int n_vs = cv0_vg_vs->op_std->ext_0;

	const struct const_Multiarray_d*const geom_coef = volume->geom_coef;
	struct Multiarray_d* xyz_vs = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){n_vs,d}); // returned
	mm_NN1C_Operator_Multiarray_d(cv0_vg_vs,geom_coef,xyz_vs,op_format,geom_coef->order,NULL,NULL);

	return (const struct const_Multiarray_d*) xyz_vs;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
