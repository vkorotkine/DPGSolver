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
#include "solution.h"

// Static function declarations ************************************************************************************* //

/// \brief Constructor for the members of a \ref Solver_Volume, excluding the base member.
static void construct_Solver_Volume
	(struct Solver_Volume* volume, ///< \ref Volume.
	 const struct Simulation* sim  ///< \ref Simulation.
	);

/// \brief Destructor for a \ref Solver_Volume.
static void destructor_Solver_Volume
	(struct Solver_Volume* volume ///< Standard.
	);

// Interface functions ********************************************************************************************** //

void construct_Solver_Volumes (struct Simulation*const sim)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next)
		construct_Solver_Volume((struct Solver_Volume*)curr,sim);
}

void destructor_Solver_Volumes (struct Intrusive_List* solver_volumes)
{
	for (const struct Intrusive_Link* curr = solver_volumes->first; curr; ) {
		struct Intrusive_Link* next = curr->next;
		destructor_Solver_Volume((struct Solver_Volume*) curr);
		curr = next;
	}
}

static void destructor_Solver_Volume (struct Solver_Volume* volume)
{
	destructor_const_Multiarray_d(volume->geom_coef);
	destructor_Multiarray_d(volume->sol_coef);
	destructor_Multiarray_d(volume->grad_coef);
	destructor_const_Multiarray_d(volume->metrics_vm);
	destructor_const_Multiarray_d(volume->metrics_vc);
	destructor_const_Multiarray_d(volume->jacobian_det_vc);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void construct_Solver_Volume (struct Solver_Volume* volume, const struct Simulation* sim)
{
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
