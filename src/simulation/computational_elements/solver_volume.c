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

#include "simulation.h"
#include "geometry.h"
#include "solution.h"

// Static function declarations ************************************************************************************* //

/** \brief Constructor for an individual \ref Solver_Volume.
 *  \return Standard. */
static struct Solver_Volume* constructor_Solver_Volume
	(struct Volume* volume,       ///< \ref Volume.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for an individual \ref Solver_Volume.
static void destructor_Solver_Volume
	(struct Solver_Volume* volume ///< Standard.
	);

// Interface functions ********************************************************************************************** //

struct Intrusive_List* constructor_Solver_Volumes (struct Simulation*const sim)
{
	struct Intrusive_List* volumes        = sim->volumes;
	struct Intrusive_List* solver_volumes = constructor_empty_IL(IL_SOLVER_VOLUME);

	for (struct Intrusive_Link* curr = volumes->first; curr; curr = curr->next)
		push_back_IL(solver_volumes,(struct Intrusive_Link*) constructor_Solver_Volume((struct Volume*)curr,sim));

	destructor_IL(volumes);
	return solver_volumes;
}

void destructor_Solver_Volumes (struct Intrusive_List* solver_volumes)
{
	for (const struct Intrusive_Link* curr = solver_volumes->first; curr; ) {
		struct Intrusive_Link* next = curr->next;
		destructor_Solver_Volume((struct Solver_Volume*) curr);
		curr = next;
	}
	destructor_IL(solver_volumes);
}

static void destructor_Solver_Volume (struct Solver_Volume* volume)
{
	destructor_Multiarray_d(volume->sol_coef);
	destructor_Multiarray_d(volume->grad_coef);
	destructor_const_Multiarray_d(volume->metrics_vg);
	destructor_const_Multiarray_d(volume->metrics_vc);
	destructor_const_Multiarray_d(volume->jacobian_det_vc);

	destructor_Volume((struct Volume*) volume);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct Solver_Volume* constructor_Solver_Volume (struct Volume* volume, const struct Simulation* sim)
{
	struct Solver_Volume* solver_volume = calloc(1,sizeof *solver_volume); // returned

	memcpy(&solver_volume->volume,volume,sizeof *volume); // shallow copy of the base.

	solver_volume->p = sim->p_s_v[0];

	solver_volume->sol_coef  = constructor_default_Multiarray_d(); // destructed
	solver_volume->grad_coef = constructor_default_Multiarray_d(); // destructed

	const_constructor_move_Multiarray_d(&solver_volume->metrics_vg,
	                                    constructor_default_Multiarray_d()); // destructed
	const_constructor_move_Multiarray_d(&solver_volume->metrics_vc,
	                                    constructor_default_Multiarray_d()); // destructed
	const_constructor_move_Multiarray_d(&solver_volume->jacobian_det_vc,
	                                    constructor_default_Multiarray_d()); // destructed

	return solver_volume;
}
