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

#include "solution_peterson.h"

#include <assert.h>
#include <math.h>
#include <string.h>

#include "macros.h"
#include "definitions_tol.h"

#include "multiarray.h"

#include "file_processing.h"
#include "simulation.h"
#include "solution.h"
#include "solution_advection.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Return a \ref Multiarray_d\* container holding the solution values at the input coordinates.
 *  \return See brief. */
static struct Multiarray_d* constructor_sol_peterson
	(const struct Simulation* sim,        ///< Defined for \ref set_sol_peterson.
	 const struct const_Multiarray_d* xyz ///< xyz coordinates at which to evaluate the solution.
	);

// Interface functions ********************************************************************************************** //

void set_sol_peterson (const struct Simulation* sim, struct Solution_Container sol_cont)
{
	const struct const_Multiarray_d* xyz = constructor_xyz_sol(sim,&sol_cont); // destructed
	struct Multiarray_d* sol = constructor_sol_peterson(sim,xyz); // destructed
	destructor_const_Multiarray_d(xyz);

	update_Solution_Container_sol(&sol_cont,sol);
	destructor_Multiarray_d(sol);
}

const struct const_Multiarray_d* constructor_const_sol_peterson
	(const struct const_Multiarray_d* xyz, const struct Simulation* sim)
{
	struct Multiarray_d* sol = constructor_sol_peterson(sim,xyz); // returned
	return (const struct const_Multiarray_d*) sol;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct Multiarray_d* constructor_sol_peterson
	(const struct Simulation* sim, const struct const_Multiarray_d* xyz)
{
	assert(sim->d == 2);

	const struct Sol_Data__Advection sol_data = get_sol_data_advection(sim);

	// Compute the solution
	const ptrdiff_t n_vs = xyz->extents[0];
	const int n_var = sim->test_case->n_var;

	struct Multiarray_d* sol = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){n_vs,n_var}); // returned

	const double* b_adv = sol_data.b_adv;
	assert((b_adv[0] == 0.0) && (b_adv[1] == 1.0)); /* Can be made flexible in future but solution below must be
	                                                 * modified. */

	const double* x = get_col_const_Multiarray_d(0,xyz);

	double* u = get_col_Multiarray_d(0,sol);
	for (int i = 0; i < n_vs; ++i)
		u[i] = sin(2.15*x[i]+0.23);

	return sol;
}
