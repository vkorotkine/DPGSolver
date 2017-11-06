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

#include "solution_advection_default.h"

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

/** \brief Version of \ref constructor_sol_fptr used for the default linear advection test cases.
 *  \return See brief. */
static struct Multiarray_d* constructor_sol_advection_default
	(const struct const_Multiarray_d* xyz, ///< See brief.
	 const struct Simulation* sim          ///< See brief.
	);

// Interface functions ********************************************************************************************** //

void set_sol_advection_default (const struct Simulation* sim, struct Solution_Container sol_cont)
{
	const char ce_type   = sol_cont.ce_type,
	           node_kind = sol_cont.node_kind;

	assert(ce_type == 'v'); // Add support for faces if necessary.

	const struct const_Multiarray_d* xyz = constructor_xyz_v(sim,sol_cont.volume,node_kind); // destructed
	struct Multiarray_d* sol             = constructor_sol_advection_default(xyz,sim);                // destructed
	destructor_const_Multiarray_d(xyz);

	update_Solution_Container_sol(&sol_cont,sol);
	destructor_Multiarray_d(sol);
}

const struct const_Multiarray_d* constructor_const_sol_advection_default
	(const struct const_Multiarray_d* xyz, const struct Simulation* sim)
{
	struct Multiarray_d* sol = constructor_sol_advection_default(xyz,sim); // returned
	return (const struct const_Multiarray_d*) sol;
}

void compute_source_advection_default (const struct Simulation* sim, struct Solver_Volume* volume)
{
	UNUSED(sim);
	UNUSED(volume);
	EXIT_ADD_SUPPORT;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief 1d version of \ref constructor_sol_advection_default.
 *  \return See brief. */
static struct Multiarray_d* constructor_sol_advection_default_1d
	(const struct const_Multiarray_d* xyz, ///< See brief.
	 const struct Simulation* sim          ///< See brief.
	);

static struct Multiarray_d* constructor_sol_advection_default
	(const struct const_Multiarray_d* xyz, const struct Simulation* sim)
{
	static bool parameters_set = false;
	static mutable_constructor_sol_fptr constructor_sol = NULL;

	if (!parameters_set) {
		parameters_set = true;
		if (sim->d == 1)
			constructor_sol = constructor_sol_advection_default_1d;
		else
			EXIT_UNSUPPORTED;
	}
	return constructor_sol(xyz,sim);
}

// Level 1 ********************************************************************************************************** //

static struct Multiarray_d* constructor_sol_advection_default_1d
	(const struct const_Multiarray_d* xyz, const struct Simulation* sim)
{
	assert(sim->d == 1);
	assert(xyz->extents[1] == 1);

	// Compute the solution
	const ptrdiff_t n_vs = xyz->extents[0];
	const int n_var = sim->test_case->n_var;

	struct Multiarray_d* sol = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){n_vs,n_var}); // returned

	const double* x = get_col_const_Multiarray_d(0,xyz);

	double* u = get_col_Multiarray_d(0,sol);
	for (int i = 0; i < n_vs; ++i)
		u[i] = sin(2.15*x[i]+0.23);

	return sol;
}
