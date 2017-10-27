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
	const char ce_type   = sol_cont.ce_type,
	           cv_type   = sol_cont.cv_type,
	           node_kind = sol_cont.node_kind;

	assert(ce_type == 'v'); // Add support for faces if necessary.

	const struct const_Multiarray_d* xyz = constructor_xyz_v(sim,sol_cont.volume,node_kind); // destructed
	struct Multiarray_d* sol             = constructor_sol_peterson(sim,xyz);                // destructed
	destructor_const_Multiarray_d(xyz);

/// \todo function for this (also change in periodic_vortex).
	assert((cv_type == 'c') || (cv_type == 'v'));
	if (cv_type == 'v') {
		assert(sol_cont.sol->data != NULL);

		sol_cont.sol->extents[0] = sol->extents[0];
		sol_cont.sol->extents[1] = sol->extents[1];
		free(sol_cont.sol->data);
		sol_cont.sol->data = sol->data;

		sol->owns_data = false;
	} else if (cv_type == 'c') {
		assert(node_kind == 's');
		compute_coef_from_val_vs(sol_cont.volume,(struct const_Multiarray_d*)sol,sol_cont.sol);
	}

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

/// \brief Container for solution data relating to "pe"terson solution.
struct Sol_Data__pe {
	struct Sol_Data__Advection sd_adv; ///< \ref Sol_Data__Advection.
};

/** \brief Return the statically allocated \ref Sol_Data__pe container.
 *  \return See brief. */
static struct Sol_Data__pe get_sol_data
	(const struct Simulation* sim ///< \ref Simulation.
	);

static struct Multiarray_d* constructor_sol_peterson
	(const struct Simulation* sim, const struct const_Multiarray_d* xyz)
{
	assert(sim->d == 2);

	const struct Sol_Data__pe sol_data = get_sol_data(sim);

	// Compute the solution
	const ptrdiff_t n_vs = xyz->extents[0];
	const int n_var = sim->test_case->n_var;

	struct Multiarray_d* sol = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){n_vs,n_var}); // returned

	const double* b_adv = ((struct Sol_Data__Advection*)&sol_data)->b_adv;
	assert((b_adv[0] == 0.0) && (b_adv[1] == 1.0)); /* Can be made flexible in future but solution below must be
	                                                 * modified. */

	const double* x = get_col_const_Multiarray_d(0,xyz);

	double* u = get_col_Multiarray_d(0,sol);
	for (int i = 0; i < n_vs; ++i)
		u[i] = sin(2.0*x[i]);

	return sol;
}

// Level 1 ********************************************************************************************************** //

/// \brief Read the required solution data into \ref Sol_Data__pe.
static void read_data_peterson
	(const char*const input_path,       ///< Defined in \ref fopen_input.
	 struct Sol_Data__pe*const sol_data ///< \ref Sol_Data__pe.
	);

/// \brief Set the remaining required solution data of \ref Sol_Data__pe based on the read values.
static void set_data_peterson
	(struct Sol_Data__pe*const sol_data ///< \ref Sol_Data__pe.
	);

static struct Sol_Data__pe get_sol_data (const struct Simulation* sim)
{
	static bool need_input = true;

	static struct Sol_Data__pe sol_data;
	if (need_input) {
		need_input = false;
		read_data_peterson(sim->input_path,&sol_data);
		set_data_peterson(&sol_data);
	}

	return sol_data;
}

// Level 2 ********************************************************************************************************** //

static void read_data_peterson (const char*const input_path, struct Sol_Data__pe*const sol_data)
{
	read_data_advection(input_path,&sol_data->sd_adv);
}

static void set_data_peterson (struct Sol_Data__pe*const sol_data)
{
	UNUSED(sol_data);
	return;
}
