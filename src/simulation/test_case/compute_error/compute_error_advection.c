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

#include "compute_error_advection.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "macros.h"

#include "element_error.h"
#include "volume.h"
#include "volume_solver.h"

#include "multiarray.h"
#include "vector.h"

#include "compute_error.h"
#include "const_cast.h"
#include "element.h"
#include "intrusive.h"
//#include "multiarray_operator.h"
//#include "operator.h"
#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Return a statically allocated `char*` holding the specific header for all of the Euler variables.
 *  \return See brief. */
static const char* compute_header_spec_advection_all
	(const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Constructor for a derived \ref Solver_Volume used to compute the exact solution.
 *  \return See brief. */
static struct Solver_Volume* constructor_Solver_Volume_exact ();

/// \brief Destructor for a derived \ref Solver_Volume used to compute the exact solution.
static void destructor_Solver_Volume_exact
	(struct Solver_Volume* s_vol_ex ///< Standard.
	);

/// \brief Set the relevant members of a duplicate \ref Solver_Volume used to compute the exact solution.
static void set_Solver_Volume_exact
	(struct Solver_Volume* s_vol_ex, ///< The partially duplicated \ref Solver_Volume.
	 struct Solver_Volume* s_vol     ///< The \ref Solver_Volume.
	);

// Interface functions ********************************************************************************************** //

struct Error_CE* constructor_Error_CE_advection_all (const struct Simulation* sim)
{
	const int n_out = 1;

/// \todo merge the common elements of the error computation functions when finished with advection.
/// \todo bubble up duplicated static functions.
	const double domain_volume = compute_domain_volume(sim);
	const char* header_spec    = compute_header_spec_advection_all(sim);
	struct Vector_d* sol_L2 = constructor_empty_Vector_d(n_out); // moved
	set_to_value_Vector_d(sol_L2,0.0);

	struct Solution_Container sol_cont =
		{ .ce_type = 'v', .cv_type = 'v', .node_kind = 'c', .volume = NULL, .face = NULL, .sol = NULL, };

	int domain_order = -1;

	struct Solver_Volume* s_vol[2] = { NULL, constructor_Solver_Volume_exact(), }; // destructed
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		s_vol[0] = (struct Solver_Volume*) curr;

		set_Solver_Volume_exact(s_vol[1],s_vol[0]);

		struct Multiarray_d* sol[2] = { NULL, NULL, };

		sol[0] = constructor_sol_v(sim,s_vol[0],sol_cont.node_kind);       // destructed
		sol[1] = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){0,0}), // destructed

		sol_cont.sol = sol[1];
		sol_cont.volume = s_vol[1];
		sim->test_case->set_sol(sim,sol_cont);

		subtract_in_place_Multiarray_d(sol[0],(const struct const_Multiarray_d*)sol[1]);

		increment_vol_errors_l2_2(sol_L2,(const struct const_Multiarray_d*)sol[0],s_vol[0]);

		for (int i = 0; i < 2; ++i)
			destructor_Multiarray_d(sol[i]);

		if (domain_order == -1)
			domain_order = s_vol[0]->p_ref;
		else
			assert(domain_order == s_vol[0]->p_ref);

	}
	destructor_Solver_Volume_exact(s_vol[1]);

	for (int i = 0; i < sol_L2->ext_0; ++i)
		sol_L2->data[i] = sqrt(sol_L2->data[i]/domain_volume);

	struct Vector_i* expected_order = constructor_empty_Vector_i(n_out); // moved
	set_to_value_Vector_i(expected_order,domain_order+1);

	struct Error_CE* error_ce = malloc(sizeof *error_ce); // returned

	const_cast_d(&error_ce->domain_volume,domain_volume);
	error_ce->sol_L2         = (const struct const_Vector_d*) sol_L2;
	error_ce->expected_order = (const struct const_Vector_i*) expected_order;
	const_cast_c1(&error_ce->header_spec,header_spec);

	return error_ce;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static const char* compute_header_spec_advection_all (const struct Simulation* sim)
{
	UNUSED(sim);
	static char header_spec[STRLEN_MAX];

	sprintf(header_spec,"%-14s","L2u");

	return header_spec;
}

static struct Solver_Volume* constructor_Solver_Volume_exact ()
{
	struct Solver_Volume* s_vol_ex = calloc(1,sizeof *s_vol_ex); // returned

	s_vol_ex->sol_coef = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){0,0}); // destructed

	return s_vol_ex;
}

static void destructor_Solver_Volume_exact (struct Solver_Volume* s_vol_ex)
{
	destructor_Multiarray_d(s_vol_ex->sol_coef);
	free(s_vol_ex);
}

static void set_Solver_Volume_exact (struct Solver_Volume* s_vol_ex, struct Solver_Volume* s_vol)
{
	struct Volume* b_vol_ex = (struct Volume*) s_vol_ex,
	             * b_vol    = (struct Volume*) s_vol;

	const_cast_b(&b_vol_ex->curved,b_vol->curved);
	const_cast_const_Element(&b_vol_ex->element,b_vol->element);

	const_cast_i(&s_vol_ex->p_ref,s_vol->p_ref);
	const_constructor_move_const_Multiarray_d(&s_vol_ex->geom_coef,s_vol->geom_coef);
}
