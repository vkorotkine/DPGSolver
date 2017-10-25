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

#include "compute_error_euler.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"

#include "element_error.h"
#include "volume.h"
#include "volume_solver.h"

#include "multiarray.h"
#include "vector.h"

#include "compute_error.h"
#include "const_cast.h"
#include "const_cast.h"
#include "element.h"
#include "intrusive.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "solution_euler.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Constructor for a derived \ref Solver_Volume used to compute the exact solution.
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

struct Error_CE* constructor_Error_CE_euler_all (const struct Simulation* sim)
{
	const int d     = sim->d;
	const int n_out = d+2+1;

	double domain_volume = 0.0;
	struct Vector_d* sol_L2 = constructor_empty_Vector_d(n_out); // moved

	struct Solution_Container sol_cont =
		{ .ce_type = 'v', .cv_type = 'v', .node_kind = 'c', .volume = NULL, .face = NULL, .sol = NULL, };

	ptrdiff_t s_extents[2] = { 0, 1 };
	struct Multiarray_d* s = constructor_move_Multiarray_d_dyn_extents('C',2,s_extents,false,NULL); // destructed

	struct Solver_Volume* s_vol[2] = { NULL, constructor_Solver_Volume_exact(), }; // destructed
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		s_vol[0] = (struct Solver_Volume*) curr;
		domain_volume += compute_volume(s_vol[0]);

		set_Solver_Volume_exact(s_vol[1],s_vol[0]);

printf("\n\n\n");
		s->extents[0] = s_vol[0]->jacobian_det_vc->extents[0];
		struct Multiarray_d* sol[2] = { NULL, NULL, };
		for (int i = 0; i < 2; ++i) {
			sol[i] = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){0,0}), // destructed

			sol_cont.sol = sol[i];
			sol_cont.volume = s_vol[i];
			sim->test_case->set_sol(sim,sol_cont);

			sol[i]->extents[1] = d+2;

print_Multiarray_d(sol[i]);
			convert_variables(sol[i],'c','p');

print_Multiarray_d(sol[i]);
			s->data = get_col_Multiarray_d(d+2,sol[i]);
			compute_entropy(s,(const struct const_Multiarray_d*)sol[i],'p');

			sol[i]->extents[1]    = d+2+1;
print_Multiarray_d(sol[i]);
		}

		for (int i = 0; i < 2; ++i)
			destructor_Multiarray_d(sol[i]);
EXIT_UNSUPPORTED;
	}
	destructor_Solver_Volume_exact(s_vol[1]);

	struct Error_CE* error_ce = malloc(sizeof *error_ce); // returned

	const_cast_d(&error_ce->domain_volume,domain_volume);
	error_ce->sol_L2 = (const struct const_Vector_d*) sol_L2;

EXIT_UNSUPPORTED;
	return error_ce;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct Solver_Volume* constructor_Solver_Volume_exact ()
{
	struct Solver_Volume* s_vol_ex = calloc(1,sizeof *s_vol_ex); // returned

	s_vol_ex->sol_coef = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){0,0}); // destructed

	return s_vol_ex;
}

static void destructor_Solver_Volume_exact (struct Solver_Volume* s_vol_ex)
{
	destructor_Multiarray_d(s_vol_ex->sol_coef);
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
