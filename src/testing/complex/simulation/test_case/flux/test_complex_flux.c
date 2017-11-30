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

#include "test_complex_flux.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "macros.h"
#include "definitions_core.h"

#include "test_complex_flux_advection.h"
#include "test_complex_flux_euler.h"
#include "test_complex_test_case.h"

#include "complex_multiarray.h"

#include "const_cast.h"
#include "flux.h"
#include "flux_advection.h"
#include "flux_euler.h"
#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Set the function pointers in \ref Flux_Input_c.
static void set_derived_Flux_Input_fptrs
	(struct Flux_Input_c* flux_i ///< \ref Flux_Input_c.
	);

// Interface functions ********************************************************************************************** //

struct Flux_Input_c* constructor_Flux_Input_c (const struct Simulation* sim)
{
	assert(sim->test_case->solver_method_curr == 'i');

	struct Flux_Input* flux_i_b = constructor_Flux_Input(sim); // destructed.

	struct Flux_Input_c* flux_i = calloc(1,sizeof *flux_i); // free.

	memcpy(flux_i,flux_i_b,sizeof(struct Flux_Input)); // shallow copy of the base.
	destructor_Flux_Input(flux_i_b);

	const_cast_b(&flux_i->has_complex_J,has_complex_Jacobians(sim->method));
	set_derived_Flux_Input_fptrs(flux_i);

	return flux_i;
}

void destructor_Flux_Input_c (struct Flux_Input_c* flux_i)
{
	free(flux_i);
}

struct Flux_c* constructor_Flux_c (const struct Flux_Input_c* flux_i)
{
	assert((flux_i->s != NULL && flux_i->s->layout == 'C') ||
	       (flux_i->g != NULL && flux_i->g->layout == 'C'));

	struct Flux_Input* flux_i_b = (struct Flux_Input*) flux_i;
	const bool* compute_member = flux_i_b->compute_member;
	const bool has_c_J = flux_i->has_complex_J;

	const int n_eq = flux_i_b->n_eq,
	          n_vr = flux_i_b->n_var;
	const ptrdiff_t n_n = ( flux_i->s != NULL ? flux_i->s->extents[0] : flux_i->g->extents[0] );

	struct mutable_Flux_c* flux = calloc(1,sizeof *flux); // destructed

	flux->f     = (compute_member[0] ?
		constructor_zero_Multiarray_c('C',3,(ptrdiff_t[]){n_n,DIM,n_eq})          : NULL); // destructed
	flux->df_ds = ((compute_member[1] && has_c_J) ?
		constructor_zero_Multiarray_c('C',4,(ptrdiff_t[]){n_n,DIM,n_eq,n_vr})     : NULL); // destructed
	flux->df_dg = ((compute_member[2] && has_c_J) ?
		constructor_zero_Multiarray_c('C',5,(ptrdiff_t[]){n_n,DIM,n_eq,n_vr,DIM}) : NULL); // destructed

	flux_i->compute_Flux(flux_i,flux);

	return (struct Flux_c*) flux;
}

void destructor_Flux_c (struct Flux_c* flux)
{
	if (flux->f != NULL)
		destructor_const_Multiarray_c(flux->f);
	if (flux->df_ds != NULL)
		destructor_const_Multiarray_c(flux->df_ds);
	if (flux->df_dg != NULL)
		destructor_const_Multiarray_c(flux->df_dg);
	free(flux);
}

void compute_Flux_c_1 (const struct Flux_Input_c* flux_i, struct mutable_Flux_c* flux)
{
	flux_i->compute_Flux_1st(flux_i,flux);
}

void compute_Flux_c_12 (const struct Flux_Input_c* flux_i, struct mutable_Flux_c* flux)
{
	flux_i->compute_Flux_1st(flux_i,flux);
	flux_i->compute_Flux_2nd(flux_i,flux);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void set_derived_Flux_Input_fptrs (struct Flux_Input_c* flux_i)
{
	struct Flux_Input* flux_i_b = (struct Flux_Input*) flux_i;

	// compute_Flux
	if (flux_i_b->compute_Flux == compute_Flux_1)
		flux_i->compute_Flux = compute_Flux_c_1;
	else if (flux_i_b->compute_Flux == compute_Flux_12)
		flux_i->compute_Flux = compute_Flux_c_12;
	else
		EXIT_UNSUPPORTED;

	// compute_Flux_1st
	if (!flux_i->has_complex_J) {
		if (flux_i_b->compute_Flux_1st == compute_Flux_advection_jacobian)
			flux_i->compute_Flux_1st = compute_Flux_c_advection;
		else if (flux_i_b->compute_Flux_1st == compute_Flux_euler_jacobian)
			flux_i->compute_Flux_1st = compute_Flux_c_euler;
		else if (flux_i_b->compute_Flux_1st == NULL)
			flux_i->compute_Flux_1st = NULL;
		else
			EXIT_UNSUPPORTED;
	} else {
		if (flux_i_b->compute_Flux_1st == compute_Flux_advection_jacobian)
			flux_i->compute_Flux_1st = compute_Flux_c_advection_jacobian;
		else if (flux_i_b->compute_Flux_1st == compute_Flux_euler_jacobian)
			flux_i->compute_Flux_1st = compute_Flux_c_euler_jacobian;
		else if (flux_i_b->compute_Flux_1st == NULL)
			flux_i->compute_Flux_1st = NULL;
		else
			EXIT_UNSUPPORTED;
	}

	// compute_Flux_2nd
// will also need to be in the `has_complex_J` switch statement.
	if (flux_i_b->compute_Flux_2nd == NULL)
		flux_i->compute_Flux_2nd = NULL;
	else
		EXIT_UNSUPPORTED;
}
