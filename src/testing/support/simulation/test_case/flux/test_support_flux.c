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

#include "test_support_flux.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "macros.h"

#include "test_support_flux_advection.h"

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

	struct Flux_Input_c* flux_i = calloc(1,sizeof *flux_i); // destructed.

	memcpy(flux_i,flux_i_b,sizeof(struct Flux_Input)); // shallow copy of the base.
	destructor_Flux_Input(flux_i_b);

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

	const int d    = flux_i_b->d,
	          n_eq = flux_i_b->n_eq;
	const ptrdiff_t n_n = ( flux_i->s != NULL ? flux_i->s->extents[0] : flux_i->g->extents[0] );

	struct mutable_Flux_c* flux = calloc(1,sizeof *flux); // destructed

	flux->f = constructor_zero_Multiarray_c('C',3,(ptrdiff_t[]){n_n,d,n_eq});

	flux_i->compute_Flux(flux_i,flux);

	return (struct Flux_c*) flux;
}

void destructor_Flux_c (struct Flux_c* flux)
{
	destructor_const_Multiarray_c(flux->f);
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
	if (flux_i_b->compute_Flux_1st == compute_Flux_advection_jacobian)
		flux_i->compute_Flux_1st = compute_Flux_c_advection;
	else if (flux_i_b->compute_Flux_1st == compute_Flux_euler_jacobian)
		EXIT_ADD_SUPPORT;
	else if (flux_i_b->compute_Flux_1st == NULL)
		flux_i->compute_Flux_1st = NULL;
	else
		EXIT_UNSUPPORTED;

	// compute_Flux_2nd
	if (flux_i_b->compute_Flux_2nd == NULL)
		flux_i->compute_Flux_2nd = NULL;
	else
		EXIT_UNSUPPORTED;
}
