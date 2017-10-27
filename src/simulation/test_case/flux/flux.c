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

#include "flux.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"

#include "multiarray.h"

#include "const_cast.h"
#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

struct Flux_Input* constructor_Flux_Input (const struct Simulation* sim)
{
	struct Flux_Input* flux_i = calloc(1,sizeof *flux_i); // returned

	const_cast_c1(&flux_i->input_path,sim->input_path);

	struct Test_Case* test_case = sim->test_case;
	const_cast_i(&flux_i->d,sim->d);
	const_cast_i(&flux_i->n_eq,test_case->n_eq);
	const_cast_i(&flux_i->n_var,test_case->n_var);

	const_cast_b(&flux_i->has_1st_order,test_case->has_1st_order);
	const_cast_b(&flux_i->has_2nd_order,test_case->has_2nd_order);

	flux_i->compute_Flux = test_case->compute_Flux;
	switch (test_case->solver_method_curr) {
	case 'e':
		flux_i->compute_member   = test_case->flux_comp_mem_e;
		flux_i->compute_Flux_1st = test_case->compute_Flux_e[0];
		flux_i->compute_Flux_2nd = test_case->compute_Flux_e[1];
		break;
	case 'i':
		flux_i->compute_member   = test_case->flux_comp_mem_i;
		flux_i->compute_Flux_1st = test_case->compute_Flux_i[0];
		flux_i->compute_Flux_2nd = test_case->compute_Flux_i[1];
		break;
	default:
		EXIT_ERROR("Unsupported: %c.\n",test_case->solver_method_curr);
		break;
	}

	return flux_i;
}

void destructor_Flux_Input (struct Flux_Input* flux_i)
{
	free(flux_i);
}

struct Flux* constructor_Flux (const struct Flux_Input* flux_i)
{
	assert((flux_i->s != NULL && flux_i->s->layout == 'C') ||
	       (flux_i->g != NULL && flux_i->g->layout == 'C'));

	const bool* compute_member = flux_i->compute_member;
	assert(compute_member[0] || compute_member[1] || compute_member[2]);

	const int d     = flux_i->d,
	          n_eq  = flux_i->n_eq,
		    n_var = flux_i->n_var;
	const ptrdiff_t n_n = ( flux_i->s != NULL ? flux_i->s->extents[0] : flux_i->g->extents[0] );

	struct mutable_Flux* flux = calloc(1,sizeof *flux); // returned

	flux->f     = (compute_member[0] ? constructor_zero_Multiarray_d('C',3,(ptrdiff_t[]){n_n,d,n_eq})         : NULL);
	flux->df_ds = (compute_member[1] ? constructor_zero_Multiarray_d('C',4,(ptrdiff_t[]){n_n,d,n_eq,n_var})   : NULL);
	flux->df_dg = (compute_member[2] ? constructor_zero_Multiarray_d('C',5,(ptrdiff_t[]){n_n,d,n_eq,n_var,d}) : NULL);

	flux_i->compute_Flux(flux_i,flux);

	return (struct Flux*) flux;
}

void destructor_Flux (struct Flux* flux)
{
	if (flux->f != NULL)
		destructor_const_Multiarray_d(flux->f);
	if (flux->df_ds != NULL)
		destructor_const_Multiarray_d(flux->df_ds);
	if (flux->df_dg != NULL)
		destructor_const_Multiarray_d(flux->df_dg);

	free(flux);
}

void compute_Flux_1 (const struct Flux_Input* flux_i, struct mutable_Flux* flux)
{
	flux_i->compute_Flux_1st(flux_i,flux);
}

void compute_Flux_12 (const struct Flux_Input* flux_i, struct mutable_Flux* flux)
{
	flux_i->compute_Flux_1st(flux_i,flux);
	flux_i->compute_Flux_2nd(flux_i,flux);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
