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
 *  \todo Clean-up (potentially template functions).
 */

#include "test_complex_flux_advection.h"

#include <assert.h>

#include "complex_multiarray_minimal.h"

#include "test_complex_flux.h"
#include "solution_advection.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void compute_Flux_c_advection (const struct Flux_Input_c* flux_i, struct mutable_Flux_c* flux)
{
	struct Flux_Input* flux_i_b = (struct Flux_Input*) flux_i;

	static bool need_input = true;
	static struct Sol_Data__Advection sol_data;
	if (need_input) {
		need_input = false;
		read_data_advection(flux_i_b->input_path,&sol_data);
	}

	int const d   = flux_i_b->d,
	          Neq = 1;
	const ptrdiff_t NnTotal = flux_i->s->extents[0];

	double complex const *const W = flux_i->s->data;
	double complex       *const F = flux->f->data;

	double complex *F_ptr[d*Neq];
	for (int eq = 0; eq < Neq; eq++)  {
	for (int dim = 0; dim < d; dim++) {
		F_ptr[eq*d+dim] = &F[(eq*d+dim)*NnTotal];
	}}

	const double* b_adv = sol_data.b_adv;
	for (int n = 0; n < NnTotal; n++) {
		for (int dim = 0; dim < d; dim++) {
			*F_ptr[dim] = b_adv[dim]*W[n];
			F_ptr[dim]++;
		}
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
