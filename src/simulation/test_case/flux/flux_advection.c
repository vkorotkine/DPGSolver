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
 *  \todo Clean-up.
 */

#include "flux_advection.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"
#include "definitions_core.h"
#include "definitions_test_case.h"

#include "multiarray.h"

#include "flux.h"
#include "solution_advection.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void compute_Flux_advection (const struct Flux_Input* flux_i, struct mutable_Flux* flux)
{
	static bool need_input = true;
	static struct Sol_Data__Advection sol_data;
	if (need_input) {
		need_input = false;
		read_data_advection(flux_i->input_path,&sol_data);
	}

	int const d   = flux_i->d,
	          Neq = 1;
	const ptrdiff_t NnTotal = flux_i->s->extents[0];

	double const *const W = flux_i->s->data;
	double       *const F = flux->f->data;

	double *F_ptr[d*Neq];
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

void compute_Flux_advection_jacobian (const struct Flux_Input* flux_i, struct mutable_Flux* flux)
{
	static bool need_input = true;
	static struct Sol_Data__Advection sol_data;
	if (need_input) {
		need_input = false;
		read_data_advection(flux_i->input_path,&sol_data);
	}

	int const d = flux_i->d;
	const ptrdiff_t NnTotal = flux_i->s->extents[0];

	double const *const W    = flux_i->s->data;
	double       *const F    = flux->f->data;
	double       *const dFdW = flux->df_ds->data;

	// Store pointers to the arrays that the data will be written into. Note: using Neq == Nvar == 1.
	double *F_ptr[d];
	if (F != NULL) {
		for (int dim = 0; dim < d; dim++)
			F_ptr[dim] = &F[dim*NnTotal];
	}

	double *dFdW_ptr[d];
	for (int dim = 0; dim < d; dim++)
		dFdW_ptr[dim] = &dFdW[dim*NnTotal];

	const double* b_adv = sol_data.b_adv;
	for (int n = 0; n < NnTotal; n++) {
		if (F != NULL) {
			for (int dim = 0; dim < d; dim++) {
				*F_ptr[dim] = b_adv[dim]*W[n];
				F_ptr[dim]++;
			}
		}

		for (int dim = 0; dim < d; dim++) {
			*dFdW_ptr[dim] = b_adv[dim];
			dFdW_ptr[dim]++;
		}
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
