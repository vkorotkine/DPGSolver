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

#include "test_complex_numerical_flux_advection.h"

#include <assert.h>
#include <stdlib.h>

#include "macros.h"
#include "definitions_test_case.h"

#include "test_complex_numerical_flux.h"

#include "complex_multiarray.h"
#include "multiarray.h"

#include "const_cast.h"
#include "flux.h"
#include "numerical_flux.h"
#include "solution_advection.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void compute_Numerical_Flux_c_advection_upwind
	(const struct Numerical_Flux_Input_c* num_flux_i, struct mutable_Numerical_Flux_c* num_flux)
{
	struct Numerical_Flux_Input* num_flux_i_b = (struct Numerical_Flux_Input*) num_flux_i;
	struct Boundary_Value_Input* bv_l_b       = (struct Boundary_Value_Input*) &num_flux_i->bv_l;

	static bool need_input = true;
	static struct Sol_Data__Advection sol_data;
	if (need_input) {
		need_input = false;
		read_data_advection(num_flux_i_b->bv_l.input_path,&sol_data);
	}

	int const d   = num_flux_i_b->bv_l.d;
	const ptrdiff_t NnTotal = num_flux_i->bv_l.s->extents[0];

	double const *const nL = bv_l_b->normals->data;

	double complex const *const WL = num_flux_i->bv_l.s->data,
	                     *const WR = num_flux_i->bv_r.s->data;

	double complex       *const nFluxNum = num_flux->nnf->data;

	const double* b_adv = sol_data.b_adv;
	for (int n = 0; n < NnTotal; n++) {
		double b_dot_n = 0.0;
		for (int dim = 0; dim < d; dim++)
			b_dot_n += b_adv[dim]*nL[n*d+dim];

		if (b_dot_n >= 0.0)
			nFluxNum[n] = b_dot_n*WL[n];
		else
			nFluxNum[n] = b_dot_n*WR[n];
	}
}

void compute_Numerical_Flux_c_advection_upwind_jacobian
	(const struct Numerical_Flux_Input_c* num_flux_i, struct mutable_Numerical_Flux_c* num_flux)
{
	struct Numerical_Flux_Input* num_flux_i_b = (struct Numerical_Flux_Input*) num_flux_i;
	struct Boundary_Value_Input* bv_l_b       = (struct Boundary_Value_Input*) &num_flux_i->bv_l;

	static bool need_input = true;
	static struct Sol_Data__Advection sol_data;
	if (need_input) {
		need_input = false;
		read_data_advection(num_flux_i_b->bv_l.input_path,&sol_data);
	}

	int const d   = num_flux_i_b->bv_l.d;
	const ptrdiff_t NnTotal = num_flux_i->bv_l.s->extents[0];

	double const *const nL = bv_l_b->normals->data;

	double complex const *const WL = num_flux_i->bv_l.s->data,
	                     *const WR = num_flux_i->bv_r.s->data;

	double complex      *const nFluxNum     = num_flux->nnf->data,
	                    *const dnFluxNumdWL = num_flux->neigh_info[0].dnnf_ds->data,
	                    *const dnFluxNumdWR = num_flux->neigh_info[1].dnnf_ds->data;
	assert(nFluxNum     != NULL);
	assert(dnFluxNumdWL != NULL);
	assert(dnFluxNumdWR != NULL);

	const double* b_adv = sol_data.b_adv;
	for (int n = 0; n < NnTotal; n++) {
		double b_dot_n = 0.0;
		for (int dim = 0; dim < d; dim++)
			b_dot_n += b_adv[dim]*nL[n*d+dim];

		if (b_dot_n >= 0.0) {
			nFluxNum[n]     = b_dot_n*WL[n];
			dnFluxNumdWL[n] = b_dot_n;
			dnFluxNumdWR[n] = 0.0;
		} else {
			nFluxNum[n]     = b_dot_n*WR[n];
			dnFluxNumdWL[n] = 0.0;
			dnFluxNumdWR[n] = b_dot_n;
		}
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
