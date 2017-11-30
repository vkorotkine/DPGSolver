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
 *  \brief Provides the templated linear advection flux functions.
 *  \todo Clean-up.
 */

#include <stddef.h>

#include "multiarray.h"

#include "flux.h"
#include "solution_advection.h"

// Static function declarations ************************************************************************************* //

#define NEQ 1 ///< Number of equations.

// Interface functions ********************************************************************************************** //

void compute_Flux_advection_T (const struct Flux_Input_T* flux_i, struct mutable_Flux_T* flux)
{
	struct Flux_Input_R* flux_i_r = (struct Flux_Input_R*) flux_i;

	static bool need_input = true;
	static struct Sol_Data__Advection sol_data;
	if (need_input) {
		need_input = false;
		read_data_advection(flux_i_r->input_path,&sol_data);
	}

	const ptrdiff_t NnTotal = flux_i->s->extents[0];

	Type const *const W = flux_i->s->data;
	Type       *const F = flux->f->data;

	Type *F_ptr[DIM*NEQ];
	for (int eq = 0; eq < NEQ; eq++)  {
	for (int dim = 0; dim < DIM; dim++) {
		F_ptr[eq*DIM+dim] = &F[(eq*DIM+dim)*NnTotal];
	}}

	const double* b_adv = sol_data.b_adv;
	for (int n = 0; n < NnTotal; n++) {
		for (int dim = 0; dim < DIM; dim++) {
			*F_ptr[dim] = b_adv[dim]*W[n];
			F_ptr[dim]++;
		}
	}
}

void compute_Flux_advection_jacobian_T (const struct Flux_Input_T* flux_i, struct mutable_Flux_T* flux)
{
	struct Flux_Input_R* flux_i_r = (struct Flux_Input_R*) flux_i;

	static bool need_input = true;
	static struct Sol_Data__Advection sol_data;
	if (need_input) {
		need_input = false;
		read_data_advection(flux_i_r->input_path,&sol_data);
	}

	const ptrdiff_t NnTotal = flux_i->s->extents[0];

	Type const *const W    = flux_i->s->data;
	Type       *const F    = flux->f->data;
	Type       *const dFdW = flux->df_ds->data;

	// Store pointers to the arrays that the data will be written into. Note: using NEQ == NVAR == 1.
	Type *F_ptr[DIM];
	if (F != NULL) {
		for (int dim = 0; dim < DIM; dim++)
			F_ptr[dim] = &F[dim*NnTotal];
	}

	Type *dFdW_ptr[DIM];
	for (int dim = 0; dim < DIM; dim++)
		dFdW_ptr[dim] = &dFdW[dim*NnTotal];

	const double* b_adv = sol_data.b_adv;
	for (int n = 0; n < NnTotal; n++) {
		if (F != NULL) {
			for (int dim = 0; dim < DIM; dim++) {
				*F_ptr[dim] = b_adv[dim]*W[n];
				F_ptr[dim]++;
			}
		}

		for (int dim = 0; dim < DIM; dim++) {
			*dFdW_ptr[dim] = b_adv[dim];
			dFdW_ptr[dim]++;
		}
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
