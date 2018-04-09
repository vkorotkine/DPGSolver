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
 *  \brief Provides the templated linear advection numerical flux functions.
 */

#include <assert.h>
#include <stddef.h>
#include <math.h>

#include "macros.h"
#include "definitions_bc.h"


#include "def_templates_multiarray.h"

#include "def_templates_numerical_flux.h"

#include "def_templates_boundary.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void compute_Numerical_Flux_T_advection_upwind
	(const struct Numerical_Flux_Input_T* num_flux_i, struct mutable_Numerical_Flux_T* num_flux)
{
	static bool need_input = true;
	static struct Sol_Data__Advection sol_data;
	if (need_input) {
		need_input = false;
		read_data_advection(&sol_data);
	}

	const ptrdiff_t NnTotal = num_flux_i->bv_l.s->extents[0];

	double const *const nL = num_flux_i->bv_l.normals->data;

	Type const *const WL = num_flux_i->bv_l.s->data,
	           *const WR = num_flux_i->bv_r.s->data;

	Type       *const nFluxNum = num_flux->nnf->data;

	const struct const_Multiarray_d*const xyz_Ma = num_flux_i->bv_l.xyz;
	const Real*const xyz[DIM] = ARRAY_DIM( get_col_const_Multiarray_R(0,xyz_Ma),
	                                       get_col_const_Multiarray_R(1,xyz_Ma),
	                                       get_col_const_Multiarray_R(2,xyz_Ma) );

	for (int n = 0; n < NnTotal; n++) {
		const Real xyz_n[DIM] = ARRAY_DIM(xyz[0][n],xyz[1][n],xyz[2][n]);
		const double*const b_adv = sol_data.compute_b_adv(xyz_n);

		double b_dot_n = 0.0;
		for (int dim = 0; dim < DIM; dim++)
			b_dot_n += b_adv[dim]*nL[n*DIM+dim];

		if (b_dot_n >= 0.0)
			nFluxNum[n] = b_dot_n*WL[n];
		else
			nFluxNum[n] = b_dot_n*WR[n];
	}
}

void compute_Numerical_Flux_T_advection_upwind_jacobian
	(const struct Numerical_Flux_Input_T* num_flux_i, struct mutable_Numerical_Flux_T* num_flux
	)
{
	static bool need_input = true;
	static struct Sol_Data__Advection sol_data;
	if (need_input) {
		need_input = false;
		read_data_advection(&sol_data);
	}

	const ptrdiff_t NnTotal = num_flux_i->bv_l.s->extents[0];

	double const *const nL = num_flux_i->bv_l.normals->data;

	Type const *const WL = num_flux_i->bv_l.s->data,
	           *const WR = num_flux_i->bv_r.s->data;

	Type       *const nFluxNum     = num_flux->nnf->data,
	           *const dnFluxNumdWL = num_flux->neigh_info[0].dnnf_ds->data,
	           *const dnFluxNumdWR = num_flux->neigh_info[1].dnnf_ds->data;
	assert(nFluxNum     != NULL);
	assert(dnFluxNumdWL != NULL);
	assert(dnFluxNumdWR != NULL);

	const struct const_Multiarray_d*const xyz_Ma = num_flux_i->bv_l.xyz;
	const Real*const xyz[DIM] = ARRAY_DIM( get_col_const_Multiarray_R(0,xyz_Ma),
	                                       get_col_const_Multiarray_R(1,xyz_Ma),
	                                       get_col_const_Multiarray_R(2,xyz_Ma) );

const int bc = num_flux_i->bv_l.bc % BC_STEP_SC;
	for (int n = 0; n < NnTotal; n++) {
		const Real xyz_n[DIM] = ARRAY_DIM(xyz[0][n],xyz[1][n],xyz[2][n]);
		const double*const b_adv = sol_data.compute_b_adv(xyz_n);

		double b_dot_n = 0.0;
		for (int dim = 0; dim < DIM; dim++)
			b_dot_n += b_adv[dim]*nL[n*DIM+dim];

		if (b_dot_n >= 0.0) {
			nFluxNum[n]     = b_dot_n*WL[n];
			dnFluxNumdWL[n] = b_dot_n;
			dnFluxNumdWR[n] = 0.0;
		} else {
switch (bc) {
case BC_SLIPWALL: {
#if TYPE_RC == TYPE_REAL
//printf("nf: % .3e % .3e % .3e % .3e % .3e\n",b_dot_n,WL[n]-WR[n],b_dot_n*(WL[n]-WR[n]),b_dot_n*WL[n],b_dot_n*WR[n]);
#endif
	const double h = num_flux_i->bv_l.h;
	const int p = 3;
	UNUSED(h);

	nFluxNum[n]     = b_dot_n*WL[n] + pow(h,p+1);
	dnFluxNumdWL[n] = b_dot_n;
	dnFluxNumdWR[n] = 0.0;
	break;
} default:
			nFluxNum[n]     = b_dot_n*WR[n];
			dnFluxNumdWL[n] = 0.0;
			dnFluxNumdWR[n] = b_dot_n;
	break;
}
		}
	}
//if (bc == BC_SLIPWALL)
//printf("\n");
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
