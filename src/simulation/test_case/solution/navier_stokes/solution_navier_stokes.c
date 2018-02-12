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

#include "solution_navier_stokes.h"

#include "multiarray.h"

#include "boundary.h"
#include "compute_error_navier_stokes.h"
#include "file_processing.h"
#include "flux_euler.h"
#include "flux_navier_stokes.h"
#include "geometry.h"
#include "geometry_parametric.h"
#include "numerical_flux_euler.h"
#include "numerical_flux_navier_stokes.h"
#include "simulation.h"
#include "solution.h"
#include "solution_euler.h"
#include "test_case.h"

#include "taylor_couette/solution_taylor_couette.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "solution_navier_stokes_T.c"

double compute_kappa_const_cp (const double mu, const double Cp, const double Pr)
{
	return mu*Cp/Pr;
}

double compute_cp_ideal_gas (const double r_s)
{
	return GAMMA/GM1*r_s;
}

void compute_viscosity
	(struct Multiarray_d* mu, const struct const_Multiarray_d* vars, const char var_type, const char*const input_path)
{
	assert(mu->extents[0] == vars->extents[0]);
	assert(mu->extents[1] == 1);
	assert(vars->extents[1]-2 == DIM);
	assert(vars->layout == 'C');

	compute_mu_fptr_d compute_mu = get_compute_mu_fptr_d(input_path);

	if (var_type != 'c')
		convert_variables((struct Multiarray_d*)vars,var_type,'c');

	const double*const rho         = get_col_const_Multiarray_d(0,vars),
	            *const rhouvw[DIM] = ARRAY_DIM( get_col_const_Multiarray_d(1,vars),
		                                      get_col_const_Multiarray_d(2,vars),
		                                      get_col_const_Multiarray_d(3,vars) ),
	            *const E           = get_col_const_Multiarray_d(vars->extents[1]-1,vars);

	const ptrdiff_t ext_0 = mu->extents[0];
	for (ptrdiff_t i = 0; i < ext_0; ++i) {
		const double rhouvw_i[DIM] = ARRAY_DIM( rhouvw[0][i], rhouvw[1][i], rhouvw[2][i] );
		mu->data[i] = compute_mu(input_path,rho[i],rhouvw_i,E[i]);
	}

	if (var_type != 'c')
		convert_variables((struct Multiarray_d*)vars,'c',var_type);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
