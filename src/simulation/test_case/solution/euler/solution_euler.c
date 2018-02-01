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

#include "solution_euler.h"

#include "multiarray.h"

#include "boundary.h"
#include "compute_error_euler.h"
#include "const_cast.h"
#include "flux_euler.h"
#include "geometry.h"
#include "geometry_parametric.h"
#include "numerical_flux_euler.h"
#include "simulation.h"
#include "solution.h"
#include "test_case.h"

#include "periodic_vortex/solution_periodic_vortex.h"
#include "supersonic_vortex/solution_supersonic_vortex.h"
#include "free_stream/solution_free_stream.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "solution_euler_T.c"

void compute_entropy (struct Multiarray_d* s, const struct const_Multiarray_d* vars, const char var_type)
{
	assert(s->extents[0] == vars->extents[0]);
	assert(s->extents[1] == 1);
	assert(vars->layout == 'C');

	if (var_type != 'p')
		convert_variables((struct Multiarray_d*)vars,var_type,'p');

	const double* rho = get_col_const_Multiarray_d(0,vars),
	            * p   = get_col_const_Multiarray_d(vars->extents[1]-1,vars);

	const ptrdiff_t ext_0 = s->extents[0];
	for (ptrdiff_t i = 0; i < ext_0; ++i)
//		s->data[i] = p[i]*pow(rho[i],-GAMMA);
		s->data[i] = log(p[i]*pow(rho[i],-GAMMA));

	if (var_type != 'p')
		convert_variables((struct Multiarray_d*)vars,'p',var_type);
}

void compute_mach (struct Multiarray_d* mach, const struct const_Multiarray_d* vars, const char var_type)
{
	assert(mach->extents[0] == vars->extents[0]);
	assert(mach->extents[1] == 1);
	assert(vars->layout == 'C');

	if (var_type != 'p')
		convert_variables((struct Multiarray_d*)vars,var_type,'p');

	const ptrdiff_t d = vars->extents[1]-2;

	const double* rho = get_col_const_Multiarray_d(0,vars),
		      * u   = get_col_const_Multiarray_d(1,vars),
		      * v   = (d > 1 ? get_col_const_Multiarray_d(2,vars) : NULL),
		      * w   = (d > 2 ? get_col_const_Multiarray_d(3,vars) : NULL),
	            * p   = get_col_const_Multiarray_d(vars->extents[1]-1,vars);

	const ptrdiff_t ext_0 = mach->extents[0];
	for (ptrdiff_t i = 0; i < ext_0; ++i) {
		double V2 = 0.0;
		switch (d) {
			case 3: V2 += w[i]*w[i]; // fallthrough
			case 2: V2 += v[i]*v[i]; // fallthrough
			case 1: V2 += u[i]*u[i]; break;
			default: EXIT_ERROR("Unsupported: %td\n",d); break;
		}
		const double c2 = GAMMA*p[i]/rho[i];
		mach->data[i] = sqrt(V2/c2);
	}

	if (var_type != 'p')
		convert_variables((struct Multiarray_d*)vars,'p',var_type);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
