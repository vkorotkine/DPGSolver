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

#include <assert.h>
#include <math.h>
#include <string.h>

#include "macros.h"
#include "definitions_intrusive.h"
#include "definitions_math.h"
#include "definitions_tol.h"


#include "def_templates_solution.h"
#include "def_templates_solution_euler.h"

#include "def_templates_multiarray.h"
#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Return a \ref Multiarray_T\* container holding the solution values at the input coordinates.
 *  \return See brief. */
static struct Multiarray_T* constructor_sol_periodic_vortex_T
	(const struct const_Multiarray_R* xyz, ///< xyz coordinates at which to evaluate the solution.
	 const struct Simulation* sim          ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void set_sol_periodic_vortex (const struct Simulation* sim, struct Solution_Container_T sol_cont)
{
	const struct const_Multiarray_R* xyz = constructor_xyz_sol(sim,&sol_cont); // destructed
	struct Multiarray_T* sol = constructor_sol_periodic_vortex_T(xyz,sim); // destructed
	destructor_const_Multiarray_R(xyz);

	update_Solution_Container_sol_T(&sol_cont,sol);
	destructor_Multiarray_T(sol);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct Multiarray_T* constructor_sol_periodic_vortex_T
	(const struct const_Multiarray_R* xyz, const struct Simulation* sim)
{
	assert(DIM >= 2);
	const struct Sol_Data__pv sol_data = get_sol_data(sim);

	// Set the coordinates of the vortex centre depending on the time.
	Real x_c = 0.0,
	     y_c = 0.0;
	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	set_xy_c(&x_c,&y_c,&sol_data,test_case->time);

	// Compute the solution
	const ptrdiff_t n_n = xyz->extents[0];
	assert(DIM == xyz->extents[1]);

	const Real* x = get_col_const_Multiarray_R(0,xyz),
	          * y = get_col_const_Multiarray_R(1,xyz);

	const int n_var = test_case->n_var;

	struct Multiarray_T* sol = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){n_n,n_var}); // returned

	Type* rho = get_col_Multiarray_T(0,sol),
	    * u   = get_col_Multiarray_T(1,sol),
	    * v   = get_col_Multiarray_T(2,sol),
	    * p   = get_col_Multiarray_T(n_var-1,sol);
	for (int i = 0; i < n_n; ++i) {
		const Real rho_inf = sol_data.rho_inf,
		           u_inf   = sol_data.u_inf,
		           v_inf   = sol_data.v_inf,
		           p_inf   = sol_data.p_inf,
		           r_v     = sol_data.r_v,
		           con     = sol_data.con;

		const Real r2 = (pow(x[i]-x_c,2.0)+pow(y[i]-y_c,2.0))/(r_v*r_v);

		rho[i] = rho_inf;
		u[i]   = u_inf - con*(y[i]-y_c)/(r_v*r_v)*exp(-0.5*r2);
		v[i]   = v_inf + con*(x[i]-x_c)/(r_v*r_v)*exp(-0.5*r2);
		p[i]   = p_inf - rho_inf*(con*con)/(2*r_v*r_v)*exp(-r2);
	}

	if (DIM == 3) {
		Type* w = get_col_Multiarray_T(3,sol);
		for (int i = 0; i < n_n; ++i)
			w[i] = 0.0;
	}
	convert_variables_T(sol,'p','c');

	return sol;
}
