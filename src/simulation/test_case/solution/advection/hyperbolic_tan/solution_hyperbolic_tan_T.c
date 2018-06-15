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
#include "definitions_core.h"
#include "definitions_solution.h"
#include "definitions_tol.h"


#include "def_templates_solution.h"
#include "def_templates_solution_advection.h"

#include "def_templates_multiarray.h"
#include "def_templates_math_functions.h"
#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Return a \ref Multiarray_T\* container holding the solution values at the input coordinates.
 *  \return See brief. */
static struct Multiarray_T* constructor_sol_hyperbolic_tan
	(const struct Simulation* sim,        ///< Defined for \ref set_sol_hyperbolic_tan_T.
	 const struct const_Multiarray_T* xyz ///< xyz coordinates at which to evaluate the solution.
	);

// Interface functions ********************************************************************************************** //

void set_sol_hyperbolic_tan_T (const struct Simulation* sim, struct Solution_Container_T sol_cont)
{
	const struct const_Multiarray_T* xyz = constructor_xyz_sol_T(sim,&sol_cont); // destructed
	struct Multiarray_T* sol = constructor_sol_hyperbolic_tan(sim,xyz); // destructed
	destructor_const_Multiarray_T(xyz);

	update_Solution_Container_sol_T(&sol_cont,sol,sim);
	destructor_Multiarray_T(sol);
}

const struct const_Multiarray_T* constructor_const_sol_hyperbolic_tan_T
	(const struct const_Multiarray_T* xyz, const struct Simulation* sim)
{
	struct Multiarray_T* sol = constructor_sol_hyperbolic_tan(sim,xyz); // returned
	return (const struct const_Multiarray_T*) sol;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct Multiarray_T* constructor_sol_hyperbolic_tan
	(const struct Simulation* sim, const struct const_Multiarray_T* xyz)
{
	assert(DIM == 1);

	const ptrdiff_t n_n = xyz->extents[0];
	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const int n_var   = test_case->n_var;
	const double time = test_case->time;

	struct Multiarray_T* sol = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){n_n,n_var}); // returned

	const Type* x = get_col_const_Multiarray_T(0,xyz);

	Type* u = get_col_Multiarray_T(0,sol);

	const struct Sol_Data__Advection_T sol_data = get_sol_data_advection_T();
	assert(sol_data.compute_b_adv == compute_b_adv_constant_T);

	const Real scale = 250.0*1e-3;
	for (int i = 0; i < n_n; ++i) {
		const double*const b_adv = sol_data.compute_b_adv(&x[i]);
		const double x_t = real_T(x[i]) - b_adv[0]*time;
		u[i] = 0.5*(1.0+tanh(scale*(x_t-20.0)));
	}

	return sol;
}

#include "undef_templates_solution.h"
#include "undef_templates_solution_advection.h"

#include "undef_templates_multiarray.h"
#include "undef_templates_math_functions.h"
#include "undef_templates_test_case.h"
