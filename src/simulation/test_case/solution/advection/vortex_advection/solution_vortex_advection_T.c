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
#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Return a \ref Multiarray_T\* container holding the solution values at the input coordinates.
 *  \return See brief. */
static struct Multiarray_T* constructor_sol_vortex_advection
	(const struct Simulation* sim,        ///< Defined for \ref set_sol_vortex_advection_T.
	 const struct const_Multiarray_R* xyz ///< xyz coordinates at which to evaluate the solution.
	);

// Interface functions ********************************************************************************************** //

void set_sol_vortex_advection_T (const struct Simulation* sim, struct Solution_Container_T sol_cont)
{
	const struct const_Multiarray_R* xyz = constructor_xyz_sol_T(sim,&sol_cont); // destructed
	struct Multiarray_T* sol = constructor_sol_vortex_advection(sim,xyz); // destructed
	destructor_const_Multiarray_R(xyz);

	update_Solution_Container_sol_T(&sol_cont,sol,sim);
	destructor_Multiarray_T(sol);
}

const struct const_Multiarray_T* constructor_const_sol_vortex_advection_T
	(const struct const_Multiarray_R* xyz, const struct Simulation* sim)
{
	struct Multiarray_T* sol = constructor_sol_vortex_advection(sim,xyz); // returned
	return (const struct const_Multiarray_T*) sol;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct Multiarray_T* constructor_sol_vortex_advection
	(const struct Simulation* sim, const struct const_Multiarray_R* xyz)
{
	assert(DIM == 2);

	static int adv_type = -1;
	static bool use_constant_solution = false;

	static bool requires_input = true;
	if (requires_input) {
		const struct Sol_Data__Advection sol_data = get_sol_data_advection();
		if (sol_data.compute_b_adv == compute_b_adv_vortex)
			adv_type = ADVECTION_TYPE_VORTEX;
		else if (sol_data.compute_b_adv == compute_b_adv_constant)
			adv_type = ADVECTION_TYPE_CONST;
		else
			EXIT_UNSUPPORTED;
		use_constant_solution = sol_data.use_constant_solution;
	}

	// Compute the solution
	const ptrdiff_t n_n = xyz->extents[0];
	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const int n_var = test_case->n_var;

	struct Multiarray_T* sol = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){n_n,n_var}); // returned

	const Real* x = get_col_const_Multiarray_R(0,xyz),
	          * y = get_col_const_Multiarray_R(1,xyz);

	Type* u = get_col_Multiarray_T(0,sol);

	const struct Sol_Data__Advection sol_data = get_sol_data_advection();
	const double scale = sol_data.u_scale;
	assert(scale != 0.0);
	switch (adv_type) {
	case ADVECTION_TYPE_VORTEX:
		if (!use_constant_solution) {
			for (int i = 0; i < n_n; ++i) {
				const Real r  = sqrt(x[i]*x[i]+y[i]*y[i]);
				u[i] = scale*sin(1.1e1*r)*cos(1.3e1*r);
			}
		} else {
			for (int i = 0; i < n_n; ++i)
				u[i] = scale;
		}
		break;
	case ADVECTION_TYPE_CONST: {
		for (int i = 0; i < n_n; ++i)
			u[i] = scale*sin(0.1*y[i])*cos(0.3*y[i]);
		break;
	} default:
		EXIT_ERROR("Unsupported: %d\n",adv_type);
		break;
	}

	return sol;
}

#include "undef_templates_solution.h"
#include "undef_templates_solution_advection.h"

#include "undef_templates_multiarray.h"
#include "undef_templates_test_case.h"
