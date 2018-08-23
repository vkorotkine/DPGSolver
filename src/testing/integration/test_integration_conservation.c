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
#include "petscsys.h"

#include "macros.h"
#include "definitions_test_case.h"
#include "definitions_tol.h"

#include "test_base.h"
#include "test_integration.h"

#include "vector.h"

#include "volume_solver.h"

#include "const_cast.h"
#include "core.h"
#include "intrusive.h"
#include "math_functions.h"
#include "simulation.h"
#include "solve.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// Flag for whether the conservation of the partially converged solution should be checked.
#define CHECK_PARTIALLY_CONVERGED false

/// \brief Check the flux imbalance magnitudes in each element and over the entire domain.
static void check_flux_imbalance
	(struct Test_Info*const test_info, ///< \ref Test_Info.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

/** \test Performs integration testing for conservation of the various solver methods (\ref
 *        test_integration_conservation.c).
 *  \return 0 on success (when the element-local and global conservation is achieved (negligible flux imbalance)).
 */
int main
	(int argc,   ///< Standard.
	 char** argv ///< Standard.
	)
{
	assert_condition_message(argc == 3,"Invalid number of input arguments");

	const char* petsc_options_name = set_petsc_options_name(argv[2]);
	PetscInitialize(&argc,&argv,petsc_options_name,PETSC_NULL);

	const char* ctrl_name = argv[1];

	struct Test_Info test_info = { .n_warn = 0, };

	struct Integration_Test_Info* int_test_info = constructor_Integration_Test_Info(ctrl_name);

	const int p  = int_test_info->p_ref[0],
	          ml = int_test_info->ml[0],
	          p_prev  = p-1,
	          ml_prev = ml-1;

	const int adapt_type = int_test_info->adapt_type;
	const char*const ctrl_name_curr = set_file_name_curr(adapt_type,p,ml,false,ctrl_name);

	struct Simulation* sim = NULL;
	structor_simulation(&sim,'c',adapt_type,p,ml,p_prev,ml_prev,ctrl_name_curr,'r',false); // destructed

	solve_for_solution(sim);
	compute_flux_imbalances(sim);

	check_flux_imbalance(&test_info,sim);

	structor_simulation(&sim,'d',adapt_type,p,ml,p_prev,ml_prev,NULL,'r',false);

	destructor_Integration_Test_Info(int_test_info);

	PetscFinalize();
	OUTPUT_SUCCESS;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void check_flux_imbalance (struct Test_Info*const test_info, const struct Simulation*const sim)
{
	UNUSED(test_info);
	bool pass = true;
	const double tol = 5e1*EPS;

	struct Test_Case* test_case = (struct Test_Case*)sim->test_case_rc->tc;

	struct Vector_d*const flux_imbalance_global = constructor_zero_Vector_d(test_case->n_var); // destructed
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		const struct Solver_Volume*const s_vol = (struct Solver_Volume*) curr;

		const struct const_Vector_d*const flux_imbalance = (struct const_Vector_d*) s_vol->flux_imbalance;

		add_to_Vector_d_d(flux_imbalance_global,flux_imbalance->data);
		if (maximum_abs_d(flux_imbalance->data,flux_imbalance->ext_0) > tol) {
			pass = false;

			const struct Volume*const vol = (struct Volume*) curr;
			printf("Local imbalance (vol: %d): \n",vol->index);
			print_const_Vector_d(flux_imbalance);
		}
	}
	if (maximum_abs_d(flux_imbalance_global->data,flux_imbalance_global->ext_0) > tol) {
		pass = false;
		printf("Global imbalance:\n");
		print_Vector_d(flux_imbalance_global);
	}
	destructor_Vector_d(flux_imbalance_global);

	assert_condition(pass);
}
