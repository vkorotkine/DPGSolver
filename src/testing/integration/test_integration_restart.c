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
#include "definitions_adaptation.h"

#include "test_base.h"
#include "test_integration.h"

#include "restart_writers.h"
#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

/** \test Performs integration testing for restarted solutions (\ref test_integration_restart.c).
 *  \return 0 on success.
 *
 *  Requiring a test case with a non-constant specified exact solution, this test:
 *  1. Uniformly refines (h/p) the input mesh several times, re-initializes the solution and outputs a restart file.
 *  2. Loops over the all degrees and all but the finest mesh level and computes the error between the restarted and
 *     exact solution.
 *
 *  The test passes if the convergence is optimal (O(h^{p+1}).
 */
int main
	(int argc,   ///< Standard.
	 char** argv ///< Standard.
	)
{
	PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

	assert_condition_message(argc == 2,"Invalid number of input arguments");

	const char* ctrl_name = argv[1];

	struct Integration_Test_Info* int_test_info = constructor_Integration_Test_Info(ctrl_name);

	const int p  = int_test_info->p_ref[0],
	          ml = int_test_info->ml[0],
	          p_prev  = p-1,
	          ml_prev = ml-1;

	const int adapt_type = int_test_info->adapt_type;
	const char*const ctrl_name_curr = set_file_name_curr(adapt_type,p,ml,false,ctrl_name);

	struct Simulation* sim = NULL;
	structor_simulation(&sim,'c',adapt_type,p,ml,p_prev,ml_prev,ctrl_name_curr,'r'); // destructed
	adapt_initial_mesh_if_required(sim);

/// \todo Update the test when initial implementation is working.
struct Test_Info test_info = { .n_warn = 0, };
test_print_warning(&test_info,"Not performing full uniform hp refinement yet.");

	const struct Test_Case*const test_case = (struct Test_Case*) sim->test_case_rc->tc;
	if (!test_case->has_analytical)
		EXIT_ERROR("This test requires the use of a test case which has an analytical solution.");

	const struct Restart_Info restart_info = { .ml = ml, .p = p, };
	output_restart(sim,&restart_info);

EXIT_UNSUPPORTED;

	structor_simulation(&sim,'d',ADAPT_0,p,ml,p_prev,ml_prev,NULL,'r');

	destructor_Integration_Test_Info(int_test_info);

	PetscFinalize();
	OUTPUT_SUCCESS;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
