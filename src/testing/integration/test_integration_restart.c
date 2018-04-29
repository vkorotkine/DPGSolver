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
#include "test_integration_convergence_support.h"

#include "adaptation.h"
#include "restart_writers.h"
#include "simulation.h"
#include "solution.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Output a restart file for the initial solution for the finest mesh level and polynomial order.
static void output_restart_finest
	(const char*const ctrl_name ///< The name of the control file.
	);

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

/* Would be better to output something significantly more accurate and run the convergence order tests normally.
 * Perhaps two input files here:
 * 1 - Fine mesh for restart writing. -> output_restart_finest
 * 2 - Relatively coarse sequence for convergence testing. -> run_convergence_order_study.
 */
	output_restart_finest(argv[1]);
	run_convergence_order_study(argc,argv,CONV_STUDY_RESTART);
EXIT_UNSUPPORTED; // Why is the convergence only optimal in entropy? THINK.

	PetscFinalize();
	OUTPUT_SUCCESS;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void output_restart_finest (const char*const ctrl_name)
{
	struct Integration_Test_Info* int_test_info = constructor_Integration_Test_Info(ctrl_name); // destructed

	const int p  = int_test_info->p_ref[0],
	          ml = int_test_info->ml[0],
	          p_prev  = p-1,
	          ml_prev = ml-1;

	const int adapt_type = int_test_info->adapt_type;
	const char*const ctrl_name_curr = set_file_name_curr(adapt_type,p,ml,false,ctrl_name);

	struct Simulation* sim = NULL;
	const char type_rc = 'r';
	structor_simulation(&sim,'c',adapt_type,p,ml,p_prev,ml_prev,ctrl_name_curr,type_rc,true); // destructed

	const struct Test_Case*const test_case = (struct Test_Case*) sim->test_case_rc->tc;
	if (!test_case->has_analytical)
		EXIT_ERROR("This test requires the use of a test case which has an analytical solution.");
	if (using_restart() && outputting_restart())
		EXIT_ERROR("%s %s","Setting the restarted solution from another restarted solution\n",
		                   "instead of the analytical solution.");

	adapt_to_maximum_refinement(sim,int_test_info);
	set_initial_solution(sim);
	if (outputting_restart())
		output_restart(sim);

	structor_simulation(&sim,'d',ADAPT_0,p,ml,p_prev,ml_prev,NULL,type_rc,false);

	destructor_Integration_Test_Info(int_test_info);
}
