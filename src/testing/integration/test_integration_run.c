
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

#include "petscsys.h"

#include "macros.h"

#include "test_base.h"
#include "test_integration.h"
#include "test_integration_convergence_support.h"

#include "solve.h"
#include "definitions_visualization.h"
#include "visualization.h"
#include "compute_error.h"

#include "core.h"

// Static function declarations ************************************************************************************* //

static void run_case (int argc, char** argv);

// Interface functions ********************************************************************************************** //

/** \test Performs integration testing for the solution convergence (\ref test_integration_run.c).
 *  \return 0 on success (when the solution converges).
 */
int main
	(int argc,   ///< Standard.
	 char** argv ///< Standard.
	)
{
	// Ensure that this is a test case with only one command line argument
	assert_condition_message((argc == 3) || (argc == 4),"Invalid number of input arguments");

	// Initialize PETSc
	const char* petsc_options_name = set_petsc_options_name(argv[2]);
	PetscInitialize(&argc,&argv,petsc_options_name,PETSC_NULL);

	run_case(argc,argv);

	PetscFinalize();
	OUTPUT_SUCCESS;
}

// Static functions ************************************************************************************************* //

static void run_case (int argc, char** argv) {

	/*
	Run the given case and output the results (Paraview)
	so they can be processed.
	
	Arguments:
		argc = The number of command line arguments
		argv = The array holding command line arguments
	
	Return:
		- 
	*/

	const char* ctrl_name = argv[1];

	struct Test_Info test_info = { .n_warn = 0, };

	// MSB: Read the test info. This will include the control file name and 
	// mesh level and order information
	struct Integration_Test_Info* int_test_info = constructor_Integration_Test_Info(ctrl_name);
	
	const int* p_ref  = int_test_info->p_ref,
	         * ml_ref = int_test_info->ml;

	struct Simulation* sim = NULL;
	
	const char type_rc = 'r';  // test for real and complex

	bool ignore_static = false;
	int ml_max = ml_ref[1];

	UNUSED(ml_max);
	UNUSED(argc);
	UNUSED(test_info);

	int ml_prev = ml_ref[0]-1,
    	p_prev  = p_ref[0]-1;

	// Run only the lowest mesh level case
	int ml = ml_ref[0];
	int p = p_ref[0];

	const int adapt_type = int_test_info->adapt_type;
	//const char*const ctrl_name_curr = set_file_name_curr(adapt_type,p,ml,false,ctrl_name);
	// MSB: Want to correct the ctrl file name because we want to add the mesh level and p
	const char*const ctrl_name_curr = set_file_name_curr(adapt_type,p,ml,true,ctrl_name);

	structor_simulation(&sim,'c',adapt_type,p,ml,p_prev,ml_prev,ctrl_name_curr,type_rc,ignore_static);
	solve_for_solution(sim);

	// Output the results in graphical format
	output_visualization(sim,VIS_GEOM_EDGES);
	output_visualization(sim,VIS_GEOM_VOLUMES);
	output_visualization(sim,VIS_NORMALS);
	output_visualization(sim,VIS_SOLUTION);

	output_error(sim);

	destructor_Integration_Test_Info(int_test_info);
}

// Level 0 ********************************************************************************************************** //


