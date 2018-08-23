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
#include "solve.h"

#include "core.h"

#include "matrix.h"
#include "matrix_print.h"
#include "multiarray.h"
#include "multiarray_math.h"
#include "multiarray_constructors.h"
#include "test_support_multiarray.h"
#include "test_support_matrix.h"

#include "optimization_case.h"
#include "adjoint.h"
#include "functionals.h"


#define L2_DIFFERENCE_MAX 5E-7


// Static function declarations ************************************************************************************* //

/** \brief Run the given case first (solve the initial flow over the 
 *	geometry) and then perform the test
 */
static void run_test (
	int argc, ///< The number of command line arguments
	char** argv  ///< The array holding command line arguments
	);

/** \brief This function will test the construction of the adjoint equation matrices.
 * 	the RHS will be only tested as the LHS is already validated in 
 *	test_integration_linearization.
 *
 *	\return True if the RHS is correct
 */
static bool test_adjoint (
	struct Optimization_Case *optimization_case,  ///< The optimization case data structure
	functional_fptr f,  ///< The functional to test sensitivites for
	functional_fptr_c f_c,  ///< The complex version of the functional being tested
	double tol  ///< Allowed tolerance for L2 norm of difference between complex step and finite difference computation
	);

// Interface functions ********************************************************************************************** //

/** \test Performs integration testing for the adjoint. Will compare the RHS to that found
 *	using finite differences.
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

	run_test(argc,argv);

	PetscFinalize();
	OUTPUT_SUCCESS;
}

// Static functions ************************************************************************************************* //

static void run_test (int argc, char** argv) {

	// =================================
	//         Solve the Flow
	// =================================

	const char* ctrl_name = argv[1];
	struct Test_Info test_info = { .n_warn = 0, };

	struct Integration_Test_Info* int_test_info = constructor_Integration_Test_Info(ctrl_name);
	
	const int* p_ref  = int_test_info->p_ref,
	         * ml_ref = int_test_info->ml;

	struct Simulation* sim = NULL;
	
	const char type_rc = 'r';
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
	const char*const ctrl_name_curr = set_file_name_curr(adapt_type,p,ml,true,ctrl_name);

	structor_simulation(&sim,'c',adapt_type,p,ml,p_prev,ml_prev,ctrl_name_curr,type_rc,ignore_static);
	solve_for_solution(sim);


	// =================================
	//        Test the Adjoint
	// =================================

	double tol = 1E-7;
	bool pass = true;

	struct Optimization_Case *optimization_case = constructor_Optimization_Case(sim);

	printf("\n======================================\n"); 
	printf("            functional_cl             \n");
	printf("======================================\n"); 
	if(!test_adjoint(optimization_case, functional_cl, functional_cl_c, tol))
		pass = false;
	
	printf("\n======================================\n"); 
	printf("          functional_cm_le            \n");
	printf("======================================\n"); 
	if(!test_adjoint(optimization_case, functional_cm_le, functional_cm_le_c, tol))
		pass = false;

	printf("\n======================================\n"); 
	printf("   functional_inverse_pressure_design \n");
	printf("======================================\n"); 
	if(!test_adjoint(optimization_case, functional_inverse_pressure_design, 
		functional_inverse_pressure_design_c, tol))
		pass = false;
	
	assert_condition(pass);

	// Clear allocated structures:
	destructor_Optimization_Case(optimization_case);
	destructor_Integration_Test_Info(int_test_info);
}


static bool test_adjoint (struct Optimization_Case *optimization_case, functional_fptr f, functional_fptr_c f_c, 
	double tol){

	// Setup and compute the sensitivities
	struct Adjoint_Data *adjoint_data_fd = constructor_Adjoint_Data(optimization_case);  // finite differences
	struct Adjoint_Data *adjoint_data_cs  = constructor_Adjoint_Data(optimization_case);  // complex step

	setup_adjoint_fd(adjoint_data_fd, optimization_case->sim, optimization_case->sim_c, f, f_c);
	setup_adjoint(adjoint_data_cs, optimization_case->sim, optimization_case->sim_c, f, f_c);
	

	// Check the differences
	bool pass = true;

	double LInf_dI_dW = diff_norm_Matrix_d(adjoint_data_cs->RHS, adjoint_data_fd->RHS, "LInf");
	bool pass_dI_dW = LInf_dI_dW < tol ? true : false;

	printf("LInf dI_dW : %e \n", LInf_dI_dW);

	if (!pass_dI_dW)
		pass = false;

	expect_condition(pass_dI_dW, "dI_dW");

	destructor_Adjoint_Data(adjoint_data_cs);
	destructor_Adjoint_Data(adjoint_data_fd);

	return pass;
}

// Level 0 ********************************************************************************************************** //





