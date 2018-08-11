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
#include <stdlib.h>
#include <stdio.h>

#include "petscsys.h"
#include "parmetis.h"

#include "simulation.h"
#include "test_base.h"
#include "test_integration.h"
#include "test_integration_convergence_support.h"
//#include "test_support.h"
//#include "test_support_matrix.h"
//#include "test_support_vector.h"
//
//#include "macros.h"
//#include "definitions_core.h"
//#include "definitions_tol.h"
//
//#include "matrix.h"
//#include "vector.h"
//
//#include "math_functions.h"

#include "macros.h"
#include "core.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

/** \test Performs integration testing for the MPI library.
 *  \return 0 on success.
 *
 *  TO DO:
 *		Figure out which tests we want to perform.
 */
int main
	(int argc,   ///< Standard.
	 char** argv ///< Standard.
	)
{
	assert_condition_message(argc == 2,"Invalid number of input arguments");

	const char* ctrl_name = argv[1];
	const char* petsc_options_name = set_petsc_options_name(argv[2]);

	PetscInitialize(&argc, &argv, petsc_options_name, PETSC_NULL);

	struct Integration_Test_Info* int_test_info = constructor_Integration_Test_Info(ctrl_name);

	PetscFinalize();
	OUTPUT_SUCCESS;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

