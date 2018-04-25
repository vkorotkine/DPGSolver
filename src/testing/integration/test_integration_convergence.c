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

#include "core.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

/** \test Performs integration testing for the solution convergence (\ref test_integration_convergence.c).
 *  \return 0 on success (when the solution converges at the expected rate).
 */
int main
	(int argc,   ///< Standard.
	 char** argv ///< Standard.
	)
{
	assert_condition_message((argc == 3) || (argc == 4),"Invalid number of input arguments");

	const char* petsc_options_name = set_petsc_options_name(argv[2]);
	PetscInitialize(&argc,&argv,petsc_options_name,PETSC_NULL);

	run_convergence_order_study(argc,argv,CONV_STUDY_SOLVE);

	PetscFinalize();
	OUTPUT_SUCCESS;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
