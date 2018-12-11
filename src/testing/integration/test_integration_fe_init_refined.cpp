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

#include <string.h>
#include "petscsys.h"

#include "test_base.h"
#include "test_integration.h"
#include "test_support.h"
#include "test_support_fe_init.h"

#include "macros.h"
#include "definitions_adaptation.h"

#include "simulation.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

/** \test Performs integration testing for the computational element initialization when the initial mesh is refined
 *        (\ref test_integration_fe_init_refined.cpp).
 *  \return 0 on success.
 *
 *  The test is otherwise identical to that of \ref test_integration_fe_init.cpp.
 */
int main
	(int argc,   ///< Standard.
	 char** argv ///< Standard.
	)
{
	PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

	assert_condition_message(argc == 2,"Invalid number of input arguments");
	const char* ctrl_name = argv[1];

	struct Test_Info test_info = { .n_warn = 0, };

	struct Simulation* sim = NULL;
	const char type_rc = 'r';
	structor_simulation(&sim,'c',ADAPT_H,0,0,-1,-1,ctrl_name,type_rc,false); // destructed

	const bool pass = compare_members_fe(&test_info,sim);

	structor_simulation(&sim,'d',ADAPT_0,0,0,-1,-1,NULL,type_rc,false);

	assert_condition(pass);
	output_warning_count(&test_info);

	PetscFinalize();
	OUTPUT_SUCCESS;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
