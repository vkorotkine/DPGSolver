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
#include <string.h>
#include "petscsys.h"

#include "macros.h"
#include "definitions_core.h"

#include "test_base.h"
#include "test_integration.h"
#include "test_support.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

/** \test Performs integration testing for the fluxes (\ref test_integration_fluxes.c).
 *  \return 0 on success. */
int main
	(int nargc,  ///< Standard.
	 char** argv ///< Standard.
	)
{
	PetscInitialize(&nargc,&argv,PETSC_NULL,PETSC_NULL);

	assert_condition_message(nargc == 3,"Invalid number of input arguments");
	const char*const test_name = argv[1],
	          *const ctrl_name = argv[2];

	struct Integration_Test_Info* int_test_info = constructor_Integration_Test_Info(ctrl_name);

	const int p          = int_test_info->p_ref[0],
	          ml         = int_test_info->ml[0],
	          adapt_type = int_test_info->adapt_type;
	assert(adapt_type == ADAPT_0);

	const char*const ctrl_name_curr = set_file_name_curr(adapt_type,p,ml,false,ctrl_name);

	const char type_rc = 'r';

	struct Simulation* sim = NULL;
	structor_simulation(&sim,'c',adapt_type,p,ml,0,0,ctrl_name_curr,type_rc); // destructed

	struct Test_Info test_info = { .n_warn = 0, };
	if (strstr(test_name,"flux_advection"))
		EXIT_ADD_SUPPORT;
	else
		EXIT_ERROR("Invalid test name: %s\n",test_name);

	structor_simulation(&sim,'d',adapt_type,p,ml,0,0,NULL,type_rc);

	output_warning_count(&test_info);

	PetscFinalize();
	OUTPUT_SUCCESS;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
