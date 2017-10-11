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
/**	\file
 */

#include "test_integration.h"

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#include "macros.h"
#include "test_integration_mesh.h"
#include "test_integration_fe_init.h"
#include "test_integration_geometry.h"
#include "test_integration_euler.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void run_tests_integration (struct Test_Info*const test_info)
{
	MPI_Init(NULL,NULL);

	printf("\n\nRunning Integration Tests:\n");
	printf("-------------------------------------------------------------------------------------------------\n\n");

	test_integration_mesh(test_info,"curved_2d_mixed.msh");
	test_integration_mesh(test_info,"straight_2d_quad_periodic.msh");

	test_integration_fe_init(test_info,"extern_mesh/TEST_curved_2d_mixed");
	test_integration_fe_init(test_info,"extern_mesh/TEST_straight_2d_quad_periodic");

	test_integration_geometry(test_info,"extern_mesh/TEST_straight_2d_quad_periodic");
	test_integration_geometry(test_info,"extern_mesh/TEST_curved_2d_mixed");

//	test_integration_euler(test_info);

	MPI_Finalize();
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

