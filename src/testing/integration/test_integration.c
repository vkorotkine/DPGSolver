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
#include <string.h>
#include <mpi.h>

#include "macros.h"
#include "definitions_alloc.h"

#include "test_integration_mesh.h"
#include "test_integration_fe_init.h"
#include "test_integration_geometry.h"
#include "test_integration_euler.h"

#include "const_cast.h"
#include "file_processing.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void run_tests_integration (struct Test_Info*const test_info)
{
	MPI_Init(NULL,NULL);

	printf("\n\nRunning Integration Tests:\n");
	printf("-------------------------------------------------------------------------------------------------\n\n");

if (0) {
	test_integration_mesh(test_info,"curved_2d_mixed.msh");
	test_integration_mesh(test_info,"straight_2d_quad_periodic.msh");

	test_integration_fe_init(test_info,"extern_mesh/TEST_curved_2d_mixed");
	test_integration_fe_init(test_info,"extern_mesh/TEST_straight_2d_quad_periodic");

	test_integration_geometry(test_info,"extern_mesh/TEST_straight_2d_quad_periodic");
//	test_integration_geometry(test_info,"extern_mesh/TEST_curved_2d_mixed");
}

	test_integration_euler(test_info);

	MPI_Finalize();
}

struct Integration_Test_Info* constructor_Integration_Test_Info (const char*const ctrl_name_full)
{
	struct Integration_Test_Info* int_test_info = malloc(sizeof *int_test_info); // returned

	// Read information
	FILE *ctrl_file = fopen_checked(ctrl_name_full);

	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),ctrl_file)) {
		if (strstr(line,"ml_range_test")) read_skip_const_i_1(line,1,int_test_info->ml,2);
		if (strstr(line,"p_range_test"))  read_skip_const_i_1(line,1,int_test_info->p_ref,2);
	}
	fclose(ctrl_file);

	const_cast_i(&int_test_info->adapt_type,compute_adapt_type(int_test_info->p_ref,int_test_info->ml));

	return int_test_info;
}

void destructor_Integration_Test_Info (struct Integration_Test_Info* int_test_info)
{
	free(int_test_info);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

