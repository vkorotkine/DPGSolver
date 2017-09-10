// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "test_integration.h"

#include <stdlib.h>
#include <stdio.h>

#include "macros.h"
#include "test_integration_mesh.h"
#include "test_integration_fe_init.h"
#include "test_integration_geometry.h"
#include "test_integration_euler.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void run_tests_integration (struct Test_Info*const test_info)
{
	printf("\n\nRunning Integration Tests:\n");
	printf("----------------------------------------------------------------------------------------------------\n\n");

	test_integration_mesh(test_info,"curved_2d_mixed.msh");
	test_integration_mesh(test_info,"straight_2d_quad_periodic.msh");

	test_integration_fe_init(test_info,"extern_mesh/TEST_curved_2d_mixed");
	test_integration_fe_init(test_info,"extern_mesh/TEST_straight_2d_quad_periodic");

	test_integration_geometry(test_info,"extern_mesh/TEST_straight_2d_quad_periodic");

//	test_integration_euler(test_info);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

