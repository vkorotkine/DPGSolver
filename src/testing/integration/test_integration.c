// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "test_integration.h"

#include <stdlib.h>
#include <stdio.h>

#include "Macros.h"
#include "test_integration_mesh.h"
#include "test_integration_euler.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void run_tests_integration (struct Test_Info*const test_info)
{
	test_integration_mesh(test_info,"curved_2d_mixed.msh");
	test_integration_mesh(test_info,"straight_2d_quad_periodic.msh");

//	test_integration_euler(test_info);
	return;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

