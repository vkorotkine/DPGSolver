// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "test_unit.h"

#include <stdlib.h>
#include <stdio.h>

#include "macros.h"

#include "test_unit_containers.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void run_tests_unit (struct Test_Info*const test_info)
{
	printf("\n\nRunning Unit Tests:\n");
	printf("----------------------------------------------------------------------------------------------------\n\n");

	test_unit_containers(test_info);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

