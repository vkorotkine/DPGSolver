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

#include "test_unit.h"

#include <stdlib.h>
#include <stdio.h>

#include "macros.h"

#include "test_base.h"
#include "test_unit_containers.h"
#include "test_unit_cubature.h"
#include "test_unit_bases.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void run_tests_unit (struct Test_Info*const test_info)
{
	printf("\n\nRunning Unit Tests:\n");
	printf("-------------------------------------------------------------------------------------------------\n\n");

	test_unit_containers(test_info);
	test_unit_cubature(test_info);
	test_print_warning(test_info,"Need to add tests for plotting nodes");
	test_unit_bases(test_info);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
