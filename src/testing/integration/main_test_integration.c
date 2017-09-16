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

#include <stdio.h>
#include <stdbool.h>

#include "test_base.h"
#include "tests_integration.h"

/** \brief Provides the main interface to run **all** currently supported integration tests.
 *	\return 0. */
int main
	(int nargc,  ///< Standard.
	 char** argv ///< Standard.
	)
{
	struct Test_Info test_info = { .nargc = nargc, .argv = argv, .n_test = 0, .n_pass = 0, .n_warn = 0, };

	printf("\n\nRunning Unit Tests:\n\n\n");
	test_info.ts = clock();

	run_tests_integration(&test_info);

	test_info.te = clock();

	output_test_info(&test_info);

	return 0;
}
