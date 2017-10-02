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
#include "test_unit.h"
#include "test_integration.h"

/** \brief Provides the main interface to run **all** currently supported tests.
 *  \return 0 on successful exit. */
int main
	(int nargc,  ///< Standard.
	 char** argv ///< Standard.
	)
{
	struct {
		bool unit, integration;
	} run_tests = { .unit = true, .integration = true, };

	struct Test_Info test_info = { .nargc = nargc, .argv = argv, .n_test = 0, .n_pass = 0, .n_warn = 0, };

	test_info.t_int.equivalence_real_complex = true;
	test_info.t_int.equivalence_algorithms   = true;
	test_info.t_int.linearization            = true;
	test_info.t_int.conv_order               = true;

	printf("\nRunning All Tests:\n");
	test_info.ts = clock();

//	PetscInitialize(&nargc,&argv,PETSC_NULL,PETSC_NULL);

	if (run_tests.unit)
		run_tests_unit(&test_info);

	if (run_tests.integration)
		run_tests_integration(&test_info);

//	PetscFinalize();

	test_info.te = clock();

	output_test_info(&test_info);

	return 0;
}
