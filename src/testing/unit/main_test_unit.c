// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/** \file
 */

#include <stdio.h>
#include <stdbool.h>

#include "test_base.h"
#include "tests_unit.h"

/** \brief Provides the main interface to run **all** currently supported unit tests.
 *	\return 0. */
int main
	(int nargc,  ///< Standard.
	 char** argv ///< Standard.
	)
{
	struct Test_Info test_info = { .nargc = nargc, .argv = argv, .n_test = 0, .n_pass = 0, .n_warn = 0, };

	printf("\n\nRunning Unit Tests:\n\n\n");
	test_info.ts = clock();

	run_tests_unit(&test_info);

	test_info.te = clock();

	output_test_info(&test_info);

	return 0;
}
