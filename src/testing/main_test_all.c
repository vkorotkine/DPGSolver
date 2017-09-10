// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/** \file
 */

#include <stdio.h>
#include <stdbool.h>

#include "test_base.h"
#include "test_unit.h"
#include "test_integration.h"

/** \brief Provides the main interface to run **all** currently supported tests.
 *	\return 0. */
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
