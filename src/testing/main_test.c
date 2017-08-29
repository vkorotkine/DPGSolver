// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/** \file
 */

#include <stdio.h>
#include <time.h>
#include <stdbool.h>

#include "test_info.h"
#include "test_integration_euler.h"

/** \brief Provides the main interface to run **all** currently supported tests.
 *	\return 0. */
int main
	(int nargc,  ///< Standard.
	 char** argv ///< Standard.
	)
{
	struct Test_Info test_info = { .nargc = nargc, .argv = argv, .n_test = 0, .n_pass = 0, .n_warn = 0, };

	struct {
		bool unit, integration;
	} run_tests = { .unit = true, .integration = true, };

	test_info.t_int.equivalence_real_complex = true;
	test_info.t_int.equivalence_algorithms   = true;
	test_info.t_int.linearization            = true;
	test_info.t_int.conv_order               = true;

//	PetscInitialize(&nargc,&argv,PETSC_NULL,PETSC_NULL);

	printf("\n\nRunning Tests:\n\n\n");
	clock_t ts = clock();

	// Unit tests
	if (run_tests.unit) {
		;
	}

	// Integration tests
	if (run_tests.integration) {
		test_integration_euler(&test_info);
	}

//	PetscFinalize();

	clock_t te = clock();


/// \todo Move to external function.
	printf("\n\nRan %d test(s) in %.4f seconds.\n",test_info.n_test,(te-ts)/(double)CLOCKS_PER_SEC);

	int n_fail = test_info.n_test - test_info.n_pass;
	if (n_fail) {
		printf("\n******** FAILED %d TEST(S) ********\n\n",n_fail);
	} else {
		printf("\nAll tests passed.\n\n");

		if (test_info.n_warn)
			printf("%d warning(s) was/were generated while running tests.\n"
			       "Scroll through test passing list and verify that all is OK.\n\n",test_info.n_warn);
	}

	return 0;
}
