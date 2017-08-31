// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/** \file
 */

#include <stdio.h>
#include <stdbool.h>

#include "Macros.h"

#include "test_base.h"
#include "test_integration_euler.h"

// Static function declarations ************************************************************************************* //

/// \brief Call unit test functions.
static void run_tests_unit
	(struct Test_Info*const test_info ///< \ref Test_Info.
	);

/// \brief Call integration test functions.
static void run_tests_integration
	(struct Test_Info*const test_info ///< \ref Test_Info.
	);

// Interface functions ********************************************************************************************** //

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

	printf("\n\nRunning Tests:\n\n\n");
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

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void run_tests_unit (struct Test_Info*const test_info)
{
	UNUSED(test_info);
	return;
}

static void run_tests_integration (struct Test_Info*const test_info)
{
	test_integration_euler(test_info);
	return;
}
