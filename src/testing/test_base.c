// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "test_base.h"

#include <stdlib.h>
#include <stdio.h>

#include "Macros.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void test_increment_and_print (struct Test_Info*const test_info, const bool pass, const char*const test_name)
{
	printf("%-80s ... ",test_name);

	++test_info->n_test;

	if (pass) {
		++test_info->n_pass;
		printf("Pass\n");
	} else {
		printf("Fail --- Fail --- Fail\n");
	}
}

void test_print_warning (struct Test_Info*const test_info, const char*const warn_name)
{
	++test_info->n_warn;

	printf("\n********************************************************************************************\n");
	printf("Warning: %s.\n",warn_name);
	printf("********************************************************************************************\n\n");
}

void test_print_failure (struct Test_Info*const test_info, const char*const fail_name)
{
	UNUSED(test_info);
	printf("\n********************************************************************************************\n");
	printf("Test Failure: %s.\n",fail_name);
	printf("********************************************************************************************\n\n");
}

void output_test_info (struct Test_Info*const test_info)
{
	printf("\n\nRan %d test(s) in %.4f seconds.\n",
	       test_info->n_test,(test_info->te-test_info->ts)/(double)CLOCKS_PER_SEC);

	int n_fail = test_info->n_test - test_info->n_pass;
	if (n_fail) {
		printf("\n******** FAILED %d TEST(S) ********\n\n",n_fail);
	} else {
		printf("\nAll tests passed.\n\n");

		if (test_info->n_warn) {
			if (test_info->n_warn == 1)
				printf("1 warning was generated while running tests.\n");
			else
				printf("%d warnings were generated while running tests.\n",test_info->n_warn);
			printf("Scroll through test passing list and verify that all is OK.\n\n");
		}
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

