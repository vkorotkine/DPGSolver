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
/**	\file
 */

#include "test_base.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "macros.h"
#include "definitions_alloc.h"

// Static function declarations ************************************************************************************* //

/** \brief Allocates memory and sets the path to the data file.
 *  \return See brief. */
static char* constructor_file_name_base ();

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
	printf("Warning: %s\n",warn_name);
	printf("********************************************************************************************\n\n");
}

void test_print_failure (struct Test_Info*const test_info, const char*const fail_name)
{
	UNUSED(test_info);
	printf("\n********************************************************************************************\n");
	printf("Test Failure: %s\n",fail_name);
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

char* constructor_file_name_unit (const char*const file_name_spec)
{
	char*const file_name = constructor_file_name_base();

	strcat(file_name,"/");
	strcat(file_name,file_name_spec);
	strcat(file_name,".data");

	return file_name;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static char* constructor_file_name_base ()
{
	char*const file_name_base = malloc(STRLEN_MAX * sizeof *file_name_base); // returned

	strcpy(file_name_base,"../testing/unit");

	return file_name_base;
}
