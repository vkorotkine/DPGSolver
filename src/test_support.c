// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "test_support.h"

#include <stdio.h>

#include "Test.h"

void test_print2(unsigned int const pass, char const *const string)
{
	/*
	 *	Purpose:
	 *		Print pass/failure message for testing functions.
	 *
	 *	Comments:
	 *
	 *	Notation:
	 *
	 *	References:
	 */

	TestDB.Ntest++;

	printf("%-80s",string);
	if (pass) {
		TestDB.Npass++;
		printf("Pass\n");
	} else {
		printf("Fail --- Fail --- Fail\n");
	}
}

void test_print_warning(char const *const string)
{
	TestDB.Nwarnings++;

	printf("\n\n********************************************************************************************\n");
	printf("Warning: %s.\n",string);
	printf("********************************************************************************************\n\n\n");
}
