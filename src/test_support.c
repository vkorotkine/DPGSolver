// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "test_support.h"

#include <stdio.h>

#include "Test.h"

void test_print2(const unsigned int pass, const char *string)
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
