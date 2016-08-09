// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "test_support.h"

void test_print(const unsigned int pass)
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

	if (pass)
		printf("Pass\n");
	else
		printf("Fail --- Fail --- Fail\n");
}
