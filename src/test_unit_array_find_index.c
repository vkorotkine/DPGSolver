// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "test_unit_array_find_index.h"

#include <stdlib.h>
#include <stdio.h>

#include "Test.h"

#include "test_support.h"
#include "array_find_index.h"

/*
 *	Purpose:
 *		Test correctness of implementation of array_find_index.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void test_unit_array_find_index(void)
{
	unsigned int pass;

	/*
	 *	array_find_indexo_ui:
	 *
	 *		Input:
	 *
	 *			A = [0 1 1 4 5 7 7 7 7 7]
	 *			val = 7
	 *
	 *		Expected output:
	 *
	 *			Idx = 5
	 *			Len = 5
	 */

	unsigned int A[10] = { 0, 1, 1, 4, 5, 7, 7, 7, 7, 7 };
	unsigned int Idx, Len;

	array_find_indexo_ui(10,A,7,&Idx,&Len);

	pass = 0;
	if (Idx == 5 && Len == 5)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("array_find_indexo_ui:                            ");
	test_print(pass);
}
