// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>

#include "test.h"
#include "functions.h"

/*
 *	Purpose:
 *		Test correctness of implementation of factorial.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void test_imp_math_factorial(void)
{
	unsigned int pass;

	/*
	 *	factorial_ull:
	 *
	 *		Input:
	 *
	 *			A = 3
	 *			B = 7
	 *			C = 12
	 *
	 *		Expected Output:
	 *
	 *			A : 6
	 *			B : 5040
	 *			C : 479001600
	 */

	unsigned int A = 3, B = 7, C = 12;

	pass = 0;
	if (factorial_ull(A) == 6         &&
	    factorial_ull(B) == 5040      &&
	    factorial_ull(C) == 479001600)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("math_factorial_ull:                              ");
	test_print(pass);

}
