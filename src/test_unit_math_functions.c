// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "test_unit_math_functions.h"

#include <stdlib.h>
#include <stdio.h>

#include "Parameters.h"
#include "Test.h"

#include "test_support.h"
#include "array_norm.h"
#include "math_functions.h"

/*
 *	Purpose:
 *		Test correctness of implementation of math functions.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void test_unit_math_factorial(void)
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

void test_unit_math_gamma(void)
{
	unsigned int pass;

	/*
	 *	gamma_d:
	 *
	 *		Input:
	 *
	 *			A = 4
	 *			B = 8
	 *			C = 16
	 *
	 *		Expected Output:
	 *
	 *			A : 6
	 *			B : 5040
	 *			C : 1307674368000
	 */

	unsigned int A_ull = 4, B_ull = 8, C_ull = 16;

	pass = 0;
	if (gamma_d(A_ull) == 6.0    &&
		gamma_d(B_ull) == 5040.0 &&
		gamma_d(C_ull) == 1307674368000.0)
		    pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("math_gamma_d:                                    ");
	test_print(pass);

	/*
	 *	gamma_d:
	 *
	 *		Input:
	 *
	 *			A = 4.0
	 *			B = 2.72
	 *			C = 3.14
	 *
	 *		Expected Output:
	 *
	 *			A : 6.0
	 *			B : 1.569638593925876
	 *			C : 2.284480633817801
	 */

	double A_d = 4.0, B_d = 2.72, C_d = 3.14;
	double errA_d, errB_d, errC_d;

	errA_d = gamma_d(A_d) - 6.0;
	errB_d = gamma_d(B_d) - 1.569638593925876;
	errC_d = gamma_d(C_d) - 2.284480633817801;

	pass = 0;
	if (array_norm_d(1,&errA_d,"Inf") < EPS    &&
	    array_norm_d(1,&errB_d,"Inf") < EPS*10 &&
	    array_norm_d(1,&errC_d,"Inf") < EPS*10)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("math_gamma_d:                                    ");
	test_print(pass);
}
