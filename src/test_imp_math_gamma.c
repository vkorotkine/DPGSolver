#include <stdlib.h>
#include <stdio.h>

#include "test.h"
#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Test correctness of implementation of gamma.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void test_imp_math_gamma(void)
{
	unsigned int pass;

	/*
	 *	gamma_ull:
	 *
	 *		Input:
	 *
	 *			A = 4
	 *			B = 8
	 *
	 *		Expected Output:
	 *
	 *			A : 6
	 *			B : 5040
	 */

	unsigned int A_ull = 4, B_ull = 8;

	pass = 0;
	if (gamma_ull(A_ull) == 6   &&
		gamma_ull(B_ull) == 5040)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("math_gamma_ull:                                  ");
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
