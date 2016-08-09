// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>

#include "test.h"
#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Test correctness of implementation of inverse.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void test_imp_matrix_inverse(void)
{
	unsigned int pass;

	/*
	 *	inverse_d:
	 *
	 *		Input:
	 *
	 *			A = [0.81 0.91 0.28
	 *			     0.91 0.63 0.55
	 *			     0.13 0.10 0.96]
	 *
	 *		Expected Output:
	 *
	 *			AInv = [-1.949472564488963  2.998315752149631 -1.149188901693112
	 *			         2.844074106905415 -2.628135803563513  0.676181189610850
	 *			        -0.032266643028100 -0.132257778565730  1.126850456519812]
	 */

	unsigned int N = 3;
	double *I3, *Aroundtrip_c, *AInv_c,
	       A[9]    = { 0.81, 0.91, 0.28, 0.91, 0.63, 0.55, 0.13, 0.10, 0.96 },
	       AInv[9] = { -1.949472564488963,  2.998315752149631, -1.149188901693112,
	                    2.844074106905415, -2.628135803563513,  0.676181189610850,
	                   -0.032266643028100, -0.132257778565730,  1.126850456519812 };

	I3 = identity_d(N); // free

	AInv_c       = inverse_d(3,3,A,I3);      // free
	Aroundtrip_c = inverse_d(N,N,AInv_c,I3); // free

	pass = 0;
	if (array_norm_diff_d(9,AInv,AInv_c,"Inf")    < EPS &&
		array_norm_diff_d(9,A,Aroundtrip_c,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("matrix_inverse_d:                                ");
	test_print(pass);

	free(I3);
	free(AInv_c);
	free(Aroundtrip_c);
}
