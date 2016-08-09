// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>

#include "test.h"
#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Test correctness of implementation of identity.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void test_imp_matrix_identity(void)
{
	unsigned int pass;

	/*
	 *	identity_d:
	 *
	 *		Input:
	 *
	 *			N = 4
	 *
	 *		Expected Output:
	 *
	 *			I = [1.0 0.0 0.0 0.0
	 *			     0.0 1.0 0.0 0.0
	 *			     0.0 0.0 1.0 0.0
	 *			     0.0 0.0 0.0 1.0]
	 */

	unsigned int N = 4;
	double *Imat, // ToBeModified (name)
	       I4[16] = {1.0, 0.0, 0.0 , 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};

	Imat = identity_d(N); // free

	pass = 0;
	if (array_norm_diff_d(16,Imat,I4,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("matrix_identity_d:                               ");
	test_print(pass);

	free(Imat);
}
