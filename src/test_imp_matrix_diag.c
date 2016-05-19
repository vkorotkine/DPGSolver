// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>

#include "test.h"
#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Test correctness of implementation of diag.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void test_imp_matrix_diag(void)
{
	unsigned int pass;

	/*
	 *	diag_d:
	 *
	 *		Input:
	 *
	 *			N = 4
	 *			x = [1.0 2.0 3.0 4.0]
	 *
	 *		Expected Output:
	 *
	 *			X = [1.0 0.0 0.0 0.0
	 *			     0.0 2.0 0.0 0.0
	 *			     0.0 0.0 3.0 0.0
	 *			     0.0 0.0 0.0 4.0]
	 */

	unsigned int N = 4;
	double *X,
	       x[4] = { 1.0, 2.0, 3.0, 4.0 },
	       X4[16] = {1.0, 0.0, 0.0 , 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 4.0};

	X = diag_d(x,N); // free

	pass = 0;
	if (array_norm_diff_d(16,X,X4,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("matrix_diag_d:                                   ");
	test_print(pass);

	free(X);
}
