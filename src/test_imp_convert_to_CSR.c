// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>

#include "test.h"
#include "database.h"
#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Test correctness of implementation of convert_to_CSR.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void test_imp_convert_to_CSR(void)
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

	unsigned int NRows = 5, NCols = 5,
	             rowIndex55[6] = {0, 3, 5, 8, 11, 13},
	             columns55[13] = {0, 1, 3, 0, 1, 2, 3, 4, 0, 2, 3, 1, 4};
	double       values55[13]  = {1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 2.0, 7.0, 8.0, -5.0},
	             A[25]  = { 1.0, -1.0,  0.0, -3.0,  0.0,
	                       -2.0,  5.0,  0.0,  0.0,  0.0,
	                        0.0,  0.0,  4.0,  6.0,  4.0,
	                       -4.0,  0.0,  2.0,  7.0,  0.0,
	                        0.0,  8.0,  0.0,  0.0, -5.0};

	struct S_OpCSR *A_sp;

	A_sp = malloc(sizeof *A_sp); // free

	convert_to_CSR_d(NRows,NCols,A,A_sp);

	pass = 0;
	if (array_norm_diff_ui(6,A_sp->rowIndex,rowIndex55,"Inf") < EPS &&
	    array_norm_diff_ui(13,A_sp->columns,columns55,"Inf")  < EPS &&
	    array_norm_diff_d(13,A_sp->values,values55,"Inf")     < EPS)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("convert_to_CSR_d:                                ");
	test_print(pass);

	free(A_sp);
}
