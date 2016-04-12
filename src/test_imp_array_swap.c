#include <stdlib.h>
#include <stdio.h>

#include "test.h"
#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Test correctness of implementation of array_swap.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */


void test_imp_array_swap(void)
{
	unsigned int pass;

	/*
	 *	array_swap_ui:
	 *
	 *		Input:
	 *
	 *			A = [ 1 2 3 4 5 6 7 8 9 ]
	 *			B = [ 9 8 7 6 5 4 3 2 1 ]
	 *
	 *		Expected Output:
	 *
	 *			NIn = 1 :
	 *				A = [ 9 2 3 4 5 6 7 8 9 ]
	 *				B = [ 1 8 7 6 5 4 3 2 1 ]
	 *
	 *			NIn = 9 (step = 1) :
	 *				A = [ 9 8 7 6 5 4 3 2 1 ]
	 *				B = [ 1 2 3 4 5 6 7 8 9 ]
	 *
	 *			NIn = 3 (step = 3) :
	 *				A = [ 9 2 3 6 5 6 3 8 9 ]
	 *				B = [ 1 8 7 4 5 4 7 2 1 ]
	 *
	 */

	unsigned int i;
	unsigned int A_ui[9], B_ui[9],
	             A_ui0[9] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 },
	             B_ui0[9] = { 9, 8, 7, 6, 5, 4, 3, 2, 1 };
	unsigned int A_ui11[9] = { 9, 2, 3, 4, 5, 6, 7, 8, 9 },
	             B_ui11[9] = { 1, 8, 7, 6, 5, 4, 3, 2, 1 },
	             A_ui91[9] = { 9, 8, 7, 6, 5, 4, 3, 2, 1 },
	             B_ui91[9] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 },
	             A_ui33[9] = { 9, 2, 3, 6, 5, 6, 3, 8, 9 },
	             B_ui33[9] = { 1, 8, 7, 4, 5, 4, 7, 2, 1 };

	// NIn = 1, step = 1
	for (i = 0; i < 9; i++) {
		A_ui[i] = A_ui0[i];
		B_ui[i] = B_ui0[i];
	}

	array_swap_ui(A_ui,B_ui,1,1);

	pass = 0;
	if (array_norm_diff_ui(9,A_ui,A_ui11,"Inf") < EPS &&
	    array_norm_diff_ui(9,B_ui,B_ui11,"Inf") < EPS)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("array_swap_ui (1,1):                             ");
	test_print(pass);

	// NIn = 9, step = 1
	for (i = 0; i < 9; i++) {
		A_ui[i] = A_ui0[i];
		B_ui[i] = B_ui0[i];
	}

	array_swap_ui(A_ui,B_ui,9,1);

	pass = 0;
	if (array_norm_diff_ui(9,A_ui,A_ui91,"Inf") < EPS &&
	    array_norm_diff_ui(9,B_ui,B_ui91,"Inf") < EPS)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("              (9,1):                             ");
	test_print(pass);

	// NIn = 9, step = 3
	for (i = 0; i < 9; i++) {
		A_ui[i] = A_ui0[i];
		B_ui[i] = B_ui0[i];
	}

	array_swap_ui(A_ui,B_ui,3,3);

	pass = 0;
	if (array_norm_diff_ui(9,A_ui,A_ui33,"Inf") < EPS &&
	    array_norm_diff_ui(9,B_ui,B_ui33,"Inf") < EPS)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("              (3,3):                             ");
	test_print(pass);

	/*
	 *	array_swap_d:
	 *
	 *		Same example as above cast to double.
	 */

	double A_d[9], A_d0[9], A_d11[9], A_d91[9], A_d33[9],
	       B_d[9], B_d0[9], B_d11[9], B_d91[9], B_d33[9];

	for (i = 0; i < 9; i++) {
		A_d0[i]  = (double) A_ui0[i];
		A_d11[i] = (double) A_ui11[i];
		A_d91[i] = (double) A_ui91[i];
		A_d33[i] = (double) A_ui33[i];

		B_d0[i]  = (double) B_ui0[i];
		B_d11[i] = (double) B_ui11[i];
		B_d91[i] = (double) B_ui91[i];
		B_d33[i] = (double) B_ui33[i];
	}

	// NIn = 1, step = 1
	for (i = 0; i < 9; i++) {
		A_d[i] = A_d0[i];
		B_d[i] = B_d0[i];
	}

	array_swap_d(A_d,B_d,1,1);

	pass = 0;
	if (array_norm_diff_d(9,A_d,A_d11,"Inf") < EPS &&
	    array_norm_diff_d(9,B_d,B_d11,"Inf") < EPS)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("array_swap_d (1,1):                              ");
	test_print(pass);

	// NIn = 9, step = 1
	for (i = 0; i < 9; i++) {
		A_d[i] = A_d0[i];
		B_d[i] = B_d0[i];
	}

	array_swap_d(A_d,B_d,9,1);

	pass = 0;
	if (array_norm_diff_d(9,A_d,A_d91,"Inf") < EPS &&
	    array_norm_diff_d(9,B_d,B_d91,"Inf") < EPS)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("             (9,1):                              ");
	test_print(pass);

	// NIn = 9, step = 3
	for (i = 0; i < 9; i++) {
		A_d[i] = A_d0[i];
		B_d[i] = B_d0[i];
	}

	array_swap_d(A_d,B_d,3,3);

	pass = 0;
	if (array_norm_diff_d(9,A_d,A_d33,"Inf") < EPS &&
	    array_norm_diff_d(9,B_d,B_d33,"Inf") < EPS)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("             (3,3):                              ");
	test_print(pass);
}
