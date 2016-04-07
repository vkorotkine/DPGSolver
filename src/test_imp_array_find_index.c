#include <stdlib.h>
#include <stdio.h>

#include "functions.h"

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

void test_imp_array_find_index(void)
{
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

	unsigned int pass = 0;

	unsigned int *A, Idx, Len;

	A = malloc(10 * sizeof *A); // free
	A[0] = 0; A[1] = 1; A[2] = 1; A[3] = 4; A[4] = 5;
	A[5] = 7; A[6] = 7; A[7] = 7; A[8] = 7; A[9] = 7;

	array_find_indexo_ui(10,A,7,&Idx,&Len);

	if (Idx == 5 && Len == 5)
		pass = 1;

	//     0         10        20        30        40        50
	printf("array_find_indexo_ui:                            ");
	test_print(pass);

	free(A);
}
