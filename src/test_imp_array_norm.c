#include <stdlib.h>
#include <stdio.h>

#include "functions.h"
#include "parameters.h"

/*
 *	Purpose:
 *		Test correctness of implementation of array_norm.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
*/

void test_imp_array_norm(void)
{
	/*
	 *	array_norm_d:
	 *
	 *		Input:
	 *
	 *			A = [1.0 2.0 3.0]
	 *
	 *		Expected output:
	 *
	 *			Inf : 3.0
	 *			L1  : 6.0
	 *			L2  : 3.741657386773941
	 */

	unsigned int pass;

	double *A, errorInf, errorL1, errorL2;

	A = malloc(3 * sizeof *A); // free
	A[0] = 1.0; A[1] = 2.0; A[2] = 3.0;

	errorInf = array_norm_d(3,A,"Inf")-3.0;
	errorL1  = array_norm_d(3,A,"L1")-6.0;
	errorL2  = array_norm_d(3,A,"L2")-3.741657386773941;

	pass = 0;
	if (array_norm_d(1,&errorInf,"Inf") < EPS)
		pass = 1;

	//     0         10        20        30        40        50
	printf("array_norm_d (Inf):                              ");
	test_print(pass);

	pass = 0;
	if (array_norm_d(1,&errorL1,"Inf")  < EPS)
		pass = 1;

	//     0         10        20        30        40        50
	printf("             (L1) :                              ");
	test_print(pass);

	pass = 0;
	if (array_norm_d(1,&errorL2,"Inf")  < EPS)
		pass = 1;

	//     0         10        20        30        40        50
	printf("             (L2) :                              ");
	test_print(pass);


	free(A);
}
