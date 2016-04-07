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

void test_imp_array_sort(void)
{
	/*
	 *	array_sort_ui:
	 *
	 *		Input:
	 *
	 *			A = [0 1 1 2 3 0 0 1 2 1
	 *			     1 3 3 5 1 0 1 0 4 3
	 *				 3 1 4 3 0 1 2 1 2 5]
	 *
	 *			Indices = [0 1 2 3 4 5 6 7 8 9]
	 *
	 *		Expected progression/output:
	 *
	 *			Unsorted:
	 *			0 1 1 2 3 0 0 1 2 1
	 *			1 3 3 5 1 0 1 0 4 3
	 *			3 1 4 3 0 1 2 1 2 5
	 *			-------------------
	 *			0 1 2 3 4 5 6 7 8 9
	 *
	 *			First row:
	 *			0 0 0 1 1 1 1 2 2 3
	 *			1 1 0 3 3 0 3 5 4 1
	 *			3 2 1 5 4 1 1 3 2 0
	 *			-------------------
	 *			0 6 5 9 2 7 1 3 8 4
	 *
	 *			Second row:
	 *			0-block               1-block               2-block               3-block (no change)
	 *			0 0 0 1 1 1 1 2 2 3   0 0 0 1 1 1 1 2 2 3   0 0 0 1 1 1 1 2 2 3   0 0 0 1 1 1 1 2 2 3
	 *			0 1 1 3 3 0 3 5 4 1   0 1 1 0 3 3 3 5 4 1   0 1 1 0 3 3 3 4 5 1   0 1 1 0 3 3 3 4 5 1
	 *			1 2 3 5 4 1 1 3 2 0   1 2 3 1 4 5 1 3 2 0   1 2 3 1 4 5 1 2 3 0   1 2 3 1 4 5 1 2 3 0
	 *			-------------------   -------------------   -------------------   -------------------
	 *			5 6 0 9 2 7 1 3 8 4   5 6 0 7 2 9 1 3 8 4   5 6 0 7 2 9 1 8 3 4   5 6 0 7 2 9 1 8 3 4
	 *
	 *			Third row:
	 *			00-block (no ch.)     01-block (no ch.)     10-block (no ch.)     13-block
	 *			0 0 0 1 1 1 1 2 2 3   0 0 0 1 1 1 1 2 2 3   0 0 0 1 1 1 1 2 2 3   0 0 0 1 1 1 1 2 2 3
	 *			0 1 1 0 3 3 3 5 4 1   0 1 1 0 3 3 3 5 4 1   0 1 1 0 3 3 3 4 5 1   0 1 1 0 3 3 3 4 5 1
	 *			1 2 3 1 4 5 1 2 3 0   1 2 3 1 4 5 1 2 3 0   1 2 3 1 4 5 1 2 3 0   1 2 3 1 1 4 5 2 3 0
	 *			-------------------   -------------------   -------------------   -------------------
	 *			5 6 0 7 2 9 1 8 3 4   5 6 0 7 2 9 1 8 3 4   5 6 0 7 2 9 1 8 3 4   5 6 0 7 1 2 9 8 3 4
	 *
	 *			24-block (no ch.)     25-block (no ch.)     31-block (no ch.)     Sorted!
	 *			0 0 0 1 1 1 1 2 2 3   0 0 0 1 1 1 1 2 2 3   0 0 0 1 1 1 1 2 2 3
	 *			0 1 1 0 3 3 3 4 5 1   0 1 1 0 3 3 3 4 5 1   0 1 1 0 3 3 3 4 5 1
	 *			1 2 3 1 1 4 5 2 3 0   1 2 3 1 1 4 5 2 3 0   1 2 3 1 1 4 5 2 3 0
	 *			-------------------   -------------------   -------------------
	 *			5 6 0 7 1 2 9 8 3 4   5 6 0 7 1 2 9 8 3 4   5 6 0 7 1 2 9 8 3 4
	 *
	 */

	unsigned int pass = 0;

	unsigned int i, j, rows, NRows = 3, NCols = 10;
	unsigned int *Ind, *Ind0, *Ind1, *Ind2, *Ind3;
	unsigned int *A, *A0, *A1, *A2, *A3;
	double errA[4], errInd[4];

	A  = malloc(NRows*NCols * sizeof *A);  // free
	A0 = malloc(NRows*NCols * sizeof *A0); // free
	A1 = malloc(NRows*NCols * sizeof *A1); // free
	A2 = malloc(NRows*NCols * sizeof *A2); // free
	A3 = malloc(NRows*NCols * sizeof *A3); // free

	A[0]  = 0; A[1]  = 1; A[2]  = 1; A[3]  = 2; A[4]  = 3; A[5]  = 0; A[6]  = 0; A[7]  = 1; A[8]  = 2; A[9]  = 1;
	A[10] = 1; A[11] = 3; A[12] = 3; A[13] = 5; A[14] = 1; A[15] = 0; A[16] = 1; A[17] = 0; A[18] = 4; A[19] = 3;
	A[20] = 3; A[21] = 1; A[22] = 4; A[23] = 3; A[24] = 0; A[25] = 1; A[26] = 2; A[27] = 1; A[28] = 2; A[29] = 5;

	A3[0]  = 0; A3[1]  = 0; A3[2]  = 0; A3[3]  = 1; A3[4]  = 1; A3[5]  = 1; A3[6]  = 1; A3[7]  = 2; A3[8]  = 2; A3[9]  = 3;
	A3[10] = 0; A3[11] = 1; A3[12] = 1; A3[13] = 0; A3[14] = 3; A3[15] = 3; A3[16] = 3; A3[17] = 4; A3[18] = 5; A3[19] = 1;
	A3[20] = 1; A3[21] = 2; A3[22] = 3; A3[23] = 1; A3[24] = 1; A3[25] = 4; A3[26] = 5; A3[27] = 2; A3[28] = 3; A3[29] = 0;

	for (i = 0; i < NRows; i++) {
	for (j = 0; j < NCols; j++) {
		A0[i*NCols+j] = A[i*NCols+j];
		if (i == 0) {
			A1[i*NCols+j] = A3[i*NCols+j];
			A2[i*NCols+j] = A3[i*NCols+j];
		} else if (i == 1) {
			A1[i*NCols+j] = A[i*NCols+j];
			A2[i*NCols+j] = A3[i*NCols+j];
		} else if (i == 2) {
			A1[i*NCols+j] = A[i*NCols+j];
			A2[i*NCols+j] = A[i*NCols+j];
		}
	}}

	array_print_ui(NRows,NCols,A,'R');
	array_print_ui(NRows,NCols,A0,'R');
	array_print_ui(NRows,NCols,A1,'R');
	array_print_ui(NRows,NCols,A2,'R');
	array_print_ui(NRows,NCols,A3,'R');

	Ind  = malloc(NCols * sizeof *Ind);  // free
	Ind0 = malloc(NCols * sizeof *Ind0); // free
	Ind1 = malloc(NCols * sizeof *Ind1); // free
	Ind2 = malloc(NCols * sizeof *Ind2); // free
	Ind3 = malloc(NCols * sizeof *Ind3); // free

	for (i = 0; i < NCols; i++) {
		Ind[i] = i;
		Ind0[i] = Ind[i];
	}

	Ind1[0] = 0; Ind1[1] = 6; Ind1[2] = 5; Ind1[3] = 9; Ind1[4] = 2;
	Ind1[5] = 7; Ind1[6] = 1; Ind1[7] = 3; Ind1[8] = 8; Ind1[9] = 4;

	Ind2[0] = 5; Ind2[1] = 6; Ind2[2] = 0; Ind2[3] = 7; Ind2[4] = 2;
	Ind2[5] = 9; Ind2[6] = 1; Ind2[7] = 8; Ind2[8] = 3; Ind2[9] = 4;

	Ind3[0] = 5; Ind3[1] = 6; Ind3[2] = 0; Ind3[3] = 7; Ind3[4] = 1;
	Ind3[5] = 2; Ind3[6] = 9; Ind3[7] = 8; Ind3[8] = 3; Ind3[9] = 4;

	rows = 0;
	errA[rows]   = 0.0;
	errInd[rows] = 0.0;
	for (i = 0; i < NRows; i++) {
	for (j = 0; j < NCols; j++) {
//		errA[rows] += 
		// PROBABLY WRITE AN ARRAY NORM DIFF FUNCTION (i.e. norm of one array minus another)
	}}




	array_sort_ui(NRows,NCols,A,Ind,'R','N');

	array_print_ui(NRows,NCols,A,'R');


/*
	errorInf = array_norm_d(3,A,"Inf")-3.0;
	errorL1  = array_norm_d(3,A,"L1")-6.0;
	errorL2  = array_norm_d(3,A,"L2")-3.741657386773941;

	if (array_norm_d(1,&errorInf,"Inf") < EPS &&
		array_norm_d(1,&errorL1,"Inf")  < EPS &&
		array_norm_d(1,&errorL2,"Inf")  < EPS)
		pass = 1;
*/

	//     0         10        20        30        40        50
	printf("array_sort_ui:                                   ");
	test_print(pass);

	free(A);
	free(A0);
	free(A1);
	free(A2);
	free(A3);

	free(Ind);
	free(Ind0);
	free(Ind1);
	free(Ind2);
	free(Ind3);
}
