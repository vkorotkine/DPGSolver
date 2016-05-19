// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>

/*
 *	Purpose:
 *		Find the pointer to the first index whose array entries correspond to the value to be found.
 *
 *	Comments:
 *		The 'o' in the function name indicates that the input array is (o)rdered. In this case, a binary search is used.
 *		Note the usage of integer division (i.e. a = b/c is analogous to a = floor(b/c))
 *
 *	Notation:
 *		A   : Input array to looked in.
 *		val : Value of A to be found.
 *
 *	References:
 *
 */

void array_find_indexo_ui(const unsigned int LenA, const unsigned int *A, const unsigned int val,
                          unsigned int *IdxF, unsigned int *LenF)
{
	unsigned int IndL, IndR, IndM, Aval;

	*LenF = 0;

	IndL = 0;
	IndR = LenA - 1;
	IndM = (IndL + IndR)/2;
	Aval = A[IndM];

	if (Aval == val) *LenF = 1;
	while (Aval != val) {
		if      (val < Aval) IndR = IndM-1;
		else if (val > Aval) IndL = IndM+1;
		IndM = (IndL + IndR)/2;
		Aval = A[IndM];

// printf("Loop: %d %d %d %d %d %d\n",IndL,IndR,IndM,val,Aval,*LenF);

		if (Aval == val) {
			*LenF = 1;
			break;
		}

		if (IndL == IndR) {
			printf("Loop break: %d %d %d %d %d %d\n",IndL,IndR,IndM,val,Aval,*LenF);
			break;
		}
	}
	if (*LenF == 0)
		printf("Did not find any matches in array_find_indexo_ui"), exit(1);

	// Find how many times the value is repeated and the first entry
	IndL = IndM;
	IndR = IndM;
	while(IndL > 0      && A[IndL-1] == Aval) IndL--;
	while(IndR < LenA-1 && A[IndR+1] == Aval) IndR++;

	*LenF = IndR - IndL + 1;
	*IdxF = IndL;

// array_print_i(1,LenA,A);
// printf("%d %d %d %d %d %d %d\n",IndL,IndR,IndM,val,Aval,*IdxF,*LenF);

}
