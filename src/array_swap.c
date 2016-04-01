#include <stdlib.h>

/*
 *	Purpose:
 *		Swap arrays.
 *
 *	Comments:
 *		Note that when copying array elements in order, it is faster to copy as in the functions below as opposed to
 *		using the indices of each array.
 *		Test if it is just as fast to swap void arrays (resulting in only a single swap function). (ToBeDeleted)
 *
 *	Notation:
 *
 *  References:
 *		ToBeModified: Add in the reference page number (stephen kochan programming in c)
 *
 */

void array_swap_i(int *arr1, int *arr2, const int NIn)
{
	int N;
	register int *tmp, *tmpP, *arr1P, *arr2P;

	tmp = malloc(NIn * sizeof *tmp); // free

	for (tmpP = tmp, arr1P = arr1, N = NIn; N--; *tmpP++ = *arr1P++)
		;
	for (arr1P = arr1, arr2P = arr2, N = NIn; N--; *arr1P++ = *arr2P++)
		;
	for (arr2P = arr2, tmpP = tmp, N = NIn; N--; *arr2P++ = *tmpP++)
		;

	free(tmp);
}

void array_swap_d(double *arr1, double *arr2, const int NIn)
{
	int N;
	register double *tmp, *tmpP, *arr1P, *arr2P;

	tmp = malloc(NIn * sizeof *tmp); // free

	for (tmpP = tmp, arr1P = arr1, N = NIn; N--; *tmpP++ = *arr1P++)
		;
	for (arr1P = arr1, arr2P = arr2, N = NIn; N--; *arr1P++ = *arr2P++)
		;
	for (arr2P = arr2, tmpP = tmp, N = NIn; N--; *arr2P++ = *tmpP++)
		;

/*
 * Compare speed with code above. Also, potentially add switch statement to handle cases of small NIn using loop
 * unrolling.
 *
	unsigned int N;
	register double *tmpP1, *tmpP2, *arr1P1, *arr1P2, *arr2P1, *arr2P2;

	tmpP1 = tmp; tmpP2 = tmp;
	arr1P1 = arr1; arr1P2 = arr1;
	arr2P1 = arr2; arr2P2 = arr2;

	for (N = NIn; N-- ; ) {
		*arr1P1++ = *tmpP1++;
		*arr1P2++ = *arr2P1++;
		*arr2P2++ = *tmpP2++;
	}
*/

	free(tmp);
}
