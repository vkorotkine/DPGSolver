#include <stdlib.h>

/*
 *	Purpose:
 *		Swap arrays.
 *
 *	Comments:
 *		Run the version of array_swap.c in the testing branch to see timing comparisons between different
 *		implementations.
 *		Important conclusions from time testing:
 *			Pointer indexing is faster than array indexing when the index incrementing represents a significant fraction
 *			of the overall operation (e.g. in the case of swapping consecutive elements of large arrays); the difference
 *			between the array indexing and pointer implementations disappears when step sizes larger than 1 are taken
 *			even though the index incrementing still represents a large fraction of the operation (Investigate).
 *			Declaring loop indices/counters as 'unsigned' and 'register' variables can result in significant speed-ups.
 *
 *	Notation:
 *
 *	References:
 *		ToBeModified: Add in the reference page number (stephen kochan programming in c)
 *		              Chapter 10. A Slight Digression about Program Optimization
 *
 */

void array_swap_ui(register unsigned int *arr1, register unsigned int *arr2, const unsigned int NIn,
                   const unsigned int stepIn)
{
	register unsigned int N, step;
	register unsigned int tmp;

	tmp   = *arr1;
	*arr1 = *arr2;
	*arr2 = tmp;

	switch (NIn) {
	case 1:
		break;
	default:
		switch (stepIn) {
		case 1:
			for (N = NIn-1; N-- ; ) {
				arr1++;
				arr2++;

				tmp   = *arr1;
				*arr1 = *arr2;
				*arr2 = tmp;
			}
			break;
		default:
			for (N = NIn-1, step = stepIn; N-- ; ) {
				arr1 += step;
				arr2 += step;

				tmp   = *arr1;
				*arr1 = *arr2;
				*arr2 = tmp;
			}
			break;
		}
		break;
	}
}

/*
void array_swap_i(register int *arr1, register int *arr2, const int NIn, const int stepIn)
{
	register unsigned int N, step;
	register int tmp;

	tmp   = *arr1;
	*arr1 = *arr2;
	*arr2 = tmp;

	switch (NIn) {
	case 1:
		break;
	default:
		switch (stepIn) {
		case 1:
			for (N = NIn-1; N-- ; ) {
				arr1++;
				arr2++;

				tmp   = *arr1;
				*arr1 = *arr2;
				*arr2 = tmp;
			}
			break;
		default:
			for (N = NIn-1, step = stepIn; N-- ; ) {
				arr1 += step;
				arr2 += step;

				tmp   = *arr1;
				*arr1 = *arr2;
				*arr2 = tmp;
			}
			break;
		}
		break;
	}
}
*/

void array_swap_d(register double *arr1, register double *arr2, const int NIn, const int stepIn)
{
	register unsigned int N, step;
	register double tmp;

	tmp   = *arr1;
	*arr1 = *arr2;
	*arr2 = tmp;

	switch (NIn) {
	case 1:
		break;
	default:
		switch (stepIn) {
		case 1:
			for (N = NIn-1; N-- ; ) {
				arr1++;
				arr2++;

				tmp   = *arr1;
				*arr1 = *arr2;
				*arr2 = tmp;
			}
			break;
		default:
			for (N = NIn-1, step = stepIn; N-- ; ) {
				arr1 += step;
				arr2 += step;

				tmp   = *arr1;
				*arr1 = *arr2;
				*arr2 = tmp;
			}
			break;
		}
		break;
	}
}
