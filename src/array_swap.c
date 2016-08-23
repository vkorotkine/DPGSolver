// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "array_swap.h"

#include <stdlib.h>
#include <complex.h>
 
/*
 *	Purpose:
 *		Swap arrays.
 *
 *	Comments:
 *		Run 'test_speed_array_swap' when profiling and use fastest option for the code.
 *		Based on testing to date, it seems that alt2 is fastest in general, notable for the double swap with step > 1.
 *		(ToBeDeleted)
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

void array_swap_d(register double *arr1, register double *arr2, const unsigned int NIn, const unsigned int stepIn)
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

void array_swap_cmplx(register double complex *arr1, register double complex *arr2, const unsigned int NIn,
                      const unsigned int stepIn)
{
	register unsigned int N, step;
	register double complex tmp;

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
