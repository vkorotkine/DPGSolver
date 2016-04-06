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
 *		              Chapter 10. A Slight Digression about Program Optimization
 *
 */

/* Replace with fastest option found after testing */
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

/* Very slow compared to combined options */
void array_swap_d(double *arr1, double *arr2, const int NIn)
{
	unsigned int N;
	register double *tmp, *tmpP, *arr1P, *arr2P;

	tmp = malloc(NIn * sizeof *tmp); // free

	for (tmpP = tmp, arr1P = arr1, N = NIn; N--; *tmpP++ = *arr1P++)
		;
	for (arr1P = arr1, arr2P = arr2, N = NIn; N--; *arr1P++ = *arr2P++)
		;
	for (arr2P = arr2, tmpP = tmp, N = NIn; N--; *arr2P++ = *tmpP++)
		;

	free(tmp);
}

void array_swap_comb_d(double *arr1, double *arr2, const int NIn)
{
	register unsigned int N;
	double *tmp;
	register double *tmpP1, *tmpP2, *arr1P1, *arr1P2, *arr2P1, *arr2P2;

	tmp = malloc(NIn * sizeof *tmp); // free

	tmpP1 = tmp; tmpP2 = tmp;
	arr1P1 = arr1; arr1P2 = arr1;
	arr2P1 = arr2; arr2P2 = arr2;

	for (N = NIn; N-- ; ) {
		*tmpP1++  = *arr1P1++;
		*arr1P2++ = *arr2P1++;
		*arr2P2++ = *tmpP2++;
	}
	free(tmp);
}

void array_swap_comb_sep_d(double *arr1, double *arr2, const int NIn)
{
	register unsigned int N;
	double *tmp;
	register double *tmpP1, *tmpP2, *arr1P1, *arr1P2, *arr2P1, *arr2P2;

	tmp = malloc(NIn * sizeof *tmp); // free

	tmpP1 = tmp; tmpP2 = tmp;
	arr1P1 = arr1; arr1P2 = arr1;
	arr2P1 = arr2; arr2P2 = arr2;

	for (N = NIn; N-- ; ) {
		*tmpP1  = *arr1P1;
		tmpP1++, arr1P1++;
		*arr1P2 = *arr2P1;
		arr1P2++, arr2P1++;
		*arr2P2 = *tmpP2;
		arr2P2++, tmpP2++;
	}
	free(tmp);
}

void array_swap_comb_sep1_d(double *arr1, double *arr2, const int NIn)
{
	register unsigned int N;
	double *tmp;
	register double *tmpP1, *tmpP2, *arr1P1, *arr1P2, *arr2P1, *arr2P2;

	tmp = malloc(NIn * sizeof *tmp); // free

	tmpP1 = tmp; tmpP2 = tmp;
	arr1P1 = arr1; arr1P2 = arr1;
	arr2P1 = arr2; arr2P2 = arr2;

	for (N = NIn; N-- ; ) {
		*tmpP1  = *arr1P1;
		tmpP1   = tmpP1+1;
		arr1P1  = arr1P1+1;

		*arr1P2 = *arr2P1;
		arr1P2  = arr1P2+1;
		arr2P1  = arr2P1+1;

		*arr2P2 = *tmpP2;
		arr2P2  = arr2P2+1;
		tmpP2   = tmpP2+1;
	}
	free(tmp);
}

void array_swap_comb_NoPointer_d(register double *arr1, register double *arr2, const int NIn)
{
	register unsigned int i, N;
	register double *tmp;

	tmp = malloc(NIn * sizeof *tmp); // free

	for (i = 0, N = NIn; N-- ; i = i+1) {
		tmp[i]  = arr1[i];
		arr1[i] = arr2[i];
		arr2[i] = tmp[i];
	}
	free(tmp);
}

void array_swap_comb_NoPointer_minMem_d(register double *arr1, register double *arr2, const int NIn)
{
	register unsigned int i, N;
	register double tmp;

	for (i = 0, N = NIn; N-- ; i = i+1) {
		tmp     = arr1[i];
		arr1[i] = arr2[i];
		arr2[i] = tmp;
	}
}

void array_swap_comb_minMem_d(register double *arr1, register double *arr2, const int NIn)
{
	register unsigned int N;
	register double tmp;

	for (N = NIn; N-- ;) {
		tmp     = *arr1;
		*arr1++ = *arr2;
		*arr2++ = tmp;
	}
}

void array_swap_comb_minMem1_d(register double *arr1, register double *arr2, const int NIn)
{
	register unsigned int N;
	register double tmp;

	for (N = NIn; N-- ; ) {
		tmp   = *arr1;

		*arr1 = *arr2;
		arr1  = arr1+1;

		*arr2 = tmp;
		arr2  = arr2+1;
	}
}

void array_swap_comb_NoPointer_minMem_arb_d(register double *arr1, register double *arr2, const int NIn, const int stepIn)
{
	register unsigned int i, N, step;
	register double tmp;

	for (i = 0, N = NIn, step = stepIn; N-- ; i = i+step) {
		tmp     = arr1[i];
		arr1[i] = arr2[i]; *(arr1+i) = *(arr2+i)
		arr2[i] = tmp;
	}
}

void array_swap_comb_minMem_arb_d(register double *arr1, register double *arr2, const int NIn, const int stepIn)
{
	register unsigned int N, step;
	register double tmp;

	for (N = NIn, step = stepIn; N-- ; ) {
		tmp   = *arr1;
		*arr1 = *arr2;
		*arr2 = tmp;

		arr1 = arr1+step;
		arr2 = arr2+step;
	}
}

void array_swap_comb_minMem_arb_switch_d(register double *arr1, register double *arr2, const int NIn, const int stepIn)
{
	register unsigned int N, step;
	register double tmp;

	switch(stepIn) {
		case 1 :
			for (N = NIn, step = stepIn; N-- ; ) {
				tmp     = *arr1;
				*arr1++ = *arr2;
				*arr2++ = tmp;
			}
			break;
		default :
			for (N = NIn, step = stepIn; N-- ; ) {
				tmp   = *arr1;
				*arr1 = *arr2;
				*arr2 = tmp;

				arr1 = arr1+step;
				arr2 = arr2+step;
			}
			break;
	}
}




#ifdef DEBUG

#include <time.h>
#include <stdio.h>

int main(int argc, char **argv)
{
	int i, iMax, NE, Ndof;
	double *a, *b;

	NE = 2e3;
	Ndof = 20;

	a = malloc(NE*Ndof *sizeof *a);
	b = malloc(NE*Ndof *sizeof *b);

	for (i = 0, iMax = NE*Ndof; iMax--; i++)
		a[i] = (double) i;
	for (i = 0, iMax = NE*Ndof; iMax--; i++)
		b[i] = (double) iMax;


	clock_t ts, te;
	double tt;

/*
	ts = clock();
	array_swap_d(a,b,NE*Ndof);
	te = clock();
	tt = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("Time Separate                    : %10.6f\n",tt);

	ts = clock();
	array_swap_comb_d(a,b,NE*Ndof);
	te = clock();
	tt = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("Time Combined                    : %10.6f\n",tt);

	ts = clock();
	array_swap_comb_sep_d(a,b,NE*Ndof);
	te = clock();
	tt = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("Time Combined (SepPointer)       : %10.6f\n",tt);

	ts = clock();
	array_swap_comb_sep1_d(a,b,NE*Ndof);
	te = clock();
	tt = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("Time Combined (Sep1Pointer)      : %10.6f\n",tt);

	ts = clock();
	array_swap_comb_NoPointer_d(a,b,NE*Ndof);
	te = clock();
	tt = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("Time Combined (NoPointer)        : %10.6f\n",tt);
*/

	printf("Testing with consecutive element swapping:\n");

	ts = clock();
	array_swap_comb_NoPointer_minMem_d(a,b,NE*Ndof);
	array_swap_comb_NoPointer_minMem_d(a,b,NE*Ndof);
	te = clock();
	tt = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("Time Combined (NoPointer,minMem) : %10.6f\n",tt);

	ts = clock();
	array_swap_comb_minMem_d(a,b,NE*Ndof);
	array_swap_comb_minMem_d(a,b,NE*Ndof);
	te = clock();
	tt = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("Time Combined (minMem)           : %10.6f\n",tt);

	ts = clock();
	array_swap_comb_minMem1_d(a,b,NE*Ndof);
	array_swap_comb_minMem1_d(a,b,NE*Ndof);
	te = clock();
	tt = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("Time Combined (minMem1)          : %10.6f\n",tt);




	printf("\nTesting with sequential element swapping:\n");

	ts = clock();
	array_swap_comb_NoPointer_minMem_arb_d(a,b,NE*Ndof,1);
	array_swap_comb_NoPointer_minMem_arb_d(a,b,NE*Ndof,1);
	te = clock();
	tt = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("Time Combined (NoPointer,minMem,1)    : %10.6f\n",tt);

	ts = clock();
	array_swap_comb_minMem_arb_d(a,b,NE*Ndof,1);
	array_swap_comb_minMem_arb_d(a,b,NE*Ndof,1);
	te = clock();
	tt = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("Time Combined (minMem,1)              : %10.6f\n",tt);

	ts = clock();
	array_swap_comb_minMem_arb_switch_d(a,b,NE*Ndof,1);
	array_swap_comb_minMem_arb_switch_d(a,b,NE*Ndof,1);
	te = clock();
	tt = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("Time Combined (minMem,switch,1)       : %10.6f\n",tt);

	ts = clock();
	array_swap_comb_NoPointer_minMem_arb_d(a,b,NE,Ndof);
	array_swap_comb_NoPointer_minMem_arb_d(a,b,NE,Ndof);
	te = clock();
	tt = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("Time Combined (NoPointer,minMem,Ndof) : %10.6f\n",tt);

	ts = clock();
	array_swap_comb_minMem_arb_d(a,b,NE,Ndof);
	array_swap_comb_minMem_arb_d(a,b,NE,Ndof);
	te = clock();
	tt = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("Time Combined (minMem,Ndof)           : %10.6f\n",tt);

	ts = clock();
	array_swap_comb_minMem_arb_switch_d(a,b,NE,Ndof);
	array_swap_comb_minMem_arb_switch_d(a,b,NE,Ndof);
	te = clock();
	tt = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("Time Combined (minMem,switch,Ndof)    : %10.6f\n",tt);

/*	Conclusions:
 *		Removing the memory initialization for tmp results in large speed-up.
 *		The pointer version is generally faster or equivalent to the array indexing version, likely because of the extra
 *		computation of incrementing i at each step. This hypothesis is supported by the fact that the timing difference
 *		between the two implementations disappears when stepIn is increased.
 *		A very large difference in timing is also observed when loop counter/incrementing variables are not defined as
 *		register variables; the unsigned qualifier also results in slight speed-up.
 */

	return 0;
}

#endif
