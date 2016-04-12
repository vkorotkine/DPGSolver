#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "test.h"
#include "functions.h"

/*
 *	Purpose:
 *		Test speed of various implementations of array_swap.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void test_speed_array_swap(void)
{
	/*
	 *	array_swap_ui:
	 *
	 *		Input:
	 *
	 *			N, step
	 *			arr1 = rand(1,N);
	 *			arr2 = rand(1,N);
	 */

	unsigned int i, N, step, Ntest;
	unsigned int *arr1, *arr2, RAND_MAXui;

	double  tcode, ttest;
	clock_t ts, te;

	RAND_MAXui = (unsigned int) RAND_MAX;

	// N = 1, step = 1
	Ntest = 1e5;
	N = 1;
	step = 1;

	arr1 = malloc(N*step * sizeof *arr1); // free
	arr2 = malloc(N*step * sizeof *arr2); // free

	for (i = 0; i < N*step; i++) {
		arr1[i] = (unsigned int)rand()/RAND_MAXui;
		arr2[i] = (unsigned int)rand()/RAND_MAXui;
	}

	// Current implementation => Test fails if not fastest
	ts = clock();
	for (i = 0; i < Ntest; i++) {
		array_swap_ui(arr1,arr2,N,step);
		array_swap_ui(arr1,arr2,N,step);
	}
	te = clock();
	tcode = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("\nTime (N = %d, step = %d): %e\n",N,step,tcode);

	// Alternate implementation 1
	ts = clock();
	for (i = 0; i < Ntest; i++) {
		array_swap1_ui(arr1,arr2,N,step);
		array_swap1_ui(arr1,arr2,N,step);
	}
	te = clock();
	ttest = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("Time (N = %d, step = %d): %e\n",N,step,ttest);

	// Alternate implementation 2
	ts = clock();
	for (i = 0; i < Ntest; i++) {
		array_swap2_ui(arr1,arr2,N,step);
		array_swap2_ui(arr1,arr2,N,step);
	}
	te = clock();
	ttest = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("Time (N = %d, step = %d): %e\n",N,step,ttest);

	free(arr1);
	free(arr2);

	// N = 1e2, step = 1
	Ntest = 1e2;
	N = 1e2;
	step = 1;

	arr1 = malloc(N*step * sizeof *arr1); // free
	arr2 = malloc(N*step * sizeof *arr2); // free

	for (i = 0; i < N*step; i++) {
		arr1[i] = (unsigned int)rand()/RAND_MAXui;
		arr2[i] = (unsigned int)rand()/RAND_MAXui;
	}

	// Current implementation => Test fails if not fastest
	ts = clock();
	for (i = 0; i < Ntest; i++) {
		array_swap_ui(arr1,arr2,N,step);
		array_swap_ui(arr1,arr2,N,step);
	}
	te = clock();
	tcode = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("\nTime (N = %d, step = %d): %e\n",N,step,tcode);

	// Alternate implementation 1
	ts = clock();
	for (i = 0; i < Ntest; i++) {
		array_swap1_ui(arr1,arr2,N,step);
		array_swap1_ui(arr1,arr2,N,step);
	}
	te = clock();
	ttest = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("Time (N = %d, step = %d): %e\n",N,step,ttest);

	// Alternate implementation 2
	ts = clock();
	for (i = 0; i < Ntest; i++) {
		array_swap2_ui(arr1,arr2,N,step);
		array_swap2_ui(arr1,arr2,N,step);
	}
	te = clock();
	ttest = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("Time (N = %d, step = %d): %e\n",N,step,ttest);

	free(arr1);
	free(arr2);

	// N = 1e3, step = 10
	Ntest = 1e2;
	N = 1e3;
	step = 10;

	arr1 = malloc(N*step * sizeof *arr1); // free
	arr2 = malloc(N*step * sizeof *arr2); // free

	for (i = 0; i < N*step; i++) {
		arr1[i] = (unsigned int)rand()/RAND_MAXui;
		arr2[i] = (unsigned int)rand()/RAND_MAXui;
	}

	// Current implementation => Test fails if not fastest
	ts = clock();
	for (i = 0; i < Ntest; i++) {
		array_swap_ui(arr1,arr2,N,step);
		array_swap_ui(arr1,arr2,N,step);
	}
	te = clock();
	tcode = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("\nTime (N = %d, step = %d): %e\n",N,step,tcode);

	// Alternate implementation 1
	ts = clock();
	for (i = 0; i < Ntest; i++) {
		array_swap1_ui(arr1,arr2,N,step);
		array_swap1_ui(arr1,arr2,N,step);
	}
	te = clock();
	ttest = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("Time (N = %d, step = %d): %e\n",N,step,ttest);

	// Alternate implementation 2
	ts = clock();
	for (i = 0; i < Ntest; i++) {
		array_swap2_ui(arr1,arr2,N,step);
		array_swap2_ui(arr1,arr2,N,step);
	}
	te = clock();
	ttest = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("Time (N = %d, step = %d): %e\n",N,step,ttest);

	free(arr1);
	free(arr2);
}
