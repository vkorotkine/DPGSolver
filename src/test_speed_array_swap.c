// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "test_speed_array_swap.h"

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
	printf("\ntest_speed_array_swap:\n\n");

	unsigned int i, N, step, Ntest;
	double  tcode, ttest;
	clock_t ts, te;

	/*
	 *	array_swap_ui:
	 *
	 *		Input:
	 *
	 *			N, step
	 *			arr1 = rand(1,N);
	 *			arr2 = rand(1,N);
	 */

	unsigned int *arr1ui, *arr2ui, RAND_MAXui;

	RAND_MAXui = (unsigned int) RAND_MAX;

	// N = 1, step = 1
	Ntest = 1e5;
	N = 1;
	step = 1;

	arr1ui = malloc(N*step * sizeof *arr1ui); // free
	arr2ui = malloc(N*step * sizeof *arr2ui); // free

	for (i = 0; i < N*step; i++) {
		arr1ui[i] = (unsigned int)rand()/RAND_MAXui;
		arr2ui[i] = (unsigned int)rand()/RAND_MAXui;
	}

	// Current implementation
	ts = clock();
	for (i = 0; i < Ntest; i++) {
		array_swap_ui(arr1ui,arr2ui,N,step);
		array_swap_ui(arr1ui,arr2ui,N,step);
	}
	te = clock();
	tcode = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("\narray_swap_ui, code (N = %d, step = %d): %e\n",N,step,tcode);

	// Alternate implementation 1
	ts = clock();
	for (i = 0; i < Ntest; i++) {
		array_swap1_ui(arr1ui,arr2ui,N,step);
		array_swap1_ui(arr1ui,arr2ui,N,step);
	}
	te = clock();
	ttest = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("               alt1 (N = %d, step = %d): %e\n",N,step,ttest);

	// Alternate implementation 2
	ts = clock();
	for (i = 0; i < Ntest; i++) {
		array_swap2_ui(arr1ui,arr2ui,N,step);
		array_swap2_ui(arr1ui,arr2ui,N,step);
	}
	te = clock();
	ttest = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("               alt2 (N = %d, step = %d): %e\n",N,step,ttest);

	free(arr1ui);
	free(arr2ui);

	// N = 1e2, step = 1
	Ntest = 1e2;
	N = 1e2;
	step = 1;

	arr1ui = malloc(N*step * sizeof *arr1ui); // free
	arr2ui = malloc(N*step * sizeof *arr2ui); // free

	for (i = 0; i < N*step; i++) {
		arr1ui[i] = (unsigned int)rand()/RAND_MAXui;
		arr2ui[i] = (unsigned int)rand()/RAND_MAXui;
	}

	// Current implementation
	ts = clock();
	for (i = 0; i < Ntest; i++) {
		array_swap_ui(arr1ui,arr2ui,N,step);
		array_swap_ui(arr1ui,arr2ui,N,step);
	}
	te = clock();
	tcode = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("\n               code (N = %d, step = %d): %e\n",N,step,tcode);

	// Alternate implementation 1
	ts = clock();
	for (i = 0; i < Ntest; i++) {
		array_swap1_ui(arr1ui,arr2ui,N,step);
		array_swap1_ui(arr1ui,arr2ui,N,step);
	}
	te = clock();
	ttest = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("               alt1 (N = %d, step = %d): %e\n",N,step,ttest);

	// Alternate implementation 2
	ts = clock();
	for (i = 0; i < Ntest; i++) {
		array_swap2_ui(arr1ui,arr2ui,N,step);
		array_swap2_ui(arr1ui,arr2ui,N,step);
	}
	te = clock();
	ttest = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("               alt2 (N = %d, step = %d): %e\n",N,step,ttest);

	free(arr1ui);
	free(arr2ui);

	// N = 1e3, step = 10
	Ntest = 1e2;
	N = 1e3;
	step = 10;

	arr1ui = malloc(N*step * sizeof *arr1ui); // free
	arr2ui = malloc(N*step * sizeof *arr2ui); // free

	for (i = 0; i < N*step; i++) {
		arr1ui[i] = (unsigned int)rand()/RAND_MAXui;
		arr2ui[i] = (unsigned int)rand()/RAND_MAXui;
	}

	// Current implementation
	ts = clock();
	for (i = 0; i < Ntest; i++) {
		array_swap_ui(arr1ui,arr2ui,N,step);
		array_swap_ui(arr1ui,arr2ui,N,step);
	}
	te = clock();
	tcode = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("\n               code (N = %d, step = %d): %e\n",N,step,tcode);

	// Alternate implementation 1
	ts = clock();
	for (i = 0; i < Ntest; i++) {
		array_swap1_ui(arr1ui,arr2ui,N,step);
		array_swap1_ui(arr1ui,arr2ui,N,step);
	}
	te = clock();
	ttest = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("               alt1 (N = %d, step = %d): %e\n",N,step,ttest);

	// Alternate implementation 2
	ts = clock();
	for (i = 0; i < Ntest; i++) {
		array_swap2_ui(arr1ui,arr2ui,N,step);
		array_swap2_ui(arr1ui,arr2ui,N,step);
	}
	te = clock();
	ttest = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("               alt2 (N = %d, step = %d): %e\n",N,step,ttest);

	free(arr1ui);
	free(arr2ui);

	/*
	 *	array_swap_d:
	 *
	 *		Input:
	 *
	 *			Same example as above cast to double.
	 */

	double *arr1d, *arr2d, RAND_MAXd;

	RAND_MAXd = (double) RAND_MAX;

	// N = 1, step = 1
	Ntest = 1e5;
	N = 1;
	step = 1;

	arr1d = malloc(N*step * sizeof *arr1d); // free
	arr2d = malloc(N*step * sizeof *arr2d); // free

	for (i = 0; i < N*step; i++) {
		arr1d[i] = (double)rand()/RAND_MAXd;
		arr2d[i] = (double)rand()/RAND_MAXd;
	}

	// Current implementation
	ts = clock();
	for (i = 0; i < Ntest; i++) {
		array_swap_d(arr1d,arr2d,N,step);
		array_swap_d(arr1d,arr2d,N,step);
	}
	te = clock();
	tcode = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("\narray_swap_d, code (N = %d, step = %d): %e\n",N,step,tcode);

	// Alternate implementation 1
	ts = clock();
	for (i = 0; i < Ntest; i++) {
		array_swap1_d(arr1d,arr2d,N,step);
		array_swap1_d(arr1d,arr2d,N,step);
	}
	te = clock();
	ttest = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("              alt1 (N = %d, step = %d): %e\n",N,step,ttest);

	// Alternate implementation 2
	ts = clock();
	for (i = 0; i < Ntest; i++) {
		array_swap2_d(arr1d,arr2d,N,step);
		array_swap2_d(arr1d,arr2d,N,step);
	}
	te = clock();
	ttest = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("              alt2 (N = %d, step = %d): %e\n",N,step,ttest);

	free(arr1d);
	free(arr2d);

	// N = 1e2, step = 1
	Ntest = 1e3;
	N = 1e2;
	step = 1;

	arr1d = malloc(N*step * sizeof *arr1d); // free
	arr2d = malloc(N*step * sizeof *arr2d); // free

	for (i = 0; i < N*step; i++) {
		arr1d[i] = (double)rand()/RAND_MAXd;
		arr2d[i] = (double)rand()/RAND_MAXd;
	}

	// Current implementation
	ts = clock();
	for (i = 0; i < Ntest; i++) {
		array_swap_d(arr1d,arr2d,N,step);
		array_swap_d(arr1d,arr2d,N,step);
	}
	te = clock();
	tcode = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("\n              code (N = %d, step = %d): %e\n",N,step,tcode);

	// Alternate implementation 1
	ts = clock();
	for (i = 0; i < Ntest; i++) {
		array_swap1_d(arr1d,arr2d,N,step);
		array_swap1_d(arr1d,arr2d,N,step);
	}
	te = clock();
	ttest = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("              alt1 (N = %d, step = %d): %e\n",N,step,ttest);

	// Alternate implementation 2
	ts = clock();
	for (i = 0; i < Ntest; i++) {
		array_swap2_d(arr1d,arr2d,N,step);
		array_swap2_d(arr1d,arr2d,N,step);
	}
	te = clock();
	ttest = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("              alt2 (N = %d, step = %d): %e\n",N,step,ttest);

	free(arr1d);
	free(arr2d);

	// N = 1e3, step = 10
	Ntest = 1e2;
	N = 1e3;
	step = 10;

	arr1d = malloc(N*step * sizeof *arr1d); // free
	arr2d = malloc(N*step * sizeof *arr2d); // free

	for (i = 0; i < N*step; i++) {
		arr1d[i] = (double)rand()/RAND_MAXd;
		arr2d[i] = (double)rand()/RAND_MAXd;
	}

	// Current implementation
	ts = clock();
	for (i = 0; i < Ntest; i++) {
		array_swap_d(arr1d,arr2d,N,step);
		array_swap_d(arr1d,arr2d,N,step);
	}
	te = clock();
	tcode = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("\n              code (N = %d, step = %d): %e\n",N,step,tcode);

	// Alternate implementation 1
	ts = clock();
	for (i = 0; i < Ntest; i++) {
		array_swap1_d(arr1d,arr2d,N,step);
		array_swap1_d(arr1d,arr2d,N,step);
	}
	te = clock();
	ttest = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("              alt1 (N = %d, step = %d): %e\n",N,step,ttest);

	// Alternate implementation 2
	ts = clock();
	for (i = 0; i < Ntest; i++) {
		array_swap2_d(arr1d,arr2d,N,step);
		array_swap2_d(arr1d,arr2d,N,step);
	}
	te = clock();
	ttest = (te-ts)/(double)CLOCKS_PER_SEC;
	printf("              alt2 (N = %d, step = %d): %e\n",N,step,ttest);

	free(arr1d);
	free(arr2d);
}
