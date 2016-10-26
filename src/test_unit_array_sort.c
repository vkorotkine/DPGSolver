// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "test_unit_array_sort.h"

#include <stdlib.h>
#include <stdio.h>

#include "mkl.h"

#include "Parameters.h"
#include "Test.h"

#include "test_support.h"
#include "array_sort.h"
#include "array_norm.h"

/*
 *	Purpose:
 *		Test correctness of implementation of array_sort.
 *
 *	Comments:
 *		The ordering of the (Ind)ices is partially arbitrary if certain columns of A are identical.
 *
 *	Notation:
 *
 *	References:
 */

void test_unit_array_sort(void)
{
	unsigned int pass;

	/*
	 *	array_sort_ui:
	 *
	 *		Input:
	 *
	 *			A = [0 1 1 2 3 0 0 1 2 1
	 *			     1 3 3 5 1 0 1 0 4 3
	 *			     3 1 4 3 0 1 2 1 2 5]
	 *
	 *			Ind = [0 1 2 3 4 5 6 7 8 9]
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
	 */

	unsigned int i;
	unsigned int A_ui[3*10],
	             A_ui0[3*10] = { 0, 1, 1, 2, 3, 0, 0, 1, 2, 1, 1, 3, 3, 5, 1, 0, 1, 0, 4, 3, 3, 1, 4, 3, 0, 1, 2, 1, 2, 5 },
	             A_ui1[3*10] = { 0, 0, 0, 1, 1, 1, 1, 2, 2, 3 },                               // Only 1st row is sorted
	             A_ui2[3*10] = { 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 0, 1, 1, 0, 3, 3, 3, 4, 5, 1 }, // Only 2nd row is sorted
	             A_ui3[3*10] = { 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 0, 1, 1, 0, 3, 3, 3, 4, 5, 1, 1, 2, 3, 1, 1, 4, 5, 2, 3, 0 };
	unsigned int Ind_ui[10],
	             Ind_ui0[10] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 },
	             Ind_ui3[10] = { 5, 6, 0, 7, 1, 2, 9, 8, 3, 4 };


	// row(s) sorted: 1
	for (i = 0; i < 30; i++) A_ui[i]   = A_ui0[i];
	for (i = 0; i < 10; i++) Ind_ui[i] = Ind_ui0[i];

	array_sort_ui(3,10,A_ui,Ind_ui,'R','N');

	pass = 0;
	if (array_norm_diff_ui(10,A_ui,A_ui1,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("array_sort_ui (1 row(s) sorted):                 ");
	test_print(pass);

	// row(s) sorted: 2
	for (i = 0; i < 30; i++) A_ui[i]   = A_ui0[i];
	for (i = 0; i < 10; i++) Ind_ui[i] = Ind_ui0[i];

	array_sort_ui(3,10,A_ui,Ind_ui,'R','N');

	pass = 0;
	if (array_norm_diff_ui(20,A_ui,A_ui2,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("              (2 row(s) sorted):                 ");
	test_print(pass);

	// row(s) sorted: 3
	for (i = 0; i < 30; i++) A_ui[i]   = A_ui0[i];
	for (i = 0; i < 10; i++) Ind_ui[i] = Ind_ui0[i];

	array_sort_ui(3,10,A_ui,Ind_ui,'R','N');

	pass = 0;
	if (array_norm_diff_ui(30,A_ui,A_ui3,"Inf") < EPS &&
		array_norm_diff_ui(10,Ind_ui,Ind_ui3,"Inf") < EPS)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("              (3 row(s) sorted):                 ");
	test_print(pass);


	// Transposed
	// row(s) sorted: 1
	for (i = 0; i < 30; i++) A_ui[i]   = A_ui0[i];
	for (i = 0; i < 10; i++) Ind_ui[i] = Ind_ui0[i];

	// Transpose, sort, transpose
	mkl_simatcopy('R','T',3,10,1.,(float *) A_ui,10,3);
	array_sort_ui(10,3,A_ui,Ind_ui,'R','T');
	mkl_simatcopy('R','T',10,3,1.,(float *) A_ui,3,10);

	pass = 0;
	if (array_norm_diff_ui(10,A_ui,A_ui1,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("array_sort_ui - Transposed (1 row(s) sorted):    ");
	test_print(pass);

	// row(s) sorted: 2
	for (i = 0; i < 30; i++) A_ui[i]   = A_ui0[i];
	for (i = 0; i < 10; i++) Ind_ui[i] = Ind_ui0[i];

	// Transpose, sort, transpose
	mkl_simatcopy('R','T',3,10,1.,(float *) A_ui,10,3);
	array_sort_ui(10,3,A_ui,Ind_ui,'R','T');
	mkl_simatcopy('R','T',10,3,1.,(float *) A_ui,3,10);

	pass = 0;
	if (array_norm_diff_ui(20,A_ui,A_ui2,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("                           (2 row(s) sorted):    ");
	test_print(pass);

	// row(s) sorted: 3
	for (i = 0; i < 30; i++) A_ui[i]   = A_ui0[i];
	for (i = 0; i < 10; i++) Ind_ui[i] = Ind_ui0[i];

	// Transpose, sort, transpose
	mkl_simatcopy('R','T',3,10,1.,(float *) A_ui,10,3);
	array_sort_ui(10,3,A_ui,Ind_ui,'R','T');
	mkl_simatcopy('R','T',10,3,1.,(float *) A_ui,3,10);

	pass = 0;
	if (array_norm_diff_ui(30,A_ui,A_ui3,"Inf") < EPS &&
		array_norm_diff_ui(10,Ind_ui,Ind_ui3,"Inf") < EPS)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("                           (3 row(s) sorted):    ");
	test_print(pass);

	unsigned int NRows = 30, NCols = 4;
	unsigned int B_ui[4*30],
	             B_ui0[4*30] = { 0, 2, 8, 9, 1, 3, 8, 9, 0, 1, 8, 9, 2, 3, 8, 9, 0, 1, 2, 3,
	                             6, 7, 8, 9, 4, 5, 8, 9, 5, 7, 8, 9, 4, 6, 8, 9, 4, 5, 6, 7,
	                             4, 5, 8, 9, 0, 1, 8, 9, 1, 5, 8, 9, 0, 4, 8, 9, 0, 1, 4, 5,
	                             2, 6, 8, 9, 3, 7, 8, 9, 2, 3, 8, 9, 6, 7, 8, 9, 2, 3, 6, 7,
	                             0, 4, 8, 9, 2, 6, 8, 9, 0, 2, 8, 9, 4, 6, 8, 9, 0, 2, 4, 6,
	                             5, 7, 8, 9, 1, 3, 8, 9, 3, 7, 8, 9, 1, 5, 8, 9, 1, 3, 5, 7},
	             B_ui1[4*30] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3,
	                             3, 4, 4, 4, 4, 4, 5, 5, 6, 6},
	             B_ui2[4*30] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3,
	                             3, 4, 4, 4, 4, 4, 5, 5, 6, 6, 1, 1, 1, 1, 2, 2, 2, 4, 4, 3,
	                             3, 3, 5, 5, 3, 3, 3, 6, 6, 7, 7, 5, 5, 5, 6, 6, 7, 7, 7, 7},
	             B_ui3[4*30] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3,
	                             3, 4, 4, 4, 4, 4, 5, 5, 6, 6, 1, 1, 1, 1, 2, 2, 2, 4, 4, 3,
	                             3, 3, 5, 5, 3, 3, 3, 6, 6, 7, 7, 5, 5, 5, 6, 6, 7, 7, 7, 7,
	                             2, 4, 8, 8, 4, 8, 8, 8, 8, 5, 8, 8, 8, 8, 6, 8, 8, 8, 8, 8,
	                             8, 6, 8, 8, 8, 8, 8, 8, 8, 8},
	             B_ui4[4*30] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3,
	                             3, 4, 4, 4, 4, 4, 5, 5, 6, 6, 1, 1, 1, 1, 2, 2, 2, 4, 4, 3,
	                             3, 3, 5, 5, 3, 3, 3, 6, 6, 7, 7, 5, 5, 5, 6, 6, 7, 7, 7, 7,
	                             2, 4, 8, 8, 4, 8, 8, 8, 8, 5, 8, 8, 8, 8, 6, 8, 8, 8, 8, 8,
	                             8, 6, 8, 8, 8, 8, 8, 8, 8, 8, 3, 5, 9, 9, 6, 9, 9, 9, 9, 7,
	                             9, 9, 9, 9, 7, 9, 9, 9, 9, 9, 9, 7, 9, 9, 9, 9, 9, 9, 9, 9};
	unsigned int IndB_ui[30],
	             IndB_ui0[30] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
	                              20, 21, 22, 23, 24, 25, 26, 27, 28, 29 },
	             IndB_ui4[30] = { 4, 14, 11, 2, 24, 0, 22, 13, 20, 29, 26, 1, 12, 28, 19, 17, 3, 15, 21, 16,
	                              27, 9, 10, 6, 23, 8, 25, 7, 5, 18};

	// Test 2
	// Store B in column-major format
	mkl_simatcopy('R','T',NRows,NCols,1.,(float *) B_ui0,NCols,NRows);

	// row(s) sorted: 1
	for (i = 0; i < NRows*NCols; i++) B_ui[i]    = B_ui0[i];
	for (i = 0; i < NRows; i++)       IndB_ui[i] = IndB_ui0[i];

	array_sort_ui(NCols,NRows,B_ui,IndB_ui,'R','N');

	pass = 0;
	if (array_norm_diff_ui(NRows*1,B_ui,B_ui1,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("array_sort_ui - test 2 (1 row(s) sorted):        ");
	test_print(pass);

	// row(s) sorted: 2
	for (i = 0; i < NRows*NCols; i++) B_ui[i]    = B_ui0[i];
	for (i = 0; i < NRows; i++)       IndB_ui[i] = IndB_ui0[i];

	array_sort_ui(NCols,NRows,B_ui,IndB_ui,'R','N');

	pass = 0;
	if (array_norm_diff_ui(NRows*2,B_ui,B_ui2,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("                       (2 row(s) sorted):        ");
	test_print(pass);

	// row(s) sorted: 3
	for (i = 0; i < NRows*NCols; i++) B_ui[i]    = B_ui0[i];
	for (i = 0; i < NRows; i++)       IndB_ui[i] = IndB_ui0[i];

	array_sort_ui(NCols,NRows,B_ui,IndB_ui,'R','N');

	pass = 0;
	if (array_norm_diff_ui(NRows*3,B_ui,B_ui3,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("                       (3 row(s) sorted):        ");
	test_print(pass);

	// row(s) sorted: 4
	for (i = 0; i < NRows*NCols; i++) B_ui[i]    = B_ui0[i];
	for (i = 0; i < NRows; i++)       IndB_ui[i] = IndB_ui0[i];

	array_sort_ui(NCols,NRows,B_ui,IndB_ui,'R','N');

	pass = 0;
	if (array_norm_diff_ui(NRows*4,B_ui,B_ui4,"Inf")       < EPS &&
	    array_norm_diff_ui(NRows*1,IndB_ui,IndB_ui4,"Inf") < EPS)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("                       (4 row(s) sorted):        ");
	test_print(pass);

	/*
	 *	array_sort_d:
	 *
	 *		Same as first example above cast to double.
	 */

	double A_d[3*10], A_d0[3*10], A_d1[3*10], A_d2[3*10], A_d3[3*10];

	for (i = 0; i < 30; i++) {
		A_d0[i] = (double) A_ui0[i];
		A_d1[i] = (double) A_ui1[i];
		A_d2[i] = (double) A_ui2[i];
		A_d3[i] = (double) A_ui3[i];
	}

	// row(s) sorted: 0
	for (i = 0; i < 30; i++) A_d[i]   = A_d0[i];
	for (i = 0; i < 10; i++) Ind_ui[i] = Ind_ui0[i];

	pass = 0;
	if (array_norm_diff_d(30,A_d,A_d0,"Inf") < EPS &&
		array_norm_diff_ui(10,Ind_ui,Ind_ui0,"Inf") < EPS)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("array_sort_d (0 row(s) sorted):                  ");
	test_print(pass);


	// row(s) sorted: 1
	for (i = 0; i < 30; i++) A_d[i]   = A_d0[i];
	for (i = 0; i < 10; i++) Ind_ui[i] = Ind_ui0[i];

	array_sort_d(3,10,A_d,Ind_ui,'R','N');

	pass = 0;
	if (array_norm_diff_d(10,A_d,A_d1,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("             (1 row(s) sorted):                  ");
	test_print(pass);

	// row(s) sorted: 2
	for (i = 0; i < 30; i++) A_d[i]   = A_d0[i];
	for (i = 0; i < 10; i++) Ind_ui[i] = Ind_ui0[i];

	array_sort_d(3,10,A_d,Ind_ui,'R','N');

	pass = 0;
	if (array_norm_diff_d(20,A_d,A_d2,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("             (2 row(s) sorted):                  ");
	test_print(pass);

	// row(s) sorted: 3
	for (i = 0; i < 30; i++) A_d[i]   = A_d0[i];
	for (i = 0; i < 10; i++) Ind_ui[i] = Ind_ui0[i];

	array_sort_d(3,10,A_d,Ind_ui,'R','N');

	pass = 0;
	if (array_norm_diff_d(30,A_d,A_d3,"Inf") < EPS &&
		array_norm_diff_ui(10,Ind_ui,Ind_ui3,"Inf") < EPS)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("             (3 row(s) sorted):                  ");
	test_print(pass);


	// Transposed
	// row(s) sorted: 0
	for (i = 0; i < 30; i++) A_d[i]   = A_d0[i];
	for (i = 0; i < 10; i++) Ind_ui[i] = Ind_ui0[i];

	// Transpose, sort, transpose
	mkl_dimatcopy('R','T',3,10,1.,A_d,10,3);
	mkl_dimatcopy('R','T',10,3,1.,A_d,3,10);

	pass = 0;
	if (array_norm_diff_d(30,A_d,A_d0,"Inf") < EPS &&
		array_norm_diff_ui(10,Ind_ui,Ind_ui0,"Inf") < EPS)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("array_sort_d - Transposed (0 row(s) sorted):     ");
	test_print(pass);

	// row(s) sorted: 1
	for (i = 0; i < 30; i++) A_d[i]   = A_d0[i];
	for (i = 0; i < 10; i++) Ind_ui[i] = Ind_ui0[i];

	// Transpose, sort, transpose
	mkl_dimatcopy('R','T',3,10,1.,A_d,10,3);
	array_sort_d(10,3,A_d,Ind_ui,'R','T');
	mkl_dimatcopy('R','T',10,3,1.,A_d,3,10);

	pass = 0;
	if (array_norm_diff_d(10,A_d,A_d1,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("                          (1 row(s) sorted):     ");
	test_print(pass);

	// row(s) sorted: 2
	for (i = 0; i < 30; i++) A_d[i]   = A_d0[i];
	for (i = 0; i < 10; i++) Ind_ui[i] = Ind_ui0[i];

	// Transpose, sort, transpose
	mkl_dimatcopy('R','T',3,10,1.,A_d,10,3);
	array_sort_d(10,3,A_d,Ind_ui,'R','T');
	mkl_dimatcopy('R','T',10,3,1.,A_d,3,10);

	pass = 0;
	if (array_norm_diff_d(20,A_d,A_d2,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("                          (2 row(s) sorted):     ");
	test_print(pass);

	// row(s) sorted: 3
	for (i = 0; i < 30; i++) A_d[i]   = A_d0[i];
	for (i = 0; i < 10; i++) Ind_ui[i] = Ind_ui0[i];

	// Transpose, sort, transpose
	mkl_dimatcopy('R','T',3,10,1.,A_d,10,3);
	array_sort_d(10,3,A_d,Ind_ui,'R','T');
	mkl_dimatcopy('R','T',10,3,1.,A_d,3,10);

	pass = 0;
	if (array_norm_diff_d(30,A_d,A_d3,"Inf") < EPS &&
		array_norm_diff_ui(10,Ind_ui,Ind_ui3,"Inf") < EPS)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("                          (3 row(s) sorted):     ");
	test_print(pass);
}
