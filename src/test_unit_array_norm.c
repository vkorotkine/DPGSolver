// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "test_unit_array_norm.h"

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

void test_unit_array_norm(void)
{
	unsigned int pass;

	/*
	 *	array_norm_ui:
	 *
	 *		Input:
	 *
	 *			A = [1 2 3]
	 *
	 *		Expected output:
	 *
	 *			Inf : 3
	 *			L1  : 6
	 */

	unsigned int A_ui[3] = { 1, 2, 3 };

	pass = 0;
	if ((array_norm_ui(3,A_ui,"Inf") - 3) < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("array_norm_ui (Inf):                             ");
	test_print(pass);

	pass = 0;
	if ((array_norm_ui(3,A_ui,"L1") - 6) < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("              (L1) :                             ");
	test_print(pass);

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

	double A_d[3] = { 1.0, 2.0, 3.0 };

	pass = 0;
	if ((array_norm_d(3,A_d,"Inf") - 3.0) < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("array_norm_d (Inf):                              ");
	test_print(pass);

	pass = 0;
	if ((array_norm_d(3,A_d,"L1") - 6.0) < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("             (L1) :                              ");
	test_print(pass);

	pass = 0;
	if ((array_norm_d(3,A_d,"L2") - 3.741657386773941) < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("             (L2) :                              ");
	test_print(pass);

	/*
	 *	array_norm_diff_ui:
	 *
	 *		Input:
	 *
	 *			A = [1 2 3]
	 *			B = [4 5 6]
	 *
	 *			Inf : 3
	 *			L1  : 9
	 *
	 */

	unsigned int A_diff_ui[3] = { 1, 2, 3 };
	unsigned int B_diff_ui[3] = { 4, 5, 6 };

	pass = 0;
	if ((array_norm_diff_ui(3,A_diff_ui,B_diff_ui,"Inf") - 3) < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("array_norm_diff_ui (Inf):                        ");
	test_print(pass);

	pass = 0;
	if ((array_norm_diff_ui(3,A_diff_ui,B_diff_ui,"L1") - 9) < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("                   (L1) :                        ");
	test_print(pass);

	/*
	 *	array_norm_diff_d:
	 *
	 *		Input:
	 *
	 *			A = [1.0 2.0 3.0]
	 *			B = [4.0 5.0 6.0]
	 *
	 *			Inf : 3.0
	 *			L1  : 9.0
	 *			L2  : 5.196152422706632
	 *
	 */

	double A_diff_d[3] = { 1.0, 2.0, 3.0 };
	double B_diff_d[3] = { 4.0, 5.0, 6.0 };

	pass = 0;
	if ((array_norm_diff_d(3,A_diff_d,B_diff_d,"Inf") - 3.0) < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("array_norm_diff_d (Inf):                         ");
	test_print(pass);

	pass = 0;
	if ((array_norm_diff_d(3,A_diff_d,B_diff_d,"L1") - 9.0) < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("                  (L1) :                         ");
	test_print(pass);

	pass = 0;
	if ((array_norm_diff_d(3,A_diff_d,B_diff_d,"L2") - 5.196152422706632) < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("                  (L2) :                         ");
	test_print(pass);
}
