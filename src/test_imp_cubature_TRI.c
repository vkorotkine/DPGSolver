#include <stdlib.h>
#include <stdio.h>
//#include <math.h>

#include "test.h"
#include "functions.h"
//#include "parameters.h"

/*
 *	Purpose:
 *		Test correctness of implementation of cubature_TRI.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 *		ToBeModified (Likely included in the functions where nodes/weights are stored and would be redundant here.
 */

void test_imp_cubature_TRI(void)
{
	unsigned int pass;

	/*
	 *	Input:
	 *
	 *		P, NodeType.
	 *
	 *	Expected output:
	 *
	 *		P = 3, AO:
	 *			rst = [ See below ]
	 *			w   = N/A
	 *			Nn  = ToBeModified
	 *
	 *		P = 3, Cools:
	 *			rst = [ See below ]
	 *			w   = [ See below ]
	 *			Nn  = ToBeModified
	 */

	unsigned int i, j, k, row;
	unsigned int Nn, P;
	double *rst, *w;

	// P = 3
	P = 3;

	unsigned int Nn3_AO = 1, Nn3_Cools = 1;
	double rst3_AO[Nn3_AO], rst3_Cools[Nn3_Cools], w3_Cools[Nn3_Cools];

	// AO
//	cubature_TRI(&rst,&w,&Nn,1,P,"AO"); // tbd

/*
	pass = 0;
	if (array_norm_diff_d(N13,rst13_GL,rst,"Inf") < EPS &&
		array_norm_diff_d(N13,w13_GL,w,"Inf")     < EPS &&
		Nn == N13)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("cubature_TP (d1, P3, GL):                        ");
	test_print(pass);
*/

	free(rst), free(w);

	// P = 4?

}
