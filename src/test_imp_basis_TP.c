#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "test.h"
#include "functions.h"
#include "parameters.h"

/*
 *	Purpose:
 *		Test correctness of implementation of basis_TP.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 *		GL  : http://mathworld.wolfram.com/Legendre-GaussQuadrature.html
 *		GLL : http://mathworld.wolfram.com/LobattoQuadrature.html
 */

void test_imp_basis_TP(void)
{
	unsigned int pass;

	/*
	 *	basis_TP (dE = 1):
	 *
	 *		Input:
	 *
	 *			P, rst, Nn, dE.
	 *
	 *		Expected output:
	 *
	 *			P = 3, rst14_GL:
	 *				ChiRef_rst = @(r) [sqrt(1/2)*1 sqrt(3/2)*r sqrt(5/2)*1/2*(3*r^2-1) sqrt(7/2)*1/2*(5*r^3-3*r)]
	 */

	unsigned int dE, Nn, P, ToReturn[4];
	unsigned int *Con;
	double *rst, *W;

	dE = 1;

	double *ChiRef_rst14_GL;

	P = 4;

	ToReturn[0] = 1; ToReturn[1] = 0; ToReturn[2] = 0; ToReturn[3] = 1;
	cubature_TP(&rst,&W,&Con,&Nn,ToReturn,P,dE,"GL"); // tbd

	ChiRef_rst14_GL = basis_TP(3,rst,Nn,dE);

	array_print_d(Nn,4,ChiRef_rst14_GL,'R');



/*
	pass = 0;
	if (array_norm_diff_d(N13,rst13_ES,rst,"Inf") < EPS &&
		array_norm_diff_ui(2*Nn,Con13,Con,"Inf")  < EPS &&
		Nn == P)
			pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("cubature_TP (d1, P3, ES):                        ");
	test_print(pass);
*/
}
