// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>

#include "test.h"
#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Test correctness of implementation of functions relating to jacobians of inviscid fluxes.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void test_imp_jacobian_fluxes_inviscid(void)
{
	unsigned int pass;

	/*
	 *	Input:
	 *
	 *	Expected Output:
	 *
	 */

	unsigned int Nn, Nel, d, Neq;
	double       *dFdW;

	Nn  = 3;
	Nel = 2;
/*
	// d = 1
	d   = 1;
	Neq = d+2;

	double W1[18] = {  1.01,  1.02,  1.03,  1.11,  1.12,  1.13,
	                   2.01, -2.02,  2.03, -2.11,  2.12, -2.13,
	                   5.01,  5.02,  5.03,  5.11,  5.12,  5.13};

	//     0         10        20        30        40        50
	printf("jacobian_flux_inviscid (d = 1):                  ");
	test_print(pass);

	free(dFdW);

	// d = 2
	d   = 2;
	Neq = d+2;

	double W2[24] = {  1.01,  1.02,  1.03,  1.11,  1.12,  1.13,
	                   2.01, -2.02,  2.03, -2.11,  2.12, -2.13,
	                   3.04, -3.05,  3.06, -3.14,  3.15, -3.16,
	                   5.01,  5.02,  5.03,  5.11,  5.12,  5.13};

	//     0         10        20        30        40        50
	printf("                       (d = 2):                  ");
	test_print(pass);

	free(dFdW);
*/
	// d = 3
	d   = 3;
	Neq = d+2;

	double W3[30] = {  1.01,  1.02,  1.03,  1.11,  1.12,  1.13,
	                   2.01, -2.02,  2.03, -2.11,  2.12, -2.13,
	                   3.04, -3.05,  3.06, -3.14,  3.15, -3.16,
	                   4.07, -4.08,  4.09, -4.17,  4.18, -4.19,
	                   5.01,  5.02,  5.03,  5.11,  5.12,  5.13};

	dFdW = malloc(Nn*Neq*Neq*d*Nel * sizeof *dFdW); // free
	jacobian_flux_inviscid(Nn,Nel,W3,dFdW,d,Neq);

unsigned int var, eq;
for (eq = 0; eq < Neq; eq++) {
//printf("\n\n\nEq: %d\n",eq);
for (var = 0; var < Neq; var++) {
//array_print_d(Nn*Nel,d,&dFdW[(Nn*Nel*d)*var+(Nn*Nel*d*Neq)*eq],'C');
}}

//	pass = 0;
//	if (array_norm_diff_d(Nn*Neq*d*Nel,F,F3,"Inf") < 10*EPS)
//		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("              (d = 3):                           ");
	test_print(pass);

	free(dFdW);
}
