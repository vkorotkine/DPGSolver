#include <stdlib.h>
#include <stdio.h>

#include "test.h"
#include "functions.h"
#include "parameters.h"

/*
 *	Purpose:
 *		Test correctness of implementation of basis_TRI.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

static double *basis_TRI(const unsigned int P, const double *rst, const unsigned int Nn)
{
	/* ChiRef_rst = @(r,s) [
	 */

	unsigned int i, j, jMax, d = 2;
	double *a, *b, *c;

	// Convert from rst to abc coordinates
	a = calloc(Nn , sizeof *a); // free
	b = calloc(Nn , sizeof *b); // free
	c = calloc(Nn , sizeof *c); // free

	rst_to_abc(Nn,d,*rst,*a,*b,*c);

	// This should be moved to the actual function, use the actual values here (i.e. (1-b)^0, (1-b)^1, etc)
	for (k = 0; k < Nn; k++) {
		for (i = 0; i <= P; i++) {
		for (j = 0, jMax = P-i; j <= jMax; j++) {


		}}
	}

	free(a);
	free(b);
	free(c);
}

void test_imp_basis_TRI(void)
{
	unsigned int pass;

	/*
	 *	Input:
	 *
	 *		P, rst, Nn.
	 *
	 *	Expected output:
	 *
	 *		P = 2:
	 *			ChiRef_rst = @(r,s) [
	 *
	 *		P = 3:
	 *			ChiRef_rst = @(r,s) [
	 */

	unsigned int d, Nn, Ns, Nbf, P, Prst;
	unsigned int *symms;
	double *rst, *w;

	d = 2;

	double *ChiRef2_code, *ChiRef2_test;

	P = 2;
	Prst = 4;

	cubature_TRI(&rst,&w,&symms,&Nn,&Ns,0,Prst,d,"G"); // free

//	ChiRef13_code = basis_TRI(P,rst,Nn,&Nbf,d); // free
	ChiRef13_test = basis_TRI(P,rst,Nn);         // free



	pass = 0;
/*	if (array_norm_diff_d(Nn*Nbf,ChiRef2_code,ChiRef2_test,"Inf") < EPS*10)
		pass = 1, TestDB.Npass++;
*/

	//     0         10        20        30        40        50
	printf("basis_TRI (P2):                                  ");
	test_print(pass);

	free(rst);
	free(symms);
//	free(ChiRef2_code);
//	free(ChiRef2_test);

}
