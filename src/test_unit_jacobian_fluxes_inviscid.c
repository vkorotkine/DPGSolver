// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "test_unit_jacobian_fluxes_inviscid.h"

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>

#include "Parameters.h"
#include "Test.h"

#include "test_support.h"
#include "array_norm.h"
#include "fluxes_inviscid_c.h"
#include "jacobian_fluxes_inviscid.h"

/*
 *	Purpose:
 *		Test correctness of implementation of functions relating to jacobians of inviscid fluxes.
 *
 *	Comments:
 *		Correctness is assessed based on comparison of derivatives compute using the complex step method.
 *
 *	Notation:
 *
 *	References:
 *		Squire(1998)-Using_Complex_Variables_to_Estimate_Derivatives_of_Real_Functions
 */

static void compute_dFdW_cs(const unsigned int Neq, const unsigned int Nn, const unsigned int Nel, const unsigned int d,
                            double *W, double *dFdW_cs)
{
	unsigned int   n, var, eq, dim, Nvar, NnTotal, IndW, IndF, InddFdW;
	double         h;
	double complex Wp[Nn*Nel*Neq], *F;

	h = EPS*EPS;

	Nvar    = Neq;
	NnTotal = Nn*Nel;

	F = malloc(NnTotal*d*Neq * sizeof *F); // free

	for (var = 0; var < Nvar; var++) {
		for (eq = 0; eq < Neq; eq++) {
			IndW = NnTotal*eq;
			for (n = 0; n < NnTotal; n++) {
				Wp[IndW+n] = W[IndW+n];
				if (eq == var)
					Wp[IndW+n] += h*I;
			}
		}
		flux_inviscid_c(Nn,Nel,Wp,F,d,Neq);

		for (eq = 0; eq < Neq; eq++) {
			IndF = NnTotal*eq;
			for (dim = 0; dim < d; dim++) {
				IndF    = (eq*d+dim)*NnTotal;
				InddFdW = ((eq*Neq+var)*d+dim)*NnTotal;
				for (n = 0; n < NnTotal; n++) {
					dFdW_cs[InddFdW+n] = cimag(F[IndF+n])/h;
				}
			}
		}
	}
	free(F);
}

void test_unit_jacobian_fluxes_inviscid(void)
{
	unsigned int pass;

	/*
	 *	Input:
	 *
	 *	Expected Output:
	 *
	 */

	unsigned int   i,  Nn, Nel, iMax, d, Neq, Nvar, NnTotal;
	double         *W, *dFdW, *dFdW_cs;

	Nn  = 3;
	Nel = 2;

	NnTotal = Nn*Nel;


	// d = 1
	d    = 1;
	Neq  = d+2;
	Nvar = Neq;

	double W1[18] = {  1.01,  1.02,  1.03,  1.11,  1.12,  1.13,
	                   2.01, -2.02,  2.03, -2.11,  2.12, -2.13,
	                   5.01,  5.02,  5.03,  5.11,  5.12,  5.13};

	dFdW    = malloc(NnTotal*d*Nvar*Neq * sizeof *dFdW);    // free
	dFdW_cs = malloc(NnTotal*d*Nvar*Neq * sizeof *dFdW_cs); // free

	W = malloc(NnTotal*Nvar * sizeof *W); // free
	for (i = 0, iMax = NnTotal*Nvar; i < iMax; i++)
		W[i] = W1[i];

	jacobian_flux_inviscid(Nn,Nel,W,dFdW,d,Neq);
	compute_dFdW_cs(Neq,Nn,Nel,d,W,dFdW_cs);

	pass = 0;
	if (array_norm_diff_d(NnTotal*d*Nvar*Neq,dFdW,dFdW_cs,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("jacobian_flux_inviscid (d = 1):                  ");
	test_print(pass);

	free(W);
	free(dFdW);
	free(dFdW_cs);


	// d = 2
	d    = 2;
	Neq  = d+2;
	Nvar = Neq;

	double W2[24] = {  1.01,  1.02,  1.03,  1.11,  1.12,  1.13,
	                   2.01, -2.02,  2.03, -2.11,  2.12, -2.13,
	                   3.04, -3.05,  3.06, -3.14,  3.15, -3.16,
	                   5.01,  5.02,  5.03,  5.11,  5.12,  5.13};

	dFdW    = malloc(NnTotal*d*Nvar*Neq * sizeof *dFdW);    // free
	dFdW_cs = malloc(NnTotal*d*Nvar*Neq * sizeof *dFdW_cs); // free

	W = malloc(NnTotal*Nvar * sizeof *W); // free
	for (i = 0, iMax = NnTotal*Nvar; i < iMax; i++)
		W[i] = W2[i];

	jacobian_flux_inviscid(Nn,Nel,W,dFdW,d,Neq);
	compute_dFdW_cs(Neq,Nn,Nel,d,W,dFdW_cs);

	pass = 0;
	if (array_norm_diff_d(NnTotal*d*Nvar*Neq,dFdW,dFdW_cs,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("                       (d = 2):                  ");
	test_print(pass);

	free(W);
	free(dFdW);
	free(dFdW_cs);


	// d = 3
	d    = 3;
	Neq  = d+2;
	Nvar = Neq;

	double W3[30] = {  1.01,  1.02,  1.03,  1.11,  1.12,  1.13,
	                   2.01, -2.02,  2.03, -2.11,  2.12, -2.13,
	                   3.04, -3.05,  3.06, -3.14,  3.15, -3.16,
	                   4.07, -4.08,  4.09, -4.17,  4.18, -4.19,
	                   5.01,  5.02,  5.03,  5.11,  5.12,  5.13};

	dFdW    = malloc(NnTotal*d*Nvar*Neq * sizeof *dFdW);    // free
	dFdW_cs = malloc(NnTotal*d*Nvar*Neq * sizeof *dFdW_cs); // free

	W = malloc(NnTotal*Nvar * sizeof *W); // free
	for (i = 0, iMax = NnTotal*Nvar; i < iMax; i++)
		W[i] = W3[i];

	jacobian_flux_inviscid(Nn,Nel,W,dFdW,d,Neq);
	compute_dFdW_cs(Neq,Nn,Nel,d,W,dFdW_cs);

	pass = 0;
	if (array_norm_diff_d(NnTotal*d*Nvar*Neq,dFdW,dFdW_cs,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("                       (d = 3):                  ");
	test_print(pass);

	free(W);
	free(dFdW);
	free(dFdW_cs);
}
