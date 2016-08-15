// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "test_unit_jacobian_fluxes_inviscid.h"

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>

#include "Parameters.h"
#include "Test.h"

#include "test_code_fluxes_inviscid.h"
#include "test_support.h"
#include "array_norm.h"
#include "fluxes_inviscid_c.h"
#include "jacobian_fluxes_inviscid.h"

/*
 *	Purpose:
 *		Test correctness of implementation of functions relating to jacobians of inviscid fluxes.
 *
 *	Comments:
 *		Correctness is assessed based on comparison with derivatives computed using the complex step method.
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
		for (dim = 0; dim < d; dim++) {
			IndF    = (eq*d+dim)*NnTotal;
			InddFdW = ((eq*Neq+var)*d+dim)*NnTotal;
			for (n = 0; n < NnTotal; n++)
				dFdW_cs[InddFdW+n] = cimag(F[IndF+n])/h;
		}}
	}
	free(F);
}

static unsigned int compare_jacobian_flux_inviscid(const unsigned int Nn, const unsigned int Nel, const unsigned int d,
                                                   const unsigned int Neq, double *W)
{
	unsigned int pass = 0;

	unsigned int NnTotal, Nvar;
	double       *dFdW, *dFdW_cs;

	NnTotal = Nn*Nel;
	Nvar    = Neq;

	dFdW    = malloc(NnTotal*d*Nvar*Neq * sizeof *dFdW);    // free
	dFdW_cs = malloc(NnTotal*d*Nvar*Neq * sizeof *dFdW_cs); // free

	jacobian_flux_inviscid(Nn,Nel,W,dFdW,d,Neq);
	compute_dFdW_cs(Neq,Nn,Nel,d,W,dFdW_cs);

	if (array_norm_diff_d(NnTotal*d*Nvar*Neq,dFdW,dFdW_cs,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	free(dFdW);
	free(dFdW_cs);

	return pass;
}

void test_unit_jacobian_fluxes_inviscid(void)
{
	unsigned int pass;

	/*
	 *	Input:
	 *		Jacobians computed using the linearized code and the complex step method.
	 *
	 *	Expected Output:
	 *		No difference between the results.
	 *
	 */

	unsigned int Nn, Nel, d, Neq;
	double       *W;

	for (d = 1; d <= 3; d++) {
		Neq = d+2;

		W    = initialize_W(&Nn,&Nel,d); // free
		pass = compare_jacobian_flux_inviscid(Nn,Nel,d,Neq,W);

		if (d == 1) {
			//     0         10        20        30        40        50
			printf("jacobian_flux_inviscid (d = %d):                  ",d);
		} else {
			printf("                       (d = %d):                  ",d);
		}
		test_print(pass);

		free(W);
	}
}
