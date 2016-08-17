// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "test_unit_jacobian_fluxes_inviscid.h"

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>

#include "Parameters.h"
#include "Macros.h"
#include "Test.h"

#include "test_code_fluxes_inviscid.h"
#include "test_support.h"
#include "array_norm.h"
#include "fluxes_inviscid_c.h"
#include "jacobian_fluxes_inviscid.h"

#include "array_print.h"

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

static void compute_dnFdW_cs(const unsigned int Nn, const unsigned int Nel, const unsigned int d,
                             const unsigned int Neq, double *WL, double *WR, double *dnFdWL_cs, double *dnFdWR_cs,
                             double *nL, const char *nFType)
{
	unsigned int   i, iMax, n, var, eq, Nvar, NnTotal, IndW, IndnF, InddnFdW;
	double         h;
	double complex WLp[Nn*Nel*Neq], WRp[Nn*Nel*Neq], *nF;

	h = EPS*EPS;

	Nvar    = Neq;
	NnTotal = Nn*Nel;

	nF = malloc(NnTotal*Neq * sizeof *nF); // free

	for (i = 0, iMax = NnTotal*Nvar; i < iMax; i++) {
		WLp[i] = WL[i];
		WRp[i] = WR[i];
	}

	for (var = 0; var < Nvar; var++) {
		// Left
		IndW = NnTotal*var;
		for (n = 0; n < NnTotal; n++)
			WLp[IndW+n] += h*I;

		if (strstr(nFType,"LF"))
			flux_LF_c(Nn,Nel,WLp,WRp,nF,nL,d,Neq);
		else if (strstr(nFType,"Roe"))
			printf("Add support Roe.\n"), EXIT_MSG;
		else
			printf("Error: Unsupported nFType.\n"), EXIT_MSG;

		for (eq = 0; eq < Neq; eq++) {
			IndnF    = eq*NnTotal;
			InddnFdW = (eq*Neq+var)*NnTotal;
			for (n = 0; n < NnTotal; n++)
				dnFdWL_cs[InddnFdW+n] = cimag(nF[IndnF+n])/h;
		}

		// Right
		for (n = 0; n < NnTotal; n++) {
			WLp[IndW+n] -= h*I;
			WRp[IndW+n] += h*I;
		}

		if (strstr(nFType,"LF"))
			flux_LF_c(Nn,Nel,WLp,WRp,nF,nL,d,Neq);
		else if (strstr(nFType,"Roe"))
			printf("Add support Roe.\n"), EXIT_MSG;
		else
			printf("Error: Unsupported nFType.\n"), EXIT_MSG;

		for (eq = 0; eq < Neq; eq++) {
			IndnF    = eq*NnTotal;
			InddnFdW = (eq*Neq+var)*NnTotal;
			for (n = 0; n < NnTotal; n++)
				dnFdWR_cs[InddnFdW+n] = cimag(nF[IndnF+n])/h;
		}

		for (n = 0; n < NnTotal; n++)
			WRp[IndW+n] -= h*I;
	}
	free(nF);
}

static unsigned int compare_jacobian_flux_Num(const unsigned int Nn, const unsigned int Nel, const unsigned int d,
                                              const unsigned int Neq, double *W, double *nL, const char *nFType)
{
	unsigned int pass = 0;

	unsigned int i, n, var, Nvar, IndWLR, CheckedAllLF;
	double       *W_ptr, *WL, *WR, *dnFdWL, *dnFdWR, *dnFdWL_cs, *dnFdWR_cs;

	if (Nel != 2)
		printf("Error: Unsupported Nel.\n"), EXIT_MSG;

	Nvar = Neq;

	WL        = malloc(Nn*Nvar     * sizeof *WL);        // free
	WR        = malloc(Nn*Nvar     * sizeof *WR);        // free
	dnFdWL    = malloc(Nn*Nvar*Neq * sizeof *dnFdWL);    // free
	dnFdWR    = malloc(Nn*Nvar*Neq * sizeof *dnFdWR);    // free
	dnFdWL_cs = malloc(Nn*Nvar*Neq * sizeof *dnFdWL_cs); // free
	dnFdWR_cs = malloc(Nn*Nvar*Neq * sizeof *dnFdWR_cs); // free

	W_ptr = W;
	for (var = 0; var < Nvar; var++) {
		IndWLR = var*Nn;
		for (n = 0; n < Nn; n++)
			WL[IndWLR+n] = *W_ptr++;
		for (n = 0; n < Nn; n++)
			WR[IndWLR+n] = *W_ptr++;
	}

	if (strstr(nFType,"LF")) {
		jacobian_flux_LF(Nn,1,WL,WR,dnFdWL,nL,d,Neq,'L');
		jacobian_flux_LF(Nn,1,WL,WR,dnFdWR,nL,d,Neq,'R');
	} else if (strstr(nFType,"Roe")) {
		printf("Add support Roe.\n"), EXIT_MSG;
	} else {
		printf("Error: Unsupported nFType.\n"), EXIT_MSG;
	}
	compute_dnFdW_cs(Nn,1,d,Neq,WL,WR,dnFdWL_cs,dnFdWR_cs,nL,nFType);

//array_print_d(Nn*Nvar,Neq,dnFdWL,'C');
//array_print_d(Nn*Nvar,Neq,dnFdWL_cs,'C');

	if (strstr(nFType,"LF")) {
		CheckedAllLF = 1;
		for (i = 0; i < 2; i++) {
			if (!TestDB.EnteredLF[i]) {
				CheckedAllLF = 0;
				break;
			}
		}
		if (CheckedAllLF &&
		    array_norm_diff_d(Nn*Nvar*Neq,dnFdWL,dnFdWL_cs,"Inf") < EPS &&
		    array_norm_diff_d(Nn*Nvar*Neq,dnFdWR,dnFdWR_cs,"Inf") < EPS)
				pass = 1, TestDB.Npass++;
	} else {
		if (array_norm_diff_d(Nn*Nvar*Neq,dnFdWL,dnFdWL_cs,"Inf") < EPS &&
		    array_norm_diff_d(Nn*Nvar*Neq,dnFdWR,dnFdWR_cs,"Inf") < EPS)
				pass = 1, TestDB.Npass++;
	}

	free(WL);
	free(WR);
	free(dnFdWL);
	free(dnFdWR);
	free(dnFdWL_cs);
	free(dnFdWR_cs);

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
	double       *W, *nL;

	for (d = 1; d <= 3; d++) {
		Neq = d+2;

		W  = initialize_W(&Nn,&Nel,d); // free
		nL = initialize_n(Nn,Nel,d);   // free

		// flux_inviscid
		pass = compare_jacobian_flux_inviscid(Nn,Nel,d,Neq,W);
		if (d == 1) printf("jacobian_flux_inviscid (d = %d):                  ",d);
		else        printf("         flux_inviscid (d = %d):                  ",d);
		test_print(pass);

		// flux_LF
		pass = compare_jacobian_flux_Num(Nn,Nel,d,Neq,W,nL,"LF");
		printf("         flux_LF              :                  ");
		test_print(pass);

		free(W);
		free(nL);
	}
}
