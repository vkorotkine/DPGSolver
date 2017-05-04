// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "test_unit_jacobian_fluxes.h"

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <stdbool.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "Test.h"

#include "test_code_fluxes.h"
#include "test_support.h"
#include "fluxes_inviscid_c.h"
#include "fluxes_viscous_c.h"
#include "jacobian_fluxes_inviscid.h"
#include "jacobian_fluxes_viscous.h"
#include "array_norm.h"
#include "array_free.h"
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

static void compute_dFdW_cs(unsigned int const Nn, unsigned int const Nel, unsigned int const d, double const *const W,
                            double const *const *const Q, double *const dFdW_cs, char const flux_type)
{
	unsigned int const Neq     = d+2,
	                   Nvar    = d+2,
	                   NnTotal = Nn*Nel;

	double complex **Q_c = NULL;
	if (Q != NULL) {
		Q_c = malloc(d * sizeof *Q_c); // free
		for (size_t dim = 0; dim < d; dim++) {
			Q_c[dim] = malloc(NnTotal*Nvar * sizeof *Q_c[dim]); // free
			for (size_t i = 0; i < NnTotal*Nvar; i++)
				Q_c[dim][i] = Q[dim][i];
		}
	}

	double complex *const Wp = malloc(NnTotal*Neq  * sizeof *Wp), // free
	               *const F = malloc(NnTotal*d*Neq * sizeof *F);  // free

	for (size_t var = 0; var < Nvar; var++) {
		double const h = EPS*EPS;
		for (size_t eq = 0; eq < Neq; eq++) {
			size_t const IndW = NnTotal*eq;
			for (size_t n = 0; n < NnTotal; n++) {
				Wp[IndW+n] = W[IndW+n];
				if (eq == var)
					Wp[IndW+n] += h*I;
			}
		}
		if (flux_type == 'I')
			flux_inviscid_c(Nn,Nel,Wp,F,d,Neq);
		else if (flux_type == 'V')
			flux_viscous_c(Nn,Nel,Wp,(double complex const *const *const) Q_c,F);
		else
			EXIT_UNSUPPORTED;

		for (size_t eq = 0; eq < Neq; eq++) {
		for (size_t dim = 0; dim < d; dim++) {
			size_t const IndF    = (eq*d+dim)*NnTotal,
			             InddFdW = ((eq*Neq+var)*d+dim)*NnTotal;
			for (size_t n = 0; n < NnTotal; n++)
				dFdW_cs[InddFdW+n] = cimag(F[IndF+n])/h;
		}}
	}

	if (Q != NULL)
		array_free2_cmplx(d,Q_c);

	free(Wp);
	free(F);
}

static void compute_dFdQ_cs(unsigned int const Nn, unsigned int const Nel, unsigned int const d, double const *const W,
                            double const *const *const Q, double *const *const dFdQ_cs)
{
	unsigned int const Neq     = d+2,
	                   Nvar    = d+2,
	                   NnTotal = Nn*Nel;

	double complex *const W_c = malloc(NnTotal*Neq * sizeof *W_c); // free
	for (size_t i = 0; i < NnTotal*Nvar; i++)
		W_c[i] = W[i];

	double complex **const Qp = malloc(d * sizeof *Qp); // free
	for (size_t dim = 0; dim < d; dim++) {
		Qp[dim] = malloc(NnTotal*Nvar * sizeof *Qp[dim]); // free
		for (size_t i = 0; i < NnTotal*Nvar; i++)
			Qp[dim][i] = Q[dim][i];
	}

	double complex *const F = malloc(NnTotal*d*Neq * sizeof *F);  // free
	for (size_t dim1 = 0; dim1 < d; dim1++) {
		for (size_t var = 0; var < Nvar; var++) {
			double const h = EPS*EPS;

			size_t const IndQ = NnTotal*var;
			for (size_t n = 0; n < NnTotal; n++)
				Qp[dim1][IndQ+n] += h*I;

			flux_viscous_c(Nn,Nel,W_c,(double complex const *const *const) Qp,F);

			for (size_t eq = 0; eq < Neq; eq++) {
			for (size_t dim = 0; dim < d; dim++) {
				size_t const IndF    = (eq*d+dim)*NnTotal,
							 InddFdQ = ((eq*Neq+var)*d+dim)*NnTotal;
				for (size_t n = 0; n < NnTotal; n++)
					dFdQ_cs[dim1][InddFdQ+n] = cimag(F[IndF+n])/h;
			}}

			for (size_t n = 0; n < NnTotal; n++)
				Qp[dim1][IndQ+n] -= h*I;
		}
	}

	if (Q != NULL)
		array_free2_cmplx(d,Qp);

	free(W_c);
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
	compute_dFdW_cs(Nn,Nel,d,W,NULL,dFdW_cs,'I');

	if (array_norm_diff_d(NnTotal*d*Nvar*Neq,dFdW,dFdW_cs,"Inf") < EPS)
		pass = 1;

	free(dFdW);
	free(dFdW_cs);

	return pass;
}

static unsigned int compare_jacobian_flux_viscous(unsigned int const Nn, unsigned int const Nel, unsigned int const d,
                                                  double *const W, double *const *const Q, char const pert_var)
{
	DB.d        = d;
	DB.Pr       = 0.72;
	DB.mu       = 1e-0;
	DB.Const_mu = 1;

	unsigned int pass = 0;

	unsigned int const NnTotal = Nn*Nel,
	                   Neq     = d+2,
	                   Nvar    = d+2;

	if (pert_var == 'W') {
		double *const dFdW    = malloc(NnTotal*d*Nvar*Neq * sizeof *dFdW),    // free
			   *const dFdW_cs = malloc(NnTotal*d*Nvar*Neq * sizeof *dFdW_cs); // free

		jacobian_flux_viscous(Nn,Nel,W,(double const *const *const) Q,dFdW,NULL);
		compute_dFdW_cs(Nn,Nel,d,W,(double const *const *const) Q,dFdW_cs,'V');

		if (array_norm_diff_d(NnTotal*d*Nvar*Neq,dFdW,dFdW_cs,"Inf") < EPS)
			pass = 1;

		free(dFdW);
		free(dFdW_cs);
	} else if (pert_var == 'Q') {
		double **const dFdQ    = malloc(d * sizeof *dFdQ),    // free
		       **const dFdQ_cs = malloc(d * sizeof *dFdQ_cs); // free
		for (size_t dim = 0; dim < d; dim++) {
			dFdQ[dim]    = malloc(NnTotal*d*Nvar*Neq * sizeof *dFdQ[dim]);    // free
			dFdQ_cs[dim] = malloc(NnTotal*d*Nvar*Neq * sizeof *dFdQ_cs[dim]); // free
		}

		jacobian_flux_viscous(Nn,Nel,W,(double const *const *const) Q,NULL,dFdQ);
		compute_dFdQ_cs(Nn,Nel,d,W,(double const *const *const) Q,dFdQ_cs);

		double diff = 0.0;
		for (size_t dim = 0; dim < d; dim++)
			diff += array_norm_diff_d(NnTotal*d*Nvar*Neq,dFdQ[dim],dFdQ_cs[dim],"Inf");

		if (diff < 1e1*EPS) {
			pass = 1;
		} else {
			bool PrintEnabled = 1;
			if (PrintEnabled) {
				printf("% .3e\n",diff);
				for (size_t dim = 0; dim < d; dim++) {
					printf("dFdQ, dim: %zu\n",dim);
					array_print_d(NnTotal*d*Nvar,Neq,dFdQ[dim],'C');
					array_print_d(NnTotal*d*Nvar,Neq,dFdQ_cs[dim],'C');
				}
			}
		}

		array_free2_d(d,dFdQ);
		array_free2_d(d,dFdQ_cs);
	} else {
		EXIT_UNSUPPORTED;
	}

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
			flux_Roe_c(Nn,Nel,WLp,WRp,nF,nL,d,Neq);
		else
			EXIT_UNSUPPORTED;

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
			flux_Roe_c(Nn,Nel,WLp,WRp,nF,nL,d,Neq);
		else
			EXIT_UNSUPPORTED;

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

	unsigned int i, n, var, Nvar, IndWLR, CheckedAllLF, CheckedAllRoe;
	double       *W_ptr, *WL, *WR, *dnFdWL, *dnFdWR, *dnFdWL_cs, *dnFdWR_cs;

	if (Nel != 2)
		EXIT_UNSUPPORTED;

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
		jacobian_flux_Roe(Nn,1,WL,WR,dnFdWL,nL,d,Neq,'L');
		jacobian_flux_Roe(Nn,1,WL,WR,dnFdWR,nL,d,Neq,'R');
	} else {
		EXIT_UNSUPPORTED;
	}
	compute_dnFdW_cs(Nn,1,d,Neq,WL,WR,dnFdWL_cs,dnFdWR_cs,nL,nFType);

	if (strstr(nFType,"LF")) {
		CheckedAllLF = 1;
		for (i = 0; i < TEST_N_LF; i++) {
			if (!TestDB.EnteredLF[i]) {
				CheckedAllLF = 0;
				break;
			}
		}
		if (CheckedAllLF &&
		    array_norm_diff_d(Nn*Nvar*Neq,dnFdWL,dnFdWL_cs,"Inf") < EPS &&
		    array_norm_diff_d(Nn*Nvar*Neq,dnFdWR,dnFdWR_cs,"Inf") < EPS)
				pass = 1;
	} else if (strstr(nFType,"Roe")) {
		CheckedAllRoe = 1;
		for (i = 0; i < TEST_N_ROE; i++) {
			if (!TestDB.EnteredRoe[i]) {
				CheckedAllRoe = 0;
				break;
			}
		}
		if (CheckedAllRoe &&
		    array_norm_diff_d(Nn*Nvar*Neq,dnFdWL,dnFdWL_cs,"Inf") < EPS &&
		    array_norm_diff_d(Nn*Nvar*Neq,dnFdWR,dnFdWR_cs,"Inf") < EPS)
				pass = 1;
	} else {
		if (array_norm_diff_d(Nn*Nvar*Neq,dnFdWL,dnFdWL_cs,"Inf") < EPS &&
		    array_norm_diff_d(Nn*Nvar*Neq,dnFdWR,dnFdWR_cs,"Inf") < EPS)
				pass = 1;
	}

	free(WL);
	free(WR);
	free(dnFdWL);
	free(dnFdWR);
	free(dnFdWL_cs);
	free(dnFdWR_cs);

	return pass;
}

void test_unit_jacobian_fluxes(void)
{
	unsigned int pass;

	char *PrintName = malloc(STRLEN_MAX * sizeof *PrintName); // free

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

	unsigned int dMin = 2, dMax = 3;
	for (d = dMin; d <= dMax; d++) {
		Neq = d+2;

		W  = initialize_W(&Nn,&Nel,d); // free
		nL = initialize_n(Nn,Nel,d);   // free

		// flux_inviscid
		pass = compare_jacobian_flux_inviscid(Nn,Nel,d,Neq,W);
		if (d == dMin)
			sprintf(PrintName,"jacobian_flux_inviscid (d = %d):",d);
		else
			sprintf(PrintName,"         flux_inviscid (d = %d):",d);
		test_print2(pass,PrintName);


		// flux_LF
		for (size_t i = 0; i < TEST_N_LF; i++)
			TestDB.EnteredLF[i] = 0;
		pass = compare_jacobian_flux_Num(Nn,Nel,d,Neq,W,nL,"LF");
		sprintf(PrintName,"         flux_LF              :");
		test_print2(pass,PrintName);


		// flux_Roe
		for (size_t i = 0; i < TEST_N_ROE; i++)
			TestDB.EnteredRoe[i] = 0;
		pass = compare_jacobian_flux_Num(Nn,Nel,d,Neq,W,nL,"Roe");
		sprintf(PrintName,"         flux_Roe             :");
		test_print2(pass,PrintName);


		// flux_viscous
		double **Q = initialize_Q(Nn,Nel,d); // free
		pass = compare_jacobian_flux_viscous(Nn,Nel,d,W,Q,'W');
		sprintf(PrintName,"         flux_viscous (W)     :");
		test_print2(pass,PrintName);

		pass = compare_jacobian_flux_viscous(Nn,Nel,d,W,Q,'Q');
		sprintf(PrintName,"                      (Q)     :");
		test_print2(pass,PrintName);

		free(W);
		free(nL);
		array_free2_d(d,Q);
	}

	free(PrintName);
}
