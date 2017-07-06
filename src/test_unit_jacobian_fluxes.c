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
#include "fluxes_structs.h"
#include "fluxes_inviscid_c.h"
#include "fluxes_viscous_c.h"
#include "jacobian_fluxes_inviscid.h"
#include "jacobian_fluxes_viscous.h"
#include "initialize_test_case.h"
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
	unsigned int const Neq     = DB.Neq,
	                   Nvar    = DB.Nvar,
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

	struct S_FLUX *const FLUXDATA = malloc(sizeof *FLUXDATA); // free
	FLUXDATA->PDE_index = DB.PDE_index;
	FLUXDATA->d   = d;
	FLUXDATA->Nn  = Nn;
	FLUXDATA->Nel = Nel;
	FLUXDATA->W_c = Wp;
	FLUXDATA->Q_c = (double complex const *const *const) Q_c;
	FLUXDATA->F_c = F;

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
			flux_inviscid_c(FLUXDATA);
		else if (flux_type == 'V')
			flux_viscous_c(FLUXDATA);
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
	free(FLUXDATA);
}

static void compute_dFdQ_cs(unsigned int const Nn, unsigned int const Nel, unsigned int const d, double const *const W,
                            double const *const *const Q, double *const *const dFdQ_cs)
{
	unsigned int const Neq     = DB.Neq,
	                   Nvar    = DB.Nvar,
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

	struct S_FLUX *const FLUXDATA = malloc(sizeof *FLUXDATA); // free
	FLUXDATA->PDE_index = DB.PDE_index;
	FLUXDATA->d   = d;
	FLUXDATA->Nn  = Nn;
	FLUXDATA->Nel = Nel;
	FLUXDATA->W_c = W_c;
	FLUXDATA->Q_c = (double complex const *const *const) Qp;
	FLUXDATA->F_c = F;

	for (size_t dim1 = 0; dim1 < d; dim1++) {
		for (size_t var = 0; var < Nvar; var++) {
			double const h = EPS*EPS;

			size_t const IndQ = NnTotal*var;
			for (size_t n = 0; n < NnTotal; n++)
				Qp[dim1][IndQ+n] += h*I;

			flux_viscous_c(FLUXDATA);

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
	free(FLUXDATA);
}

static unsigned int compare_jacobian_flux_inviscid(const unsigned int Nn, const unsigned int Nel, double *W)
{
	unsigned int pass = 0;

	unsigned int const d       = DB.d,
	                   Neq     = DB.Neq,
	                   Nvar    = DB.Nvar,
	                   NnTotal = Nn*Nel;

	double *const dFdW    = malloc(NnTotal*d*Nvar*Neq * sizeof *dFdW);    // free
	double *const dFdW_cs = malloc(NnTotal*d*Nvar*Neq * sizeof *dFdW_cs); // free

	struct S_FLUX *FLUXDATA = malloc(sizeof *FLUXDATA); // free
	FLUXDATA->PDE_index = DB.PDE_index;
	FLUXDATA->d    = d;
	FLUXDATA->Nn   = Nn;
	FLUXDATA->Nel  = Nel;
	FLUXDATA->W    = W;
	FLUXDATA->F    = NULL;
	FLUXDATA->dFdW = dFdW;

	jacobian_flux_inviscid(FLUXDATA);
	compute_dFdW_cs(Nn,Nel,d,W,NULL,dFdW_cs,'I');
	free(FLUXDATA);

	if (array_norm_diff_d(NnTotal*d*Nvar*Neq,dFdW,dFdW_cs,"Inf") < EPS) {
		pass = 1;
	} else {
		array_print_d(NnTotal,d*Nvar*Neq,dFdW,'C');
		array_print_d(NnTotal,d*Nvar*Neq,dFdW_cs,'C');
	}

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

	struct S_FLUX *const FLUXDATA = malloc(sizeof *FLUXDATA); // free
	FLUXDATA->PDE_index = DB.PDE_index;
	FLUXDATA->d   = d;
	FLUXDATA->Nn  = Nn;
	FLUXDATA->Nel = Nel;
	FLUXDATA->W   = W;
	FLUXDATA->Q   = (double const *const *const) Q;
	FLUXDATA->F   = NULL;

	if (pert_var == 'W') {
		double *const dFdW    = malloc(NnTotal*d*Nvar*Neq * sizeof *dFdW),    // free
			   *const dFdW_cs = malloc(NnTotal*d*Nvar*Neq * sizeof *dFdW_cs); // free

		FLUXDATA->dFdW = dFdW;
		FLUXDATA->dFdQ = NULL;
		jacobian_flux_viscous(FLUXDATA);
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

		FLUXDATA->dFdW = NULL;
		FLUXDATA->dFdQ = dFdQ;
		jacobian_flux_viscous(FLUXDATA);
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
	free(FLUXDATA);

	return pass;
}

static void compute_dnFdW_cs(unsigned int const Nn, unsigned int const Nel, double const *const nL,
                             double const *const XYZ, double *WL, double *WR, double *dnFdWL_cs, double *dnFdWR_cs)
{
	unsigned int const d       = DB.d,
	                   Neq     = DB.Neq,
	                   Nvar    = DB.Nvar,
	                   NnTotal = Nn*Nel;


	double const h = EPS*EPS;

	double complex WLp[Nn*Nel*Neq], WRp[Nn*Nel*Neq];
	for (size_t i = 0, iMax = NnTotal*Nvar; i < iMax; i++) {
		WLp[i] = WL[i];
		WRp[i] = WR[i];
	}

	double complex *const nF = malloc(NnTotal*Neq * sizeof *nF); // free

	struct S_NUMERICALFLUX *const NUMFLUXDATA = malloc(sizeof *NUMFLUXDATA); // free

	NUMFLUXDATA->NumFluxInviscid_index = DB.InviscidFluxType;
	NUMFLUXDATA->d   = d;
	NUMFLUXDATA->Nn  = Nn;
	NUMFLUXDATA->Nel = Nel;
	NUMFLUXDATA->nL  = nL;
	NUMFLUXDATA->XYZ = XYZ;

	NUMFLUXDATA->WL_c       = WLp;
	NUMFLUXDATA->WR_c       = WRp;
	NUMFLUXDATA->nFluxNum_c = nF;

	for (size_t var = 0; var < Nvar; var++) {
		// Left
		size_t const IndW = NnTotal*var;
		for (size_t n = 0; n < NnTotal; n++)
			WLp[IndW+n] += h*I;

		flux_num_inviscid_c(NUMFLUXDATA);

		for (size_t eq = 0; eq < Neq; eq++) {
			size_t const IndnF    = eq*NnTotal,
			             InddnFdW = (eq*Neq+var)*NnTotal;
			for (size_t n = 0; n < NnTotal; n++)
				dnFdWL_cs[InddnFdW+n] = cimag(nF[IndnF+n])/h;
		}

		// Right
		for (size_t n = 0; n < NnTotal; n++) {
			WLp[IndW+n] -= h*I;
			WRp[IndW+n] += h*I;
		}

		flux_num_inviscid_c(NUMFLUXDATA);

		for (size_t eq = 0; eq < Neq; eq++) {
			size_t const IndnF    = eq*NnTotal,
			             InddnFdW = (eq*Neq+var)*NnTotal;
			for (size_t n = 0; n < NnTotal; n++)
				dnFdWR_cs[InddnFdW+n] = cimag(nF[IndnF+n])/h;
		}

		for (size_t n = 0; n < NnTotal; n++)
			WRp[IndW+n] -= h*I;
	}
	free(nF);
	free(NUMFLUXDATA);
}

static unsigned int compare_jacobian_flux_Num(unsigned int const Nn, unsigned int const Nel, double const *const nL,
                                              double const *const XYZ, double const *const W)
{
	if (Nel != 2)
		EXIT_UNSUPPORTED;

	unsigned int pass = 0;

	unsigned int const d    = DB.d,
	                   Nvar = DB.Nvar,
	                   Neq  = DB.Neq;

	unsigned int i, IndWLR, CheckedAllLF, CheckedAllRoe;

	double *const WL        = malloc(Nn*Nvar     * sizeof *WL),        // free
	       *const WR        = malloc(Nn*Nvar     * sizeof *WR),        // free
	       *const dnFdWL    = malloc(Nn*Nvar*Neq * sizeof *dnFdWL),    // free
	       *const dnFdWR    = malloc(Nn*Nvar*Neq * sizeof *dnFdWR),    // free
	       *const dnFdWL_cs = malloc(Nn*Nvar*Neq * sizeof *dnFdWL_cs), // free
	       *const dnFdWR_cs = malloc(Nn*Nvar*Neq * sizeof *dnFdWR_cs); // free

	double const *W_ptr = W;
	for (size_t var = 0; var < Nvar; var++) {
		IndWLR = var*Nn;
		for (size_t n = 0; n < Nn; n++)
			WL[IndWLR+n] = *W_ptr++;
		for (size_t n = 0; n < Nn; n++)
			WR[IndWLR+n] = *W_ptr++;
	}

	struct S_NUMERICALFLUX *const NUMFLUXDATA = malloc(sizeof *NUMFLUXDATA); // free

	NUMFLUXDATA->NumFluxInviscid_index = DB.InviscidFluxType;
	NUMFLUXDATA->d   = d;
	NUMFLUXDATA->Nn  = Nn;
	NUMFLUXDATA->Nel = 1;
	NUMFLUXDATA->nL  = nL;
	NUMFLUXDATA->XYZ = XYZ;

	NUMFLUXDATA->WL           = WL;
	NUMFLUXDATA->WR           = WR;
	NUMFLUXDATA->nFluxNum     = NULL;
	NUMFLUXDATA->dnFluxNumdWL = dnFdWL;
	NUMFLUXDATA->dnFluxNumdWR = dnFdWR;

	jacobian_flux_num_inviscid(NUMFLUXDATA);
	free(NUMFLUXDATA);

	compute_dnFdW_cs(Nn,1,nL,XYZ,WL,WR,dnFdWL_cs,dnFdWR_cs);

	if (DB.InviscidFluxType == FLUX_LF) {
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
	} else if (DB.InviscidFluxType ==  FLUX_ROE) {
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
	 *	Notation:
	 *		FiType : (F)lux (i)nviscid (Type)
	 *
	 */

	unsigned int Nn, Nel, d;
	double       *W, *nL;

	set_memory_test_jacobians('A');

	unsigned int NFiTypes = 0;
	char         **FiType;
	set_FiTypes(&NFiTypes,&FiType); // free

	unsigned int NFNumTypes = 0;
	char         **FNumType;
	set_FNumTypes(&NFNumTypes,&FNumType); // free

	strcpy(DB.MeshType,"ToBeCurved"); // Meshes are not used for this test (and thus need not exist)

	test_print_warning("Advection test case currently used has b = constant");

	unsigned int dMin = 2, dMax = 3;
	for (d = dMin; d <= dMax; d++) {
		W   = initialize_W(&Nn,&Nel,d); // free
		nL  = initialize_n(Nn,Nel,d);   // free
		double const *const XYZ = initialize_XYZ(Nn,Nel,d); // free

		// flux_inviscid
		for (size_t i = 0; i < NFiTypes; i++) {
			set_parameters_test_flux_inviscid(FiType[i],d);
			initialize_test_case_parameters();

			pass = compare_jacobian_flux_inviscid(Nn,Nel,W);
			if (i == 0) {
				if (d == dMin)
					sprintf(PrintName,"jacobian_flux_%s (d = %d):",FiType[i],d);
				else
					sprintf(PrintName,"              %s (d = %d):",FiType[i],d);
			} else {
				sprintf(PrintName,"              %s        :",FiType[i]);
			}
			test_print2(pass,PrintName);

			free(DB.SolverType); // Initialized in "initialize_test_case_parameters"
		}

		// flux_num_inviscid
		for (size_t i = 0; i < NFNumTypes; i++) {
			set_parameters_test_flux_Num(FNumType[i],d);
			initialize_test_case_parameters();

			reset_entered_test_num_flux(FNumType[i]);
			pass = compare_jacobian_flux_Num(Nn,Nel,nL,XYZ,W);
			if (i == 0)
				sprintf(PrintName,"              Num_%s      :",FNumType[i]);
			else
				sprintf(PrintName,"                  %s      :",FNumType[i]);
			test_print2(pass,PrintName);

			free(DB.SolverType); // Initialized in "initialize_test_case_parameters"
		}

		// flux_viscous
		DB.Neq = DB.Nvar = d+2;
		DB.PDE_index = PDE_NAVIERSTOKES;

		double **Q = initialize_Q(Nn,Nel,d); // free
		pass = compare_jacobian_flux_viscous(Nn,Nel,d,W,Q,'W');
		sprintf(PrintName,"              viscous (W)      :");
		test_print2(pass,PrintName);

		pass = compare_jacobian_flux_viscous(Nn,Nel,d,W,Q,'Q');
		sprintf(PrintName,"                      (Q)      :");
		test_print2(pass,PrintName);

		free(W);
		free(nL);
		free((double *) XYZ);
		array_free2_d(d,Q);
	}

	set_memory_test_jacobians('F');

	array_free2_c(NFiTypes,FiType);
	array_free2_c(NFNumTypes,FNumType);

	free(PrintName);
}
