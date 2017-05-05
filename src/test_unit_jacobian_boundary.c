// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "test_unit_jacobian_boundary.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <complex.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "Test.h"

#include "test_code_fluxes.h"
#include "test_code_boundary_conditions.h"
#include "test_support.h"
#include "boundary_conditions.h"
#include "jacobian_boundary_conditions.h"
#include "boundary_conditions_c.h"
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

static void compute_dWdW_cs(const unsigned int Neq, const unsigned int Nn, const unsigned int Nel, const unsigned int d,
                            double *W, double *dWdW_cs, double *nL, double *XYZ, const char *BType)
{
	unsigned int   i, iMax, n, var, var2, Nvar, NnTotal, IndW, IndWB, InddWdW;
	double         h;
	double complex Wp[Nn*Nel*Neq], *WB;

	h = EPS*EPS;

	Nvar    = Neq;
	NnTotal = Nn*Nel;

	WB = malloc(NnTotal*Nvar * sizeof *WB); // free

	for (i = 0, iMax = NnTotal*Nvar; i < iMax; i++)
		Wp[i] = W[i];

	for (var = 0; var < Nvar; var++) {
		IndW = NnTotal*var;
		for (n = 0; n < NnTotal; n++)
			Wp[IndW+n] += h*I;

		struct S_BC *const BCdata = malloc(sizeof *BCdata); // free

		BCdata->d = d;
		BCdata->Nn = Nn;
		BCdata->Nel = Nel;

		BCdata->XYZ  = XYZ;
		BCdata->nL   = nL;

		BCdata->WL_c = Wp;
		BCdata->QL_c = NULL;

		BCdata->WB_c = WB;
		BCdata->QB_c = NULL;

		set_BC_from_BType(BCdata,BType);
		correct_XYZ_for_exact_normal(BCdata,BType);
		compute_boundary_values_c(BCdata);
		free(BCdata);

		for (var2 = 0; var2 < Nvar; var2++) {
			IndWB = NnTotal*var2;
			InddWdW = (var*Nvar+var2)*NnTotal;
			for (n = 0; n < NnTotal; n++)
				dWdW_cs[InddWdW+n] = cimag(WB[IndWB+n])/h;
		}

		for (n = 0; n < NnTotal; n++)
			Wp[IndW+n] -= h*I;
	}
	free(WB);
}

static unsigned int compare_jacobian_boundary(const unsigned int Nn, const unsigned int Nel, const unsigned int d,
                                              double *W, double **Q, double *nL, double *XYZ, const char *BType)
{
	unsigned int pass = 0;

	unsigned int const Nvar    = DB.Nvar,
	                   Neq     = DB.Neq,
	                   NnTotal = Nn*Nel;

	double *dWdW, *dWdW_cs;

	dWdW    = malloc(NnTotal*Nvar*Neq * sizeof *dWdW);    // free
	dWdW_cs = malloc(NnTotal*Nvar*Neq * sizeof *dWdW_cs); // free

	struct S_BC *const BCdata = malloc(sizeof *BCdata); // free

	BCdata->d = d;
	BCdata->Nn = Nn;
	BCdata->Nel = Nel;

	BCdata->XYZ = XYZ;
	BCdata->nL  = nL;
	BCdata->WL  = W;
	BCdata->QL  = (double const *const *const) Q;

	BCdata->dWBdWL = dWdW;
	BCdata->dQBdWL = NULL;

	set_BC_from_BType(BCdata,BType);
	correct_XYZ_for_exact_normal(BCdata,BType);
	compute_jacobian_boundary_values(BCdata);

	compute_dWdW_cs(Neq,Nn,Nel,d,W,dWdW_cs,nL,XYZ,BType);

	free(BCdata);

	if (array_norm_diff_d(NnTotal*Nvar*Neq,dWdW,dWdW_cs,"Inf") < 10*EPS) {
		pass = 1;
	} else {
		array_print_d(NnTotal*Nvar,Neq,dWdW,'C');
		array_print_d(NnTotal*Nvar,Neq,dWdW_cs,'C');
		printf("%d %d %d % .3e\n",NnTotal,Nvar,DB.Nvar,array_norm_diff_d(NnTotal*Nvar*Neq,dWdW,dWdW_cs,"Inf"));
	}

	bool CheckedAll = 0;
	check_entered_test_boundary_conditions(&CheckedAll,BType);

	if (!CheckedAll) {
		pass = 0;
		printf("Did not check all boundary condition settings for boundary %s\n",BType);
	}

	free(dWdW);
	free(dWdW_cs);

	return pass;
}

void test_unit_jacobian_boundary(void)
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

	test_print_warning("Not testing dQBdW (Currently unused in the code)");

	unsigned int i, Nn, Nel, d, Neq;
	double       *W, **Q, *nL, *XYZ;

	set_memory_test_boundary_conditions('a');

	unsigned int NBTypes = 0;
	char         **BType;
	set_BTypes(&NBTypes,&BType); // free

	strcpy(DB.MeshType,"ToBeCurved"); // Meshes are not used for this test (and thus need not exist)

	unsigned int const dMin = 2, dMax = 3;
	for (d = dMin; d <= dMax; d++) {
	for (i = 0; i < NBTypes; i++) {
		set_parameters_test_boundary_conditions(BType[i],d);
		initialize_test_case_parameters();

		reset_entered_test_boundary_conditions(BType[i]);

		Neq = d+2;

		W    = initialize_W(&Nn,&Nel,d); // free
		Q    = initialize_Q(Nn,Nel,d);   // free
		nL   = initialize_n(Nn,Nel,d);   // free
		XYZ  = initialize_XYZ(Nn,Nel,d); // free

		if (strstr(BType[i],"BackPressure"))
			update_values_BackPressure(Nn,Nel,W,nL,d);
		pass = compare_jacobian_boundary(Nn,Nel,d,W,Q,nL,XYZ,BType[i]);

		if (i == 0) {
			if (d == dMin)
				sprintf(PrintName,"jacobian_boundary_%s (d = %d):",BType[i],d);
			else
				sprintf(PrintName,"                  %s (d = %d):",BType[i],d);
		} else {
			sprintf(PrintName,"                  %s        :",BType[i]);
		}
		test_print2(pass,PrintName);

		free(W);
		array_free2_d(d,Q);
		free(nL);
		free(XYZ);

		free(DB.SolverType); // Initialized in "initialize_test_case_parameters"
	}}

	set_memory_test_boundary_conditions('f');

	for (i = 0; i < NBTypes; i++)
		free(BType[i]);

	free(PrintName);
}
