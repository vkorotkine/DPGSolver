// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "test_unit_jacobian_boundary.h"

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "Test.h"

#include "test_code_fluxes_inviscid.h"
#include "test_support.h"
#include "array_norm.h"
#include "jacobian_boundary_conditions.h"
#include "boundary_conditions_c.h"
#include "initialize_test_case.h"

#include "array_print.h"
#include "boundary_conditions.h" // ToBeDeleted

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

		if (strstr(BType,"SlipWall"))
			boundary_SlipWall_c(Nn,Nel,Wp,WB,nL,d);
		else if (strstr(BType,"Riemann"))
			boundary_Riemann_c(Nn,Nel,XYZ,Wp,NULL,WB,nL,d);
		else if (strstr(BType,"BackPressure"))
			boundary_BackPressure_c(Nn,Nel,Wp,WB,nL,d,Neq);
		else
			printf("Error: Unsupported BType.\n"), EXIT_MSG;

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
                                              const unsigned int Neq, double *W, double *nL, double *XYZ,
                                              const char *BType)
{
	unsigned int pass = 0;

	unsigned int NnTotal, Nvar, i, CheckedAll;
	double       *dWdW, *dWdW_cs;

	NnTotal = Nn*Nel;
	Nvar    = Neq;

	dWdW    = malloc(NnTotal*Nvar*Neq * sizeof *dWdW);    // free
	dWdW_cs = malloc(NnTotal*Nvar*Neq * sizeof *dWdW_cs); // free

	if (strstr(BType,"SlipWall"))
		jacobian_boundary_SlipWall(Nn,Nel,W,dWdW,nL,d,Neq);
	else if (strstr(BType,"Riemann"))
		jacobian_boundary_Riemann(Nn,Nel,XYZ,W,NULL,dWdW,nL,d,Neq);
	else if (strstr(BType,"BackPressure"))
		jacobian_boundary_BackPressure(Nn,Nel,W,dWdW,nL,d,Neq);
	else
		printf("Error: Unsupported BType.\n"), EXIT_MSG;

	compute_dWdW_cs(Neq,Nn,Nel,d,W,dWdW_cs,nL,XYZ,BType);

	if (strstr(BType,"Riemann")) {
		CheckedAll = 1;
		for (i = 0; i < 4; i++) {
			if (!TestDB.EnteredRiemann[i]) {
				CheckedAll = 0;
				break;
			}
		}
//		array_print_ui(1,4,TestDB.EnteredRiemann,'R');
		if (CheckedAll && array_norm_diff_d(NnTotal*Nvar*Neq,dWdW,dWdW_cs,"Inf") < 10*EPS)
			pass = 1, TestDB.Npass++;
	} else if (strstr(BType,"BackPressure")) {
		CheckedAll = 1;
		for (i = 0; i < 2; i++) {
			if (!TestDB.EnteredBackPressure[i]) {
				CheckedAll = 0;
				break;
			}
		}
		array_print_ui(1,2,TestDB.EnteredBackPressure,'R');
		if (CheckedAll && array_norm_diff_d(NnTotal*Nvar*Neq,dWdW,dWdW_cs,"Inf") < 10*EPS)
			pass = 1, TestDB.Npass++;
		else {
			printf("%d % .3e\n",CheckedAll,array_norm_diff_d(NnTotal*Nvar*Neq,dWdW,dWdW_cs,"Inf"));
			array_print_d(Neq*Nvar,NnTotal,dWdW,'R');
			array_print_d(Neq*Nvar,NnTotal,dWdW_cs,'R');
		}
	} else {
		if (array_norm_diff_d(NnTotal*Nvar*Neq,dWdW,dWdW_cs,"Inf") < EPS)
			pass = 1, TestDB.Npass++;
	}


	free(dWdW);
	free(dWdW_cs);

	return pass;
}

static void update_values(const unsigned int Nn, const unsigned int Nel, double *W, double *nL, const unsigned int d)
{
	unsigned int n, NnTotal, dim, var;

	NnTotal = Nn*Nel;

	if (NnTotal != 6)
		printf("Error: Unsupported.\n"), EXIT_MSG;

	if (d == 3) {
		unsigned int FinalIndices[6] = {0,1,2,3,4,5};
		for (n = 0; n < NnTotal; n++) {
			if (FinalIndices[n] != n) {
				for (dim = 0; dim < d; dim++)
					nL[n*d+dim] *= -1.0;
				if (FinalIndices[n] == n+1) {
					for (var = 1; var < d+1; var++)
						W[n+var*NnTotal] *= 0.5;
				}
			}
		}
	} else if (d == 2) {
		unsigned int FinalIndices[6] = {0,2,0,4,0,6};
//		unsigned int FinalIndices[6] = {0,0,0,0,0,0};
		for (n = 0; n < NnTotal; n++) {
			if (FinalIndices[n] != n) {
				for (dim = 0; dim < d; dim++)
					nL[n*d+dim] *= -1.0;
				if (FinalIndices[n] == n+1) {
					for (var = 1; var < d+1; var++)
						W[n+var*NnTotal] *= 0.5;
				}
			}
		}
	} else {
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}

}

void test_unit_jacobian_boundary(void)
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

printf("\nWarning: boundary_Riemann is currently not being tested for d = 1.\n\n"); TestDB.Nwarnings++;
// This requires a case where a d = 1 Riemann BC is supported.

	unsigned int NBTypes = 3;

	char         *BType[NBTypes], *TestCase;
	unsigned int i, j, Nn, Nel, d, Neq, dMin[NBTypes], dMax[NBTypes];
	double       *W, *nL, *XYZ;

	TestCase = malloc(STRLEN_MAX * sizeof *TestCase); // free
	strcpy(TestCase,"SupersonicVortex");
	DB.TestCase = TestCase;

	initialize_test_case_parameters();

	for (i = 0; i < NBTypes; i++)
		BType[i] = malloc(STRLEN_MIN * sizeof *BType[i]); // free

	strcpy(BType[0],"SlipWall    ");
	strcpy(BType[1],"Riemann     ");
	strcpy(BType[2],"BackPressure");

	dMin[0] = 1; dMax[0] = 3;
	dMin[1] = 2; dMax[1] = 3;
	dMin[2] = 2; dMax[2] = 2;

	for (i = 0; i < NBTypes; i++) {
	for (d = dMin[i]; d <= dMax[i]; d++) {
		if (strstr(BType[i],"Riemann")) {
			for (j = 0; j < 4; j++)
				TestDB.EnteredRiemann[j] = 0;
		} else if (strstr(BType[i],"BackPressure")) {
			for (j = 0; j < 2; j++)
				TestDB.EnteredBackPressure[j] = 0;
		}

		Neq = d+2;

		W    = initialize_W(&Nn,&Nel,d); // free
		nL   = initialize_n(Nn,Nel,d);   // free
		XYZ  = initialize_XYZ(Nn,Nel,d); // free
		if (strstr(BType[i],"BackPressure"))
			update_values(Nn,Nel,W,nL,d);
		pass = compare_jacobian_boundary(Nn,Nel,d,Neq,W,nL,XYZ,BType[i]);

		if (d == dMin[i]) {
			//     0         10        20        30        40        50
			printf("jacobian_boundary_%s (d = %d):          ",BType[i],d);
		} else {
			printf("                               (d = %d):          ",d);
		}
		test_print(pass);

		free(W);
		free(nL);
		free(XYZ);
	}}

	free(DB.TestCase);

	for (i = 0; i < NBTypes; i++)
		free(BType[i]);

	free(DB.SolverType); // From initialize_test_case_parameters
}
