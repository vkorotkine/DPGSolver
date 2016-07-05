// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Initialize solution on all VOLUMEs for each test case.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

struct S_OPERATORS {
	unsigned int NvnS;
	double       *I_vG_vS, *ChiInvS_vS;
};

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME)
{
	// Standard datatypes
	unsigned int P, type, curved;
	struct S_ELEMENT *ELEMENT_OPS;

	P = VOLUME->P;
	type   = VOLUME->type;
	curved = VOLUME->curved;

	ELEMENT_OPS = get_ELEMENT_type(type);

	OPS->NvnS       = ELEMENT_OPS->NvnS[P];
	OPS->ChiInvS_vS = ELEMENT_OPS->ChiInvS_vS[P][P][0];
	if (!curved) {
		OPS->I_vG_vS = ELEMENT_OPS->I_vGs_vS[1][P][0];
	} else {
		OPS->I_vG_vS = ELEMENT_OPS->I_vGc_vS[P][P][0];
	}
}

void initialize_test_case(void)
{
	// Initialize DB Parameters
	char         *TestCase = DB.TestCase;
	unsigned int d         = DB.d,
	             Testing   = DB.Testing;

	DB.Nvar = d+2;
	DB.Neq  = d+2;

	unsigned int Nvar = DB.Nvar;

	DB.OutputInterval = 1e3;

	// Standard datatypes
	char         *SolverType;
	unsigned int DOF0 = 0;
	unsigned int NvnS;
	double       *XYZ_vS, *U, *s, *W, *What;

	struct S_OPERATORS *OPS;
	struct S_VOLUME    *VOLUME;

	// silence
	SolverType = NULL;

	if (strstr(TestCase,"dSphericalBump") != NULL) {
		; // initialize_dSphericalBump();
	} else if (strstr(TestCase,"GaussianBump") != NULL) {
		; // initialize_GaussianBump();
	} else if (strstr(TestCase,"PeriodicVortex") != NULL) {
		SolverType = malloc(STRLEN_MIN * sizeof *SolverType); // keep
		strcpy(SolverType,"Explicit");

//		DB.Xc = -DB.PeriodL*0.05;
		DB.Xc =  0.0;
		DB.Yc =  0.0;
		DB.Rc =  0.2;
//		DB.PeriodFraction = 0.1;
		DB.PeriodFraction = 1.0;

		DB.MInf = 0.5;
//		DB.MInf = 0.0;
		DB.pInf = 1.0;
		DB.TInf = 1.0;
		DB.Rg   = 1.0;

		DB.Cscale = 0.1;

		DB.uInf   = DB.MInf*sqrt(GAMMA*DB.Rg*DB.TInf);
		DB.vInf   = 0.1*EPS;
		DB.wInf   = 0.1*EPS;
		DB.VInf   = sqrt(DB.uInf*DB.uInf+DB.vInf*DB.vInf+DB.wInf*DB.wInf);

		if (fabs(DB.VInf) < 10*EPS)
			DB.FinalTime = 0.3;
		else
			DB.FinalTime = DB.PeriodFraction*DB.PeriodL/DB.VInf;
	} else if (strstr(TestCase,"PolynomialBump") != NULL) {
		; // initialize_PolynomialBump();
	} else if (strstr(TestCase,"SupersonicVortex") != NULL) {
		// Standard datatypes
		double       pIn, cIn;

		SolverType = malloc(STRLEN_MIN * sizeof *SolverType); // keep
		strcpy(SolverType,"Explicit");
//		strcpy(SolverType,"Implicit");

		DB.rIn  = 1.0;
//		rOut = 1.384;

		DB.MIn = 2.25;

		DB.rhoIn = 1.0;
		pIn   = pow(DB.rhoIn,GAMMA)/GAMMA;

		cIn = sqrt(GAMMA*pIn/DB.rhoIn);
		DB.VIn = cIn*DB.MIn/DB.rIn;

		DB.FinalTime  = 1e10;
	}
	DB.SolverType = SolverType;

	if (strstr(TestCase,"PeriodicVortex")   != NULL ||
	    strstr(TestCase,"SupersonicVortex") != NULL) {

		OPS = malloc(sizeof *OPS); // free
		for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
			init_ops(OPS,VOLUME);
//printf("%d\n",VOLUME->indexg);

			NvnS         = OPS->NvnS;
			VOLUME->NvnS = NvnS;

			What         = malloc(NvnS*Nvar * sizeof *What); // keep
			VOLUME->What = What;
			VOLUME->RES  = calloc(NvnS*Nvar , sizeof *(VOLUME->RES)); // keep

			XYZ_vS = malloc(NvnS*d    * sizeof *XYZ_vS); // free

			mm_CTN_d(NvnS,d,VOLUME->NvnG,OPS->I_vG_vS,VOLUME->XYZ,XYZ_vS);

			U   = malloc(NvnS*NVAR3D * sizeof *U); // free
			W   = malloc(NvnS*Nvar   * sizeof *W); // free
			s   = malloc(NvnS*1      * sizeof *s); // free

			compute_exact_solution(NvnS,XYZ_vS,U,s,0);
			free(s);

			convert_variables(U,W,3,d,NvnS,1,'p','c');
			mm_CTN_d(NvnS,Nvar,NvnS,OPS->ChiInvS_vS,W,What);

			free(XYZ_vS);
			free(U);
			free(W);
		}
		free(OPS);
	}

	for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next)
		DOF0 += VOLUME->NvnS;
	DB.DOF0 = DOF0;

	if (Testing) {
		// Output initial solution to paraview
		output_to_paraview("ZTest_Sol_Init");
	}
}
