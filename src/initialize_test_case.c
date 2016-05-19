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
	double       *I_vG_vS;
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

	OPS->NvnS = ELEMENT_OPS->NvnS[P];
	if (!curved) {
		OPS->I_vG_vS = ELEMENT_OPS->I_vGs_vS[P];
	} else {
		OPS->I_vG_vS = ELEMENT_OPS->I_vGc_vS[P];
	}
}

static void initialize_PeriodicVortex(void)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	DB.Nvar = d+2;
	DB.Neq  = d+2;

	unsigned int Nvar = DB.Nvar,
	             Neq  = DB.Neq;

	// Standard datatypes
	unsigned int NvnS;
	double       Xc, Yc, Rc, PeriodL, MInf, pInf, TInf, Rg, Cscale, uInf, vInf, wInf, VInf,
	             *I_vG_vS, *What, *W;

	struct S_OPERATORS *OPS;
	struct S_VOLUME *VOLUME;

	Xc = -0.1;
	Yc =  0.0;
	Rc =  0.2;
	PeriodL = 2.0;

	MInf = 0.5;
	pInf = 1.0;
	TInf = 1.0;
	Rg   = 1.0;

	Cscale = 0.1;

	uInf = MInf*sqrt(GAMMA*Rg*TInf);
	vInf = 0.0;
	wInf = 0.0;
	VInf = sqrt(pow(uInf,2.0)+pow(vInf,2.0)+pow(wInf,2.0));

	for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {

		init_ops(OPS,VOLUME);

		NvnS    = OPS->NvnS;
		I_vG_vS = OPS->I_vG_vS;

		What   = malloc(NvnS*Nvar * sizeof *What);   // keep
		XYZ_vS = malloc(NvnS*d    * sizeof *XYZ_vS); // free

		mm_CTN_d(NvnS,d,VOLUME->NvnG,I_vG_vS,VOLUME->XYZ,XYZ_vS);

		X_vS = &XYZ_vS[0*NvnS];
		Y_vS = &XYZ_vS[1*NvnS];

		r2 = malloc(NvnS * sizeof *r2); // free
		for (n = 0; n < NvnS; n++)
			r2[n] = (pow(X_vS[n]-Xc,2.0)+pow(Y_vS[n]-Yc,2.0))/pow(Rc,2.0);

// set memory for u v w
		u[n] = uInf - C*(Y_vS[n]-Yc)/pow(Rc,2.0)*exp(-0.5*r2[n]);
		v[n] = vInf + C*(X_vS[n]-Xc)/pow(Rc,2.0)*exp(-0.5*r2[n]);
		w[n] = 0.0;

		free(XYZ_vS);
		free(r2);
	}

}

void initialize_test_case(void)
{
	// Initialize DB Parameters
	char *TestCase = DB.TestCase;

	if (strstr(TestCase,"dSphericalBump") != NULL)
		; // initialize_dSphericalBump();
	else if (strstr(TestCase,"GaussianBump") != NULL)
		; // initialize_GaussianBump();
	else if (strstr(TestCase,"PeriodicVortex") != NULL)
		initialize_PeriodicVortex();
	else if (strstr(TestCase,"PolynomialBump") != NULL)
		; // initialize_PolynomialBump();
	else if (strstr(TestCase,"SupersonicVortex") != NULL)
		; // initialize_SupersoncVortex();
}
