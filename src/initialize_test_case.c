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
		OPS->I_vG_vS = ELEMENT_OPS->I_vGs_vS[P][P][0];
	} else {
		OPS->I_vG_vS = ELEMENT_OPS->I_vGc_vS[P][P][0];
	}
}

static void initialize_PeriodicVortex(void)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;
	unsigned int Nvar = DB.Nvar;

	// Standard datatypes
	char         *SolverType;
	unsigned int n, NvnS;
	double       Xc, Yc, Rc, PeriodL, PeriodFraction, MInf, pInf, TInf, Rg, Cscale, uInf, vInf, wInf, VInf, rhoInf,
	             *XYZ_vS, *X_vS, *Y_vS, *r2, C,
	             *rho, *u, *v, *w, *p, *U, *What, *W;

	struct S_OPERATORS *OPS;
	struct S_VOLUME *VOLUME;

	SolverType = malloc(STRLEN_MIN * sizeof *SolverType); // keep
	strcpy(SolverType,"Explicit");

	Xc = -0.1;
//	Xc =  0.0;
	Yc =  0.0;
	Rc =  0.2;
	PeriodL        = 2.0;
	PeriodFraction = 0.1;
//	PeriodFraction = 1.0;

	MInf = 0.5;
//	MInf = 0.0;
	pInf = 1.0;
	TInf = 1.0;
	Rg   = 1.0;

	Cscale = 0.1;

	uInf   = MInf*sqrt(GAMMA*Rg*TInf);
	vInf   = 0.0;
	wInf   = 0.0;
	VInf   = sqrt(uInf*uInf+vInf*vInf+wInf*wInf);
	rhoInf = pInf/(Rg*TInf);

	C = Cscale*VInf;

	for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
		OPS = malloc(sizeof *OPS); // free
		init_ops(OPS,VOLUME);

		NvnS         = OPS->NvnS;
		VOLUME->NvnS = NvnS;

		What         = malloc(NvnS*Nvar * sizeof *What); // keep
		VOLUME->What = What;
		VOLUME->RES  = calloc(NvnS*Nvar , sizeof *(VOLUME->RES)); // keep

		XYZ_vS = malloc(NvnS*d    * sizeof *XYZ_vS); // free

		mm_CTN_d(NvnS,d,VOLUME->NvnG,OPS->I_vG_vS,VOLUME->XYZ,XYZ_vS);

		X_vS = &XYZ_vS[0*NvnS];
		Y_vS = &XYZ_vS[1*NvnS];

		r2 = malloc(NvnS * sizeof *r2); // free
		for (n = 0; n < NvnS; n++)
			r2[n] = (pow(X_vS[n]-Xc,2.0)+pow(Y_vS[n]-Yc,2.0))/(Rc*Rc);

        U   = malloc(NvnS*NVAR3D * sizeof *U); // free
		W   = malloc(NvnS*Nvar   * sizeof *W); // free

		rho = &U[NvnS*0];
		u   = &U[NvnS*1];
		v   = &U[NvnS*2];
		w   = &U[NvnS*3];
		p   = &U[NvnS*4];

		for (n = 0; n < NvnS; n++) {
			u[n]   = uInf - C*(Y_vS[n]-Yc)/(Rc*Rc)*exp(-0.5*r2[n]);
			v[n]   = vInf + C*(X_vS[n]-Xc)/(Rc*Rc)*exp(-0.5*r2[n]);
			w[n]   = wInf;
			p[n]   = pInf - rhoInf*(C*C)/(2*Rc*Rc)*exp(-r2[n]);
			rho[n] = rhoInf;
		}

		convert_variables(U,W,3,d,NvnS,1,'p','c');

/*
array_print_d(VOLUME->NvnG,d,XYZ_vS,'C');
array_print_d(NvnS,Nvar,W,'C');
exit(1);
*/

		mm_CTN_d(NvnS,Nvar,NvnS,OPS->ChiInvS_vS,W,What);

		free(XYZ_vS);
		free(r2);
		free(U);
		free(W);

		free(OPS);
	}

	DB.Xc = Xc;
	DB.Yc = Yc;
	DB.Rc = Rc;

	DB.MInf   = MInf;
	DB.pInf   = pInf;
	DB.TInf   = TInf;
	DB.VInf   = VInf;
	DB.uInf   = uInf;
	DB.vInf   = vInf;
	DB.wInf   = wInf;
	DB.Rg     = Rg;
	DB.Cscale = Cscale;

	DB.PeriodL        = PeriodL;
	DB.PeriodFraction = PeriodFraction;
	if (fabs(VInf) < 10*EPS)
		DB.FinalTime = 0.3;
	else
		DB.FinalTime = PeriodFraction*PeriodL/VInf;

	DB.SolverType     = SolverType;
}

static void initialize_SupersonicVortex(void)
{
	// Initialize DB Parameters
	unsigned int d    = DB.d,
	             Nvar = DB.Nvar;

	// Standard datatypes
	char         *SolverType;
	unsigned int i, NvnS;
	double       rIn, MIn, rhoIn, pIn, cIn, VIn, uIn, vIn, wIn, UIn[5], WIn[5], Vt,
	             *What, *XYZ_vS, *X_vS, *Y_vS, *r, *t, *rho, *u, *v, *w, *p, *U, *W;


	struct S_OPERATORS *OPS;
	struct S_VOLUME *VOLUME;

	SolverType = malloc(STRLEN_MIN * sizeof *SolverType); // keep
	strcpy(SolverType,"Explicit");
//	strcpy(SolverType,"Implicit");

	rIn  = 1.0;
//	rOut = 1.384;

	MIn = 2.25;

	rhoIn = 1.0;
	pIn   = pow(rhoIn,GAMMA)/GAMMA;

	cIn = sqrt(GAMMA*pIn/rhoIn);
	VIn = cIn*MIn/rIn;
	uIn = VIn;
	vIn = 0.0;
	wIn = 0.0;

	UIn[0] = rhoIn; UIn[1] = uIn; UIn[2] = vIn; UIn[3] = wIn; UIn[4] = pIn;
	convert_variables(UIn,WIn,3,3,1,1,'p','c');

	OPS = malloc(sizeof *OPS); // free
	for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
		init_ops(OPS,VOLUME);

		NvnS         = OPS->NvnS;
		VOLUME->NvnS = NvnS;

		What         = malloc(NvnS*Nvar * sizeof *What); // keep
		VOLUME->What = What;
		VOLUME->RES  = calloc(NvnS*Nvar , sizeof *(VOLUME->RES)); // keep

		XYZ_vS = malloc(NvnS*d    * sizeof *XYZ_vS); // free

		mm_CTN_d(NvnS,d,VOLUME->NvnG,OPS->I_vG_vS,VOLUME->XYZ,XYZ_vS);

		X_vS = &XYZ_vS[0*NvnS];
		Y_vS = &XYZ_vS[1*NvnS];

		r = malloc(NvnS * sizeof *r); // free
		t = malloc(NvnS * sizeof *t); // free
		for (i = 0; i < NvnS; i++) {
			r[i] = sqrt(X_vS[i]*X_vS[i]+Y_vS[i]*Y_vS[i]);
			t[i] = atan2(Y_vS[i],X_vS[i]);
		}

        U   = malloc(NvnS*NVAR3D * sizeof *U); // free
		W   = malloc(NvnS*Nvar   * sizeof *W); // free

		rho = &U[NvnS*0];
		u   = &U[NvnS*1];
		v   = &U[NvnS*2];
		w   = &U[NvnS*3];
		p   = &U[NvnS*4];

		for (i = 0; i < NvnS; i++) {
			rho[i] = rhoIn*pow(1.0+0.5*GM1*MIn*MIn*(1.0-pow(rIn/r[i],2.0)),1.0/GM1);
			p[i]   = pow(rho[i],GAMMA)/GAMMA;

			Vt = -MIn*cIn/r[i];
			u[i] = -sin(t[i])*Vt;
			v[i] =  cos(t[i])*Vt;
			w[i] =  0.0;
		}
		convert_variables(U,W,3,d,NvnS,1,'p','c');
		mm_CTN_d(NvnS,Nvar,NvnS,OPS->ChiInvS_vS,W,What);

		free(XYZ_vS);
		free(r);
		free(t);
		free(U);
		free(W);
	}
	free(OPS);

	DB.rIn   = rIn;
	DB.MIn   = MIn;
	DB.rhoIn = rhoIn;
	DB.VIn   = VIn;

	DB.SolverType = SolverType;
	DB.FinalTime  = 1e10;
}

void initialize_test_case(void)
{
	// Initialize DB Parameters
	char         *TestCase = DB.TestCase;
	unsigned int Testing   = DB.Testing;

	DB.Nvar = DB.d+2;
	DB.Neq  = DB.d+2;

	DB.OutputInterval = 1e3;

	// Standard datatypes
	unsigned int DOF0 = 0;
	struct S_VOLUME *VOLUME;

	if (strstr(TestCase,"dSphericalBump") != NULL)
		; // initialize_dSphericalBump();
	else if (strstr(TestCase,"GaussianBump") != NULL)
		; // initialize_GaussianBump();
	else if (strstr(TestCase,"PeriodicVortex") != NULL)
		initialize_PeriodicVortex();
	else if (strstr(TestCase,"PolynomialBump") != NULL)
		; // initialize_PolynomialBump();
	else if (strstr(TestCase,"SupersonicVortex") != NULL)
		initialize_SupersonicVortex();

	for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next)
		DOF0 += VOLUME->NvnS;
	DB.DOF0 = DOF0;

	if (Testing) {
		// Output initial solution to paraview
		output_to_paraview("ZTest_Sol_Init");
	}
}
