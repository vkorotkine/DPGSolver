// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "exact_solutions.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"

/*
 *	Purpose:
 *		Provide functions related to exact solutions.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void compute_exact_solution(const unsigned int Nn, double *XYZ, double *UEx, const unsigned int solved)
{
	// Initialize DB Parameters
	char *TestCase = DB.TestCase;

	// Standard datatypes
	unsigned int i;
	double       *X, *Y, *rhoEx, *uEx, *vEx, *wEx, *pEx;

	rhoEx = &UEx[Nn*0];
	uEx   = &UEx[Nn*1];
	vEx   = &UEx[Nn*2];
	wEx   = &UEx[Nn*3];
	pEx   = &UEx[Nn*4];

	X = &XYZ[0*Nn];
	Y = &XYZ[1*Nn];

	// Perhaps modify TestCase for Test_L2_proj and Test_update_h to make this cleaner, also in other functions (ToBeDeleted)
	if (strstr(TestCase,"PeriodicVortex") ||
	    strstr(TestCase,"Test_L2_proj")   ||
	    strstr(TestCase,"Test_update_h")) {
		// Initialize DB Parameters
		double       pInf           = DB.pInf,
					 TInf           = DB.TInf,
					 VInf           = DB.VInf,
					 uInf           = DB.uInf,
					 vInf           = DB.vInf,
					 wInf           = DB.wInf,
					 Rg             = DB.Rg,
					 Cscale         = DB.Cscale,
					 PeriodL        = DB.PeriodL,
					 PeriodFraction = DB.PeriodFraction,
					 Xc0            = DB.Xc,
					 Yc             = DB.Yc,
					 Rc             = DB.Rc;

		// Standard datatypes
		double DistTraveled, Xc, rhoInf, r2, C;

		rhoInf = pInf/(Rg*TInf);
		C      = Cscale*VInf;

		if (solved) {
			DistTraveled = PeriodL*PeriodFraction;
			Xc = Xc0 + DistTraveled;
			while (Xc > 0.5*PeriodL)
				Xc -= PeriodL;
		} else {
			Xc = Xc0;
		}

		for (i = 0; i < Nn; i++) {
			r2 = (pow(X[i]-Xc,2.0)+pow(Y[i]-Yc,2.0))/(Rc*Rc);

			uEx[i]   = uInf - C*(Y[i]-Yc)/(Rc*Rc)*exp(-0.5*r2);
			vEx[i]   = vInf + C*(X[i]-Xc)/(Rc*Rc)*exp(-0.5*r2);
			wEx[i]   = wInf;
			pEx[i]   = pInf - rhoInf*(C*C)/(2*Rc*Rc)*exp(-r2);
			rhoEx[i] = rhoInf;
		}
	} else if (strstr(TestCase,"SupersonicVortex") ||
	           strstr(TestCase,"Test_linearization")) {
		// Initialize DB Parameters
		double rIn   = DB.rIn,
		       MIn   = DB.MIn,
		       rhoIn = DB.rhoIn,
		       VIn   = DB.VIn;

		// Standard datatypes
		double r, t, Vt;

		for (i = 0; i < Nn; i++) {
			r = sqrt(X[i]*X[i]+Y[i]*Y[i]);
			t = atan2(Y[i],X[i]);

			rhoEx[i] = rhoIn*pow(1.0+0.5*GM1*MIn*MIn*(1.0-pow(rIn/r,2.0)),1.0/GM1);
			pEx[i]   = pow(rhoEx[i],GAMMA)/GAMMA;

			Vt = -VIn/r;
			uEx[i] = -sin(t)*Vt;
			vEx[i] =  cos(t)*Vt;
			wEx[i] =  0.0;
		}
	} else if (strstr(TestCase,"Poisson")) {
		for (i = 0; i < Nn; i++) {
			UEx[i] = sin(PI*X[i])*sin(PI*Y[i]);
		}
	} else {
		printf("Error: Unsupported TestCase.\n"), EXIT_MSG;
	}
}

void compute_exact_gradient(const unsigned int Nn, double *XYZ, double *QEx)
{
	// Initialize DB Parameters
	char         *TestCase = DB.TestCase;
	unsigned int d         = DB.d;

	// Standard datatypes
	unsigned int i;
	double       *X, *Y;

	if (d != 2)
		printf("Error: Unsupported d.\n"), EXIT_MSG;

	X = &XYZ[0*Nn];
	Y = &XYZ[1*Nn];

	if (strstr(TestCase,"Poisson")) {
		for (i = 0; i < Nn; i++) {
			QEx[Nn*0+i] = PI*cos(PI*X[i])*sin(PI*Y[i]);
			QEx[Nn*1+i] = PI*sin(PI*X[i])*cos(PI*Y[i]);
		}
	} else {
		printf("Error: Unsupported TestCase.\n"), EXIT_MSG;
	}
}

void compute_source(const unsigned int Nn, double *XYZ, double *source)
{
	/*
	 *	Purpose:
	 *		Computes source terms.
	 */

	// Initialize DB Parameters
	char         *TestCase = DB.TestCase;
	unsigned int Neq       = DB.Neq;

	// Standard datatypes
	unsigned int n, eq;
	double       *X, *Y;

	if (strstr(TestCase,"Poisson")) {
		X = &XYZ[Nn*0];
		Y = &XYZ[Nn*1];

		for (eq = 0; eq < Neq; eq++) {
			for (n = 0; n < Nn; n++)
				source[eq*Nn+n] = -2*PI*PI*sin(X[n])*sin(Y[n]);
		}
	} else {
		printf("Error: Unsupported TestCase.\n"), EXIT_MSG;
	}
}
