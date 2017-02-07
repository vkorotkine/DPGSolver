// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "boundary_conditions.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"

#include "variable_functions.h"

/*
 *	Purpose:
 *		Compute boundary conditions for given input parameters: XYZ, WL, WOut, nL
 *
 *	Comments:
 *		If this is found to be slow after profiling, separate cases for different dimensions to avoid unnecessary
 *		operations (ToBeDeleted).
 *
 *	Notation:
 *
 *	References:
 */

void boundary_Riemann(const unsigned int Nn, const unsigned int Nel, double *XYZ, double *WL, double *WOut, double *WB,
                      double *nL, const unsigned int d)
{
	/*
	 *	Comments:
	 *		WOut is not used for all test cases.
	 */

	// Initialize DB Parameters
	char         *TestCase = DB.TestCase;
	double       rIn       = DB.rIn,
	             MIn       = DB.MIn,
	             rhoIn     = DB.rhoIn,
	             VIn       = DB.VIn;

	// Standard datatypes
	unsigned int i, j, Indn, NnTotal;
	double       *rhoL, *uL, *vL, *wL, *pL, cL, *VnL, sL, *rhoR, *uR, *vR, *wR, *pR, cR, *VnR, sR, *UL, *UR,
	             *rhoB, *uB, *vB, *wB, *pB, *UB,
	             *X, *Y, *n, r, t, Vt, RL, RR, Vn, c, ut, vt, wt;

	// silence
	UL = WOut;

	NnTotal = Nn*Nel;

	UL = malloc(NnTotal*NVAR3D * sizeof *UL); // free
	UR = malloc(NnTotal*NVAR3D * sizeof *UR); // free
	UB = malloc(NnTotal*NVAR3D * sizeof *UB); // free

	VnL = malloc(NnTotal*1 * sizeof *VnL); // free
	VnR = malloc(NnTotal*1 * sizeof *VnR); // free

	rhoL = &UL[NnTotal*0];
	uL   = &UL[NnTotal*1];
	vL   = &UL[NnTotal*2];
	wL   = &UL[NnTotal*3];
	pL   = &UL[NnTotal*4];

	rhoR = &UR[NnTotal*0];
	uR   = &UR[NnTotal*1];
	vR   = &UR[NnTotal*2];
	wR   = &UR[NnTotal*3];
	pR   = &UR[NnTotal*4];

	rhoB = &UB[NnTotal*0];
	uB   = &UB[NnTotal*1];
	vB   = &UB[NnTotal*2];
	wB   = &UB[NnTotal*3];
	pB   = &UB[NnTotal*4];

	X = &XYZ[NnTotal*0];
	Y = &XYZ[NnTotal*1];

	n = calloc(NnTotal*DMAX , sizeof *n); // free
	for (i = 0; i < NnTotal; i++) {
	for (j = 0; j < d; j++) {
		n[i*DMAX+j] = nL[i*d+j];
	}}

	// Inner VOLUME
	convert_variables(WL,UL,d,DMAX,Nn,Nel,'c','p');
	for (i = 0; i < NnTotal; i++) {
		Indn = i*DMAX;
		VnL[i] = n[Indn  ]*uL[i]+n[Indn+1]*vL[i]+n[Indn+2]*wL[i];
	}

	// Outer VOLUME
	if (strstr(TestCase,"SupersonicVortex")) {
		// Use the exact solution for the Outer VOLUME
		for (i = 0; i < NnTotal; i++) {
			r = sqrt(X[i]*X[i]+Y[i]*Y[i]);
			t = atan2(Y[i],X[i]);

			rhoR[i] = rhoIn*pow(1.0+0.5*GM1*MIn*MIn*(1.0-pow(rIn/r,2.0)),1.0/GM1);
			pR[i]   = pow(rhoR[i],GAMMA)/GAMMA;

			Vt = -VIn/r;
			uR[i] = -sin(t)*Vt;
			vR[i] =  cos(t)*Vt;
			wR[i] =  0.0;

			Indn = i*DMAX;
//			VnR[i] = n[Indn  ]*uR[i]+n[Indn+1]*vR[i]+n[Indn+2]*wR[i];
			VnR[i] = n[Indn  ]*uR[i]+n[Indn+1]*vR[i]; // wR == 0
		}
	} else if (strstr(TestCase,"InviscidChannel")) {
		// Use exact uniform channel solution for outer VOLUME
		for (i = 0; i < NnTotal; i++) {
			rhoR[i] = DB.rhoInf;
			pR[i]   = DB.pInf;
			uR[i]   = DB.MInf*DB.cInf;
			vR[i]   = 0.0;
			wR[i]   = 0.0;

			Indn = i*DMAX;
			VnR[i] = n[Indn]*uR[i]; // vR == wR == 0
		}
	} else {
		printf("TestCase: %s\n",TestCase);
		printf("Error: Unsupported TestCase.\n"), EXIT_MSG;
	}

	for (i = 0; i < NnTotal; i++) {
		Indn = i*DMAX;
		cL = sqrt(GAMMA*pL[i]/rhoL[i]);
		cR = sqrt(GAMMA*pR[i]/rhoR[i]);

		// Riemann invariants
		RL = VnL[i] + 2.0/GM1*cL;
		RR = VnR[i] - 2.0/GM1*cR;

		Vn = 0.5*(RL+RR);
		c  = 0.25*GM1*(RL-RR);

		if (fabs(Vn) >= c) { // Supersonic
			if (Vn < 0.0) { // Inlet
				rhoB[i] = rhoR[i];
				uB[i]   = uR[i];
				vB[i]   = vR[i];
				wB[i]   = wR[i];
				pB[i]   = pR[i];
			} else {         // Outlet
				rhoB[i] = rhoL[i];
				uB[i]   = uL[i];
				vB[i]   = vL[i];
				wB[i]   = wL[i];
				pB[i]   = pL[i];
			}
		} else {                   // Subsonic
			if (Vn < 0.0) { // Inlet
				sR = sqrt(pR[i]/pow(rhoR[i],GAMMA));

				ut = uR[i] - VnR[i]*n[Indn  ];
				vt = vR[i] - VnR[i]*n[Indn+1];
				wt = wR[i] - VnR[i]*n[Indn+2];

				rhoB[i] = pow(1.0/GAMMA*c*c/(sR*sR),1.0/GM1);
			} else {         // Outlet
				sL = sqrt(pL[i]/pow(rhoL[i],GAMMA));

				ut = uL[i] - VnL[i]*n[Indn  ];
				vt = vL[i] - VnL[i]*n[Indn+1];
				wt = wL[i] - VnL[i]*n[Indn+2];

				rhoB[i] = pow(1.0/GAMMA*c*c/(sL*sL),1.0/GM1);
			}
			uB[i] = Vn*n[Indn  ] + ut;
			vB[i] = Vn*n[Indn+1] + vt;
			wB[i] = Vn*n[Indn+2] + wt;

			pB[i] = 1.0/GAMMA*c*c*rhoB[i];
		}
	}
	convert_variables(UB,WB,3,d,Nn,Nel,'p','c');

	free(UL);
	free(UR);
	free(UB);
	free(VnL);
	free(VnR);
	free(n);
}

void boundary_SlipWall(const unsigned int Nn, const unsigned int Nel, double *WL, double *WB, double *nL,
                       const unsigned int d)
{
	// Standard datatypes
	unsigned int i, NnTotal, IndE;
	double *rhoL, *rhouL, *rhovL, *rhowL, *EL, *rhoB, *rhouB, *rhovB, *rhowB, *EB, rhoVL;

	NnTotal = Nn*Nel;
	IndE = d+1;

	rhoL  = &WL[NnTotal*0];
	rhouL = &WL[NnTotal*1];
	EL    = &WL[NnTotal*IndE];

	rhoB  = &WB[NnTotal*0];
	rhouB = &WB[NnTotal*1];
	EB    = &WB[NnTotal*IndE];

	for (i = 0; i < NnTotal; i++) {
		rhoB[i] = rhoL[i];
		EB[i]   = EL[i];
	}

	if (d == 3) {
		rhovL = &WL[NnTotal*2];
		rhowL = &WL[NnTotal*3];

		rhovB = &WB[NnTotal*2];
		rhowB = &WB[NnTotal*3];

		for (i = 0; i < NnTotal; i++) {
			rhoVL = nL[i*d  ]*rhouL[i]+nL[i*d+1]*rhovL[i]+nL[i*d+2]*rhowL[i];

			rhouB[i] = rhouL[i]-2.0*rhoVL*nL[i*d  ];
			rhovB[i] = rhovL[i]-2.0*rhoVL*nL[i*d+1];
			rhowB[i] = rhowL[i]-2.0*rhoVL*nL[i*d+2];
		}
	} else if (d == 2) {
		rhovL = &WL[NnTotal*2];

		rhovB = &WB[NnTotal*2];

		for (i = 0; i < NnTotal; i++) {
			rhoVL = nL[i*d  ]*rhouL[i]+nL[i*d+1]*rhovL[i];

			rhouB[i] = rhouL[i]-2.0*rhoVL*nL[i*d  ];
			rhovB[i] = rhovL[i]-2.0*rhoVL*nL[i*d+1];
		}
	} else if (d == 1) {
		for (i = 0; i < NnTotal; i++) {
			rhoVL = nL[i*d  ]*rhouL[i];

			rhouB[i] = rhouL[i]-2.0*rhoVL*nL[i*d  ];
		}
	}
}
