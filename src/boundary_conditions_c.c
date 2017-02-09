// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "boundary_conditions_c.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <complex.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "Test.h"

#include "variable_functions_c.h"

/*
 *	Purpose:
 *		Identical to fluxes_inviscid using complex variables (for complex step verification).
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void boundary_Riemann_c(const unsigned int Nn, const unsigned int Nel, double *XYZ, double complex *WL,
                        double complex *WOut, double complex *WB, double *nL, const unsigned int d)
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
	unsigned int   i, j, Indn, NnTotal;
	double         *X, *Y, *n, r, t;
	double complex *rhoL, *uL, *vL, *wL, *pL, cL, *VnL, sL, *rhoR, *uR, *vR, *wR, *pR, cR, *VnR, sR, *UL, *UR,
	               *rhoB, *uB, *vB, *wB, *pB, *UB,
	               Vt, RL, RR, Vn, c, ut, vt, wt;

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
	convert_variables_c(WL,UL,d,DMAX,Nn,Nel,'c','p');
	for (i = 0; i < NnTotal; i++) {
		Indn = i*DMAX;
		VnL[i] = n[Indn  ]*uL[i]+n[Indn+1]*vL[i]+n[Indn+2]*wL[i];
	}

	// Outer VOLUME
	if (strstr(TestCase,"SupersonicVortex") ||
	    strstr(TestCase,"Test_linearization")) {
		// Use the exact solution for the Outer VOLUME
		for (i = 0; i < NnTotal; i++) {
			r = sqrt(X[i]*X[i]+Y[i]*Y[i]);
			t = atan2(Y[i],X[i]);

			rhoR[i] = rhoIn*pow(1.0+0.5*GM1*MIn*MIn*(1.0-pow(rIn/r,2.0)),1.0/GM1);
			pR[i]   = cpow(rhoR[i],GAMMA)/GAMMA;

			Vt = -VIn/r;
			uR[i] = -sin(t)*Vt;
			vR[i] =  cos(t)*Vt;
			wR[i] =  0.0;

			Indn = i*DMAX;
//			VnR[i] = n[Indn  ]*uR[i]+n[Indn+1]*vR[i]+n[Indn+2]*wR[i];
			VnR[i] = n[Indn  ]*uR[i]+n[Indn+1]*vR[i]; // wR == 0
		}
	} else {
		printf("Error: Unsupported TestCase.\n"), EXIT_MSG;
	}

	for (i = 0; i < NnTotal; i++) {
		Indn = i*DMAX;
		cL = csqrt(GAMMA*pL[i]/rhoL[i]);
		cR = csqrt(GAMMA*pR[i]/rhoR[i]);

		// Riemann invariants
		RL = VnL[i] + 2.0/GM1*cL;
		RR = VnR[i] - 2.0/GM1*cR;

		Vn = 0.5*(RL+RR);
		c  = 0.25*GM1*(RL-RR);

		if (cabs(Vn) >= cabs(c)) { // Supersonic
			if (creal(Vn) < 0.0) { // Inlet
				TestDB.EnteredRiemann[0]++;

				rhoB[i] = rhoR[i];
				uB[i]   = uR[i];
				vB[i]   = vR[i];
				wB[i]   = wR[i];
				pB[i]   = pR[i];
			} else {         // Outlet
				TestDB.EnteredRiemann[1]++;

				rhoB[i] = rhoL[i];
				uB[i]   = uL[i];
				vB[i]   = vL[i];
				wB[i]   = wL[i];
				pB[i]   = pL[i];
			}
		} else {                   // Subsonic
			if (creal(Vn) < 0.0) { // Inlet
				TestDB.EnteredRiemann[2]++;

				sR = csqrt(pR[i]/cpow(rhoR[i],GAMMA));

				ut = uR[i] - VnR[i]*n[Indn  ];
				vt = vR[i] - VnR[i]*n[Indn+1];
				wt = wR[i] - VnR[i]*n[Indn+2];

				rhoB[i] = cpow(1.0/GAMMA*c*c/(sR*sR),1.0/GM1);
			} else {         // Outlet
				TestDB.EnteredRiemann[3]++;

				sL = csqrt(pL[i]/cpow(rhoL[i],GAMMA));

				ut = uL[i] - VnL[i]*n[Indn  ];
				vt = vL[i] - VnL[i]*n[Indn+1];
				wt = wL[i] - VnL[i]*n[Indn+2];

				rhoB[i] = cpow(1.0/GAMMA*c*c/(sL*sL),1.0/GM1);
			}
			uB[i] = Vn*n[Indn  ] + ut;
			vB[i] = Vn*n[Indn+1] + vt;
			wB[i] = Vn*n[Indn+2] + wt;

			pB[i] = 1.0/GAMMA*c*c*rhoB[i];
		}
	}
	convert_variables_c(UB,WB,3,d,Nn,Nel,'p','c');

	free(UL);
	free(UR);
	free(UB);
	free(VnL);
	free(VnR);
	free(n);
}

void boundary_SlipWall_c(const unsigned int Nn, const unsigned int Nel, double complex *WL, double complex *WB,
                         double *nL, const unsigned int d)
{
	// Standard datatypes
	unsigned int   i, NnTotal, IndE;
	double complex *rhoL, *rhouL, *rhovL, *rhowL, *EL, *rhoB, *rhouB, *rhovB, *rhowB, *EB, rhoVL;

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

void boundary_BackPressure_c(const unsigned int Nn, const unsigned int Nel, double complex *WL, double complex *WB, double *nL,
                             const unsigned int d, const unsigned int Neq)
{
	// Standard datatypes
	unsigned int   n, NnTotal, eq, var, Nvar, IndW;
	double complex *rhoL_ptr, *rhouL_ptr, *rhovL_ptr, *rhowL_ptr, *EL_ptr,
	               rhoL, rhoL_inv, uL, vL, wL, EL, VL, V2L, pL, pInf, rhoB, cL, c2L, VnL,
	               *WL_ptr[Neq], *WB_ptr[Neq];
	double         n1, n2, n3, *n_ptr;

	NnTotal = Nn*Nel;
	Nvar    = Neq;

	for (eq = 0; eq < Neq; eq++) {
		WL_ptr[eq] = &WL[eq*NnTotal];
		WB_ptr[eq] = &WB[eq*NnTotal];
	}

	double complex zeros[NnTotal];

	for (n = 0; n < NnTotal; n++)
		zeros[n] = 0.0;

	rhoL_ptr  = WL_ptr[0];
	rhouL_ptr = WL_ptr[1];
	EL_ptr    = WL_ptr[d+1];

	n_ptr = nL;

	if (d == 3) {
		rhovL_ptr = WL_ptr[2];
		rhowL_ptr = WL_ptr[3];
	} else if (d == 2) {
		rhovL_ptr = WL_ptr[2];
		rhowL_ptr = zeros;
	} else if (d == 1) {
		rhovL_ptr = zeros;
		rhowL_ptr = zeros;
	}

	for (n = 0; n < NnTotal; n++) {
		IndW = 0;

		// Inner VOLUME
		rhoL     = *rhoL_ptr++;
		rhoL_inv = 1.0/rhoL;

		uL   = (*rhouL_ptr++)*rhoL_inv;
		vL   = (*rhovL_ptr++)*rhoL_inv;
		wL   = (*rhowL_ptr++)*rhoL_inv;
		EL   = *EL_ptr++;

		V2L = uL*uL+vL*vL+wL*wL;
		VL  = sqrt(V2L);

		pL  = GM1*(EL-0.5*rhoL*V2L);

		n1 = *n_ptr++;
		if      (d == 3) { n2 = *n_ptr++; n3 = *n_ptr++; }
		else if (d == 2) { n2 = *n_ptr++; n3 = 0.0;      }
		else if (d == 1) { n2 = 0.0;      n3 = 0.0;      }

		VnL = uL*n1+vL*n2+wL*n3;

		c2L = GAMMA*pL/rhoL;
		cL  = sqrt(c2L);

		if (creal(VnL) < 0.0) // Inlet
			printf("Warning: Velocity Inflow in boundary_BackPressure_c.\n");

		if (cabs(VL) >= cabs(cL)) { // Supersonic
			TestDB.EnteredPT[0]++;
			for (var = 0; var < Nvar; var++) {
				*WB_ptr[IndW] = *WL_ptr[IndW];
				IndW++;
			}
		} else {
			TestDB.EnteredPT[1]++;
			pInf = DB.pInf;

			rhoB = GAMMA*pInf/c2L;

			*WB_ptr[IndW++] = rhoB;
			*WB_ptr[IndW++] = uL*rhoB;
			if (d == 3) {
				*WB_ptr[IndW++] = vL*rhoB;
				*WB_ptr[IndW++] = wL*rhoB;
			} else if (d == 2) {
				*WB_ptr[IndW++] = vL*rhoB;
			}
			// Note: Using VL for the boundary
			*WB_ptr[IndW++] = pInf/GM1+0.5*rhoB*V2L;
		}

		for (var = 0; var < Nvar; var++) {
			WL_ptr[var]++;
			WB_ptr[var]++;
		}
	}
}
