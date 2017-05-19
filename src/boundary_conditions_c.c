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

#include "boundary_conditions.h"
#include "variable_functions_c.h"
#include "exact_solutions.h"

/*
 *	Purpose:
 *		Identical to boundary_conditions using complex variables (for complex step verification).
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void set_BC_from_BType(struct S_BC *const BCdata, char const *const BType)
{
	BCdata->ComputeQ = 0;
	if (strstr(BType,"SlipWall")) {
		BCdata->BC = BC_SLIPWALL;
	} else if (strstr(BType,"Riemann")) {
		BCdata->BC = BC_RIEMANN;
	} else if (strstr(BType,"BackPressure")) {
		BCdata->BC = BC_BACKPRESSURE;
	} else if (strstr(BType,"Total_TP")) {
		BCdata->BC = BC_TOTAL_TP;
	} else if (strstr(BType,"SupersonicIn")) {
		BCdata->BC = BC_SUPERSONIC_IN;
	} else if (strstr(BType,"SupersonicOut")) {
		BCdata->BC = BC_SUPERSONIC_OUT;
	} else if (strstr(BType,"NoSlip_Dirichlet")) {
		BCdata->BC = BC_NOSLIP_T;
		BCdata->ComputeQ = 1;
	} else if (strstr(BType,"NoSlip_Adiabatic")) {
		BCdata->BC = BC_NOSLIP_ADIABATIC;
		BCdata->ComputeQ = 1;
	} else if (strstr(BType,"Poisson_Dirichlet")) {
		BCdata->BC = BC_DIRICHLET;
		BCdata->ComputeQ = 1;
	} else if (strstr(BType,"Poisson_Neumann")) {
		BCdata->BC = BC_NEUMANN;
		BCdata->ComputeQ = 1;
	} else if (strstr(BType,"Advection_Inflow")) {
		BCdata->BC = BC_INFLOW;
	} else if (strstr(BType,"Advection_Outflow")) {
		BCdata->BC = BC_OUTFLOW;
	} else {
		EXIT_UNSUPPORTED;
	}
}

void correct_XYZ_for_exact_normal(struct S_BC *const BCdata, char const *const BType)
{
	/*
	 *	Purpose:
	 *		Ensure that the coordinates are on the annular boundary (necessary when 'ExactNormal' is enabled)
	 */

	if (!EXACT_NORMAL)
		return;

	if (strstr(BType,"SlipWall")) {
		unsigned int const Nn  = BCdata->Nn,
		                   Nel = BCdata->Nel;

		double *const XYZ = (double *const) BCdata->XYZ;

		unsigned int NnTotal = Nn*Nel;

		for (size_t n = 0; n < NnTotal; n++) {
			double const x = XYZ[0*NnTotal+n],
			             y = XYZ[1*NnTotal+n],
			             r = sqrt(x*x+y*y);

			XYZ[0*NnTotal+n] = x/r;
			XYZ[1*NnTotal+n] = y/r;
		}
	}
}

static void boundary_Advection_c         (struct S_BC *const BCdata);
static void boundary_Poisson_c           (struct S_BC *const BCdata);
static void boundary_Riemann_c           (struct S_BC *const BCdata);
static void boundary_SlipWall_c          (struct S_BC *const BCdata);
static void boundary_BackPressure_c      (struct S_BC *const BCdata);
static void boundary_Total_TP_c          (struct S_BC *const BCdata);
static void boundary_SupersonicInflow_c  (struct S_BC *const BCdata);
static void boundary_SupersonicOutflow_c (struct S_BC *const BCdata);
static void boundary_NoSlip_Dirichlet_c  (struct S_BC *const BCdata);
static void boundary_NoSlip_Adiabatic_c  (struct S_BC *const BCdata);

void compute_boundary_values_c(struct S_BC *const BCdata)
{
	switch((BCdata->BC) % BC_STEP_SC) {
		case BC_RIEMANN:          boundary_Riemann_c(BCdata);           break;
		case BC_SLIPWALL: {
			if (EXACT_SLIPWALL) {
				unsigned int const d    = BCdata->d,
				                   Nvar = d+2,
				                   Nn   = BCdata->Nn,
				                   Nel  = BCdata->Nel;

				unsigned int const NnTotal = Nn*Nel;

				double *WB = BCdata->WB;
				BCdata->WB = malloc(NnTotal*Nvar * sizeof *(BCdata->WB)); // free
				compute_exact_boundary_solution(BCdata);

				for (size_t i = 0; i < NnTotal*Nvar; i++)
					BCdata->WB_c[i] = BCdata->WB[i];

				free(BCdata->WB);
				BCdata->WB = WB;
			} else if (EXACT_NORMAL) {
				double const *const nL = BCdata->nL;

				BCdata->nL = compute_exact_boundary_normal(BCdata);
				boundary_SlipWall_c(BCdata);

				free((double *) BCdata->nL);
				BCdata->nL = nL;
			} else {
				boundary_SlipWall_c(BCdata);
			}
			break;
		}
		case BC_BACKPRESSURE:     boundary_BackPressure_c(BCdata);      break;
		case BC_TOTAL_TP:         boundary_Total_TP_c(BCdata);          break;
		case BC_SUPERSONIC_IN:    boundary_SupersonicInflow_c(BCdata);  break;
		case BC_SUPERSONIC_OUT:   boundary_SupersonicOutflow_c(BCdata); break;
		case BC_NOSLIP_T:         boundary_NoSlip_Dirichlet_c(BCdata);  break;
		case BC_NOSLIP_ADIABATIC: boundary_NoSlip_Adiabatic_c(BCdata);  break;
		case BC_DIRICHLET:        // fallthrough
		case BC_NEUMANN:          boundary_Poisson_c(BCdata);           break;
		case BC_INFLOW:           // fallthrough
		case BC_OUTFLOW:          boundary_Advection_c(BCdata);         break;
		default:
			printf("%d\n",BCdata->BC);
			EXIT_UNSUPPORTED;
			break;
	}
}

static void get_boundary_values_c(const double X, const double Y, double complex *const rho, double complex *const u,
                                  double complex *const v, double complex *const w, double complex *const p)
{
	// Initialize DB Parameters
	char   *TestCase = DB.TestCase;


	if (strstr(TestCase,"SupersonicVortex")) {
		// Use the exact solution for the Outer VOLUME
		double rIn   = DB.rIn,
	           MIn   = DB.MIn,
	           rhoIn = DB.rhoIn,
	           VIn   = DB.VIn;

		double         r, t;
		double complex Vt;

		r = sqrt(X*X+Y*Y);
		t = atan2(Y,X);

		*rho = rhoIn*pow(1.0+0.5*GM1*MIn*MIn*(1.0-pow(rIn/r,2.0)),1.0/GM1);
		*p   = cpow(*rho,GAMMA)/GAMMA;

		Vt = -VIn/r;
		*u = -sin(t)*Vt;
		*v =  cos(t)*Vt;
		*w =  0.0;
	} else if (strstr(TestCase,"InviscidChannel")) {
		// Use exact uniform channel solution for outer VOLUME
		*rho = DB.rhoInf;
		*p   = DB.pInf;
		*u   = DB.MInf*DB.cInf;
		*v   = 0.0;
		*w   = 0.0;
	} else if (strstr(TestCase,"SubsonicNozzle")) {
		if (fabs(Y) < EPS) { // Inflow
			*rho = DB.rhoInf;
			*p   = DB.pInf;
			*u   = 0.0;
			*v   = DB.MInf*DB.cInf;
			*w   = 0.0;
		} else if (fabs(X) < EPS) { // Outflow
			printf("Error: Use BackPressure BC here as the outlet state is not known.\n"), EXIT_MSG;
		} else {
			printf("Error: Unsupported.\n"), EXIT_MSG;
		}
	} else if (strstr(TestCase,"PrandtlMeyer")) {
		// Use supersonic inflow solution
		*rho = DB.rhoIn;
		*p   = DB.pIn;
		*u   = DB.VIn;
		*v   = 0.0;
		*w   = 0.0;
	} else {
		printf("TestCase: %s\n",TestCase);
		printf("Error: Unsupported TestCase.\n"), EXIT_MSG;
	}
}

static void boundary_Advection_c(struct S_BC *const BCdata)
{
	unsigned int const BC_index = (BCdata->BC) % BC_STEP_SC;
	if (!(BC_index == BC_INFLOW || BC_index == BC_OUTFLOW))
		EXIT_UNSUPPORTED;

	unsigned int const Nn      = BCdata->Nn,
	                   Nel     = BCdata->Nel,
	                   NnTotal = Nn*Nel;

	double *const XYZB = compute_XYZ_boundary(BCdata); // free

	double complex const *const WL = BCdata->WL_c;
	double complex       *const WB = BCdata->WB_c;

	if (BC_index == BC_INFLOW) {
		double *const WBr = malloc(NnTotal * sizeof *WBr); // free
		compute_exact_solution(NnTotal,XYZB,WBr,0);
		for (size_t n = 0; n < NnTotal; n++)
			WB[n] = WBr[n];
		free(WBr);
	} else if (BC_index == BC_OUTFLOW) {
		for (size_t n = 0; n < NnTotal; n++)
			WB[n] = WL[n];
	}
}

static void boundary_Poisson_c(struct S_BC *const BCdata)
{
	unsigned int const BC_index = (BCdata->BC) % BC_STEP_SC;
	if (!(BC_index == BC_DIRICHLET || BC_index == BC_NEUMANN))
		EXIT_UNSUPPORTED;

	unsigned int const d       = BCdata->d,
	                   Nn      = BCdata->Nn,
	                   Nel     = BCdata->Nel,
	                   NnTotal = Nn*Nel;

	if (d <= 1)
		EXIT_UNSUPPORTED;

	double *const XYZB = compute_XYZ_boundary(BCdata); // free

	double complex const *const WL = BCdata->WL_c,
	                     *const *const QL = BCdata->QL_c;
	double complex       *const WB = BCdata->WB_c,
	                     *const *const QB = BCdata->QB_c;

	if (BC_index == BC_DIRICHLET) {
		double *const WBr = malloc(NnTotal * sizeof *WBr); // free
		compute_exact_solution(NnTotal,XYZB,WBr,0);
		for (size_t n = 0; n < NnTotal; n++) {
			WB[n]  = 2.0*WBr[n];
			WB[n] -= WL[n];
		}
		free(WBr);

		if (!(QL == NULL)) {
			for (size_t dim = 0; dim < d; dim++) {
				for (size_t n = 0; n < NnTotal; n++)
					QB[dim][n] = QL[dim][n];
			}
		}
	} else if (BC_index == BC_NEUMANN) {
		for (size_t n = 0; n < NnTotal; n++)
			WB[n] = WL[n];

		if (!(QL == NULL)) {
			double *const QB_1D = malloc(NnTotal*d * sizeof *QB_1D); // free
			compute_exact_gradient(NnTotal,XYZB,QB_1D);
			for (size_t dim = 0; dim < d; dim++) {
				for (size_t n = 0; n < NnTotal; n++)
					QB[dim][n] = -QL[dim][n] + 2.0*QB_1D[NnTotal*dim+n];
			}
			free(QB_1D);
		}
	}

	free(XYZB);
}
static void boundary_Riemann_c(struct S_BC *const BCdata)
{
	unsigned int const d    = BCdata->d,
	                   Nn   = BCdata->Nn,
	                   Nel  = BCdata->Nel;

	double const *const XYZ = BCdata->XYZ,
	             *const nL  = BCdata->nL;

	double complex const *const WL = BCdata->WL_c;
	double complex       *const WB = BCdata->WB_c;

	// Initialize DB Parameters
	char         *TestCase = DB.TestCase;
	double       rIn       = DB.rIn,
	             MIn       = DB.MIn,
	             rhoIn     = DB.rhoIn,
	             VIn       = DB.VIn;

	// Standard datatypes
	unsigned int   i, j, Indn, NnTotal;
	double         *n, r, t;
	double complex *rhoL, *uL, *vL, *wL, *pL, cL, *VnL, sL, *rhoR, *uR, *vR, *wR, *pR, cR, *VnR, sR, *UL, *UR,
	               *rhoB, *uB, *vB, *wB, *pB, *UB,
	               Vt, RL, RR, Vn, c, ut, vt, wt;
	const double   *X, *Y;

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

static void boundary_SlipWall_c(struct S_BC *const BCdata)
{
	unsigned int const d    = BCdata->d,
	                   Nn   = BCdata->Nn,
	                   Nel  = BCdata->Nel;

	double const *const nL  = BCdata->nL;

	double complex const *const WL = BCdata->WL_c;
	double complex       *const WB = BCdata->WB_c;

	unsigned int         i, NnTotal, IndE;
	double complex       *rhoB, *rhouB, *rhovB, *rhowB, *EB, rhoVL;
	const double complex *rhoL, *rhouL, *rhovL, *rhowL, *EL;

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

static void boundary_BackPressure_c(struct S_BC *const BCdata)
{
	unsigned int const d    = BCdata->d,
	                   Nvar = d+2,
	                   Nn   = BCdata->Nn,
	                   Nel  = BCdata->Nel;

//	double const *const nL  = BCdata->nL;

	double complex const *const WL = BCdata->WL_c;
	double complex       *const WB = BCdata->WB_c;

	// Standard datatypes
	unsigned int   n, NnTotal, eq, var, IndW;
	double complex rhoL, rhoL_inv, uL, vL, wL, EL, VL, V2L, pL, rhoB, cL, c2L, *WB_ptr[Nvar];
	const double complex *rhoL_ptr, *rhouL_ptr, *rhovL_ptr, *rhowL_ptr, *EL_ptr, *WL_ptr[Nvar];

	// silence
	rhovL_ptr = rhowL_ptr = NULL;

	NnTotal = Nn*Nel;

	for (eq = 0; eq < Nvar; eq++) {
		WL_ptr[eq] = &WL[eq*NnTotal];
		WB_ptr[eq] = &WB[eq*NnTotal];
	}

	double complex zeros[NnTotal];

	for (n = 0; n < NnTotal; n++)
		zeros[n] = 0.0;

	rhoL_ptr  = WL_ptr[0];
	rhouL_ptr = WL_ptr[1];
	EL_ptr    = WL_ptr[d+1];

//	double const n_ptr = nL;

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
		VL  = csqrt(V2L);

		pL  = GM1*(EL-0.5*rhoL*V2L);

		c2L = GAMMA*pL/rhoL;
		cL  = csqrt(c2L);

//		n1 = *n_ptr++;
//		if      (d == 3) { n2 = *n_ptr++; n3 = *n_ptr++; }
//		else if (d == 2) { n2 = *n_ptr++; n3 = 0.0;      }
//		else if (d == 1) { n2 = 0.0;      n3 = 0.0;      }

//		double complex const VnL = uL*n1+vL*n2+wL*n3;
//		if (creal(VnL) < 0.0) // Inlet
//			printf("Warning: Velocity Inflow in boundary_BackPressure_c.\n");

		if (cabs(VL) >= cabs(cL)) { // Supersonic
			TestDB.EnteredBackPressure[0]++;
			for (var = 0; var < Nvar; var++) {
				*WB_ptr[IndW] = *WL_ptr[IndW];
				IndW++;
			}
		} else {
			TestDB.EnteredBackPressure[1]++;
			double pBack = DB.pBack;

			if (pBack < EPS)
				printf("Error: Initialize pBack.\n"), EXIT_MSG;

			rhoB = GAMMA*pBack/c2L;

			*WB_ptr[IndW++] = rhoB;
			*WB_ptr[IndW++] = uL*rhoB;
			if (d == 3) {
				*WB_ptr[IndW++] = vL*rhoB;
				*WB_ptr[IndW++] = wL*rhoB;
			} else if (d == 2) {
				*WB_ptr[IndW++] = vL*rhoB;
			}
			// Note: Using VL for the boundary
			*WB_ptr[IndW++] = pBack/GM1+0.5*rhoB*V2L;
		}

		for (var = 0; var < Nvar; var++) {
			WL_ptr[var]++;
			WB_ptr[var]++;
		}
	}
}

static void boundary_Total_TP_c(struct S_BC *const BCdata)
{
	unsigned int const d    = BCdata->d,
	                   Nvar = d+2,
	                   Nn   = BCdata->Nn,
	                   Nel  = BCdata->Nel;

	double const *const nL  = BCdata->nL;

	double complex const *const WL = BCdata->WL_c;
	double complex       *const WB = BCdata->WB_c;

	// Initialize DB Parameters
	double Rg      = DB.Rg,
	       p_Total = DB.p_Total,
	       T_Total = DB.T_Total;

	// Standard datatypes
	unsigned int   NnTotal;
	const double complex *rhoL_ptr, *rhouL_ptr, *rhovL_ptr, *rhowL_ptr = NULL, *EL_ptr, *WL_ptr[Nvar];
	double complex *WB_ptr[Nvar];
	const double   *n_ptr;

	NnTotal = Nn*Nel;

	for (size_t var = 0; var < Nvar; var++) {
		WL_ptr[var] = &WL[var*NnTotal];
		WB_ptr[var] = &WB[var*NnTotal];
	}

	double complex zeros[NnTotal];

	for (size_t n = 0; n < NnTotal; n++)
		zeros[n] = 0.0;

	rhoL_ptr  = WL_ptr[0];
	rhouL_ptr = WL_ptr[1];
	rhovL_ptr = WL_ptr[2];
	EL_ptr    = WL_ptr[d+1];

	if (d == 3) {
		rhowL_ptr = WL_ptr[3];
	} else if (d == 2) {
		rhowL_ptr = zeros;
	}

	n_ptr = nL;

	for (size_t n = 0; n < NnTotal; n++) {
		unsigned int   IndW = 0;
		double complex rhoL, rhoL_inv, uL, vL, wL, EL, V2L, pL, cL, HL, n1, n2, n3, VnL, RL;

		// silence
		n3 = 0.0;

		rhoL = *rhoL_ptr++;
		rhoL_inv = 1.0/rhoL;

		uL   = (*rhouL_ptr++)*rhoL_inv;
		vL   = (*rhovL_ptr++)*rhoL_inv;
		wL   = (*rhowL_ptr++)*rhoL_inv;
		EL   = *EL_ptr++;

		V2L = uL*uL+vL*vL+wL*wL;

		pL  = GM1*(EL-0.5*rhoL*V2L);
		cL  = csqrt(GAMMA*pL/rhoL);

		HL = (EL+pL)*rhoL_inv;

		n1 = *n_ptr++;
		n2 = *n_ptr++;
		if (d == 3)
			n3 = *n_ptr++;

		VnL = uL*n1+vL*n2+wL*n3;

		RL = VnL + 2.0/GM1*cL;

		// Solve for c
		double complex aQ, bQ, cQ, term1, term2, cM, cP, c, Vn, M, T, p, rho, u, v, w, E;

		aQ =  1.0 + 2.0/GM1;
		bQ = -2.0*RL;
		cQ =  0.5*GM1*(RL*RL - 2.0*HL);

		term1 = -bQ/(2.0*aQ);
		term2 = csqrt(bQ*bQ-4.0*aQ*cQ)/(2.0*aQ);

		cM = term1-term2;
		cP = term1+term2;

		// c = max(cM,cP)
		if (cabs(cM) > cabs(cP))
			c = cM;
		else
			c = cP;

		Vn = RL - 2.0/GM1*c;

		M = Vn/c;

		T = T_Total/(1+0.5*GM1*M*M);
		p = p_Total*cpow(T/T_Total,GAMMA/GM1);

		rho = p/(Rg*T);
		u   = Vn*n1;
		v   = Vn*n2;
		w   = Vn*n3;

		E   = p/GM1+0.5*rho*(u*u+v*v+w*w);

		*WB_ptr[IndW++] = rho;
		*WB_ptr[IndW++] = rho*u;
		*WB_ptr[IndW++] = rho*v;
		if (d == 3)
			*WB_ptr[IndW++] = rho*w;
		*WB_ptr[IndW++] = E;

		for (size_t var = 0; var < Nvar; var++)
			WB_ptr[var]++;
	}
}

static void boundary_SupersonicInflow_c(struct S_BC *const BCdata)
{
	unsigned int const d    = BCdata->d,
	                   Nvar = d+2,
	                   Nn   = BCdata->Nn,
	                   Nel  = BCdata->Nel;

	double const *const XYZ = BCdata->XYZ;

	double complex       *const WB = BCdata->WB_c;

	unsigned int   NnTotal;
	const double   *X_ptr, *Y_ptr;
	double complex *WB_ptr[Nvar];

	NnTotal = Nn*Nel;

	X_ptr = &XYZ[NnTotal*0];
	Y_ptr = &XYZ[NnTotal*1];

	for (size_t var = 0; var < Nvar; var++)
		WB_ptr[var] = &WB[(var)*NnTotal];

	for (size_t n = 0; n < NnTotal; n++) {
		unsigned int   IndW = 0;
		double         X, Y;
		double complex rhoR, uR, vR, wR, pR, V2R, ER;

		X = *X_ptr++;
		Y = *Y_ptr++;

		get_boundary_values_c(X,Y,&rhoR,&uR,&vR,&wR,&pR);

		V2R = uR*uR+vR*vR+wR*wR;
		ER = pR/GM1+0.5*rhoR*V2R;

		*WB_ptr[IndW++] = rhoR;
		*WB_ptr[IndW++] = rhoR*uR;
		*WB_ptr[IndW++] = rhoR*vR;

		if (d == 3)
			*WB_ptr[IndW++] = rhoR*wR;

		*WB_ptr[IndW++] = ER;

		for (size_t var = 0; var < Nvar; var++)
			WB_ptr[var]++;
	}
}

static void boundary_SupersonicOutflow_c(struct S_BC *const BCdata)
{
	unsigned int const d    = BCdata->d,
	                   Nvar = d+2,
	                   Nn   = BCdata->Nn,
	                   Nel  = BCdata->Nel;

	double complex const *const WL = BCdata->WL_c;
	double complex       *const WB = BCdata->WB_c;

	unsigned int const NnTotal = Nn*Nel;

	for (size_t i = 0, iMax = NnTotal*Nvar; i < iMax; i++)
		WB[i] = WL[i];
}

static void boundary_NoSlip_Dirichlet_c(struct S_BC *const BCdata)
{
	unsigned int const d       = BCdata->d,
	                   Nvar    = d+2,
	                   Nn      = BCdata->Nn,
	                   Nel     = BCdata->Nel,
	                   NnTotal = Nn*Nel;

	double const *const XYZ = BCdata->XYZ;

	double complex const *const        WL = BCdata->WL_c,
	                     *const *const QL = BCdata->QL_c;
	double complex       *const        WB = BCdata->WB_c,
	                     *const *const QB = BCdata->QB_c;

	double const *X_ptr = &XYZ[NnTotal*0],
	             *Y_ptr = &XYZ[NnTotal*1];

	double complex const *WL_ptr[Nvar];
	double complex       *WB_ptr[Nvar];
	for (size_t var = 0; var < Nvar; var++) {
		WL_ptr[var] = &WL[var*NnTotal];
		WB_ptr[var] = &WB[var*NnTotal];
	}

	double complex zeros[NnTotal];
	for (size_t n = 0; n < NnTotal; n++)
		zeros[n] = 0.0;

	double complex const *rhoL_ptr  = WL_ptr[0],
	                     *rhouL_ptr = WL_ptr[1],
	                     *rhovL_ptr = WL_ptr[2],
	                     *rhowL_ptr            ,
	                     *EL_ptr    = WL_ptr[d+1];

	if (d == 3) {
		rhowL_ptr = WL_ptr[3];
	} else if (d == 2) {
		rhowL_ptr = zeros;
	} else {
		EXIT_UNSUPPORTED;
	}

	for (size_t n = 0; n < NnTotal; n++) {
		double complex const rhoL     = *rhoL_ptr++,
		                     rhoL_inv = 1.0/rhoL,

		                     uL   = (*rhouL_ptr++)*rhoL_inv,
		                     vL   = (*rhovL_ptr++)*rhoL_inv,
		                     wL   = (*rhowL_ptr++)*rhoL_inv,
		                     EL   = *EL_ptr++,

		                     V2L = uL*uL+vL*vL+wL*wL,
		                     pL  = GM1*(EL-0.5*rhoL*V2L);

		bool   ApplyExtraBC = 0;
		double uB = 0.0, vB = 0.0, wB = 0.0, TB = 0.0;
		if (strstr(DB.TestCase,"TaylorCouette")) {
			double const X  = X_ptr[n],
			             Y  = Y_ptr[n],
			             t  = atan2(Y,X),
			             Vt = DB.omega*DB.rIn;
			uB = -sin(t)*Vt;
			vB =  cos(t)*Vt;
			wB =  0.0;
			TB =  DB.TIn;

			ApplyExtraBC = 1;
		} else {
			EXIT_UNSUPPORTED;
		}

		if (!ApplyExtraBC) {
			EXIT_UNSUPPORTED;
			printf("%f %f\n",creal(pL),TB); // Update if used (Entropy variables) (ToBeModified)
		} else {
			if (!strstr(DB.TestCase,"TaylorCouette") &&
			    !strstr(DB.TestCase,"PlaneCouette")) {
				EXIT_UNSUPPORTED;
			}

			size_t IndW = 0;
			double const rhoB = DB.rhoIn,
			             pB   = DB.pIn;
			*WB_ptr[IndW++] = -rhoL    + 2.0*rhoB;
			*WB_ptr[IndW++] = -rhoL*uL + 2.0*rhoB*uB;
			*WB_ptr[IndW++] = -rhoL*vL + 2.0*rhoB*vB;
			if (d == 3)
				*WB_ptr[IndW++] = -rhoL*wL + 2.0*rhoB*wB;
			*WB_ptr[IndW++] = -EL + 2.0*(pB/GM1+0.5*rhoB*(uB*uB+vB*vB+wB*wB));
		}

		for (size_t var = 0; var < Nvar; var++)
			WB_ptr[var]++;
	}

	// Set QB == QL (if necessary)
	if (QL == NULL)
		return;

	for (size_t dim = 0; dim < d; dim++) {
		for (size_t n = 0; n < NnTotal*Nvar; n++) {
			QB[dim][n] = QL[dim][n];
		}
	}
}

static void boundary_NoSlip_Adiabatic_c(struct S_BC *const BCdata)
{
	unsigned int const d    = BCdata->d,
	                   Nvar = d+2,
	                   Nn   = BCdata->Nn,
	                   Nel  = BCdata->Nel;

	double complex const *const        WL = BCdata->WL_c,
	                     *const *const QL = BCdata->QL_c;

	double complex       *const        WB = BCdata->WB_c,
	                     *const *const QB = BCdata->QB_c;

	unsigned int const NnTotal = Nn*Nel;

	double complex const *WL_ptr[Nvar];
	double complex       *WB_ptr[Nvar];

	for (size_t var = 0; var < Nvar; var++) {
		WL_ptr[var] = &WL[var*NnTotal];
		WB_ptr[var] = &WB[var*NnTotal];
	}

	double complex zeros[NnTotal];
	for (size_t n = 0; n < NnTotal; n++)
		zeros[n] = 0.0;

	double complex const *rhoL_ptr  = WL_ptr[0],
	                     *rhouL_ptr = WL_ptr[1],
	                     *rhovL_ptr = WL_ptr[2],
	                     *rhowL_ptr            ,
	                     *EL_ptr    = WL_ptr[d+1];

	if (d == 3) {
		rhowL_ptr = WL_ptr[3];
	} else if (d == 2) {
		rhowL_ptr = zeros;
	} else {
		EXIT_UNSUPPORTED;
	}

	for (size_t n = 0; n < NnTotal; n++) {
		double complex const rhoL     = *rhoL_ptr++,
		                     rhoL_inv = 1.0/rhoL,

		                     uL   = (*rhouL_ptr++)*rhoL_inv,
		                     vL   = (*rhovL_ptr++)*rhoL_inv,
		                     wL   = (*rhowL_ptr++)*rhoL_inv,
		                     EL   = *EL_ptr++;

		double complex u = 0.0, v = 0.0, w = 0.0;
		if (strstr(DB.TestCase,"PlaneCouette") ||
		    strstr(DB.TestCase,"TaylorCouette")) {
			; // Do nothing
		} else {
			EXIT_UNSUPPORTED;
		}

		size_t IndW = 0;
		*WB_ptr[IndW++] = rhoL;
		*WB_ptr[IndW++] = -rhoL*uL + 2.0*rhoL*u;
		*WB_ptr[IndW++] = -rhoL*vL + 2.0*rhoL*v;
		if (d == 3)
			*WB_ptr[IndW++] = -rhoL*wL + 2.0*rhoL*w;
		*WB_ptr[IndW++] = EL;

		for (size_t var = 0; var < Nvar; var++)
			WB_ptr[var]++;
	}

	if (QL == NULL)
		return;

	// Set QB == QL but set the Energy equation component of the numerical flux to 0. Alternatively, set QB here such
	// that the computed numerical flux is zero (More expensive and more difficult, but more general). (ToBeModified)
	for (size_t dim = 0; dim < d; dim++) {
		for (size_t n = 0; n < NnTotal*Nvar; n++)
			QB[dim][n] = QL[dim][n];
	}
}
