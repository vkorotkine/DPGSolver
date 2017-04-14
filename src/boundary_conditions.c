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
 *		Carlson(2011)-Inflow/Outflow_Boundary_Conditions_with_Application_to_FUN3D (NASA/TM–2011-217181)
 *		Toro(2009)-Riemann_Solvers_and_Numerical_Methods_for_Fluid_Dynamics
 */

void get_boundary_values(const double X, const double Y, double *const rho, double *const u, double *const v,
                         double *const w, double *const p)
{
	// Initialize DB Parameters
	char   *TestCase = DB.TestCase;


	if (strstr(TestCase,"SupersonicVortex")) {
		// Use the exact solution for the Outer VOLUME
		double rIn   = DB.rIn,
	           MIn   = DB.MIn,
	           rhoIn = DB.rhoIn,
	           VIn   = DB.VIn;
		double r, t, Vt;

		r = sqrt(X*X+Y*Y);
		t = atan2(Y,X);

		*rho = rhoIn*pow(1.0+0.5*GM1*MIn*MIn*(1.0-pow(rIn/r,2.0)),1.0/GM1);
		*p   = pow(*rho,GAMMA)/GAMMA;

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

void boundary_Riemann(const unsigned int Nn, const unsigned int Nel, const double *const XYZ, const double *const WL,
                      double *const WOut, double *const WB, const double *const nL, const unsigned int d)
{
	/*
	 *	Comments:
	 *		WOut is not used for all test cases.
	 *
	 *	References:
	 *		Carlson(2011): 2.2 (Note typo in eq. (14))
	 */

	// Standard datatypes
	unsigned int i, j, Indn, NnTotal;
	double       *rhoL, *uL, *vL, *wL, *pL, cL, *VnL, sL, *rhoR, *uR, *vR, *wR, *pR, cR, *VnR, sR, *UL, *UR,
	             *rhoB, *uB, *vB, *wB, *pB, *UB,
	             *n, RL, RR, Vn, c, ut, vt, wt;
	const double *X, *Y;

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

	// Inner/Outer VOLUME
	convert_variables(WL,UL,d,DMAX,Nn,Nel,'c','p');
	for (i = 0; i < NnTotal; i++) {
		Indn = i*DMAX;
		VnL[i] = n[Indn  ]*uL[i]+n[Indn+1]*vL[i]+n[Indn+2]*wL[i];

		get_boundary_values(X[i],Y[i],&rhoR[i],&uR[i],&vR[i],&wR[i],&pR[i]);
		VnR[i] = n[Indn  ]*uR[i]+n[Indn+1]*vR[i]+n[Indn+2]*wR[i];
	}

	for (i = 0; i < NnTotal; i++) {
		Indn = i*DMAX;
		cL = sqrt(GAMMA*pL[i]/rhoL[i]);
		cR = sqrt(GAMMA*pR[i]/rhoR[i]);

		// Riemann invariants
		RL = VnL[i] + 2.0/GM1*cL; // Outgoing
		RR = VnR[i] - 2.0/GM1*cR; // Incoming

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

void boundary_SlipWall(const unsigned int Nn, const unsigned int Nel, const double *const WL, double *const WB,
                       const double *const nL, const unsigned int d)
{
	/*
	 *	Comments:
	 *		Specifying equality of the total energy is equivalent to specifying equality of pressure because the
	 *		magnitude of the velocity is constant.
	 */

	// Standard datatypes
	unsigned int i, NnTotal, IndE;
	double       *rhoB, *rhouB, *rhovB, *rhowB, *EB, rhoVL;
	const double *rhoL, *rhouL, *rhovL, *rhowL, *EL;

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

void boundary_BackPressure(const unsigned int Nn, const unsigned int Nel, const double *const WL, double *const WB,
                           const double *const nL, const unsigned int d, const unsigned int Nvar)
{
	/*
	 *	Purpose:
	 *		Impose back Pressure (outflow) boundary condition.
	 *
	 *	References:
	 *		Carlson(2011): 2.4
	 */

	// Standard datatypes
	unsigned int n, NnTotal, var, IndW;
	double       rhoL, rhoL_inv, uL, vL, wL, EL, VL, V2L, pL, pBack, rhoB, cL, c2L, VnL, n1, n2, n3, *WB_ptr[Nvar];
	const double *rhoL_ptr, *rhouL_ptr, *rhovL_ptr, *rhowL_ptr, *EL_ptr, *WL_ptr[Nvar], *n_ptr;


	// silence
	n2 = n3 = 0;
	rhovL_ptr = rhowL_ptr = NULL;

	NnTotal = Nn*Nel;

	for (var = 0; var < Nvar; var++) {
		WL_ptr[var] = &WL[var*NnTotal];
		WB_ptr[var] = &WB[var*NnTotal];
	}

	double zeros[NnTotal];

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

		if (VnL < 0.0) // Inlet
			printf("\nWarning: Velocity Inflow in boundary_BackPressure.\n");

		if (fabs(VL) >= cL) { // Supersonic
//		if (1 || fabs(VL) >= cL) { // Supersonic
			for (var = 0; var < Nvar; var++) {
				*WB_ptr[IndW] = *WL_ptr[IndW];
				IndW++;
			}
		} else {
			pBack = DB.pBack;

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

void boundary_Total_TP(const unsigned int Nn, const unsigned int Nel, const double *const XYZ, const double *const WL,
                       double *const WB, const double *const nL, const unsigned int d, const unsigned int Nvar)
{
	/*
	 *	Purpose:
	 *		Impose total (P)ressure/(T)emperature (inflow) boundary condition.
	 *
	 *	Comments:
	 *		eq. (38/47) in Carlson(2011) implies that the velocity should be normal to the boundary. As the direction of
	 *		the flow velocity cannot be known, this implies that this boundary condition is not physically correct...
	 *
	 *	References:
	 *		Carlson(2011): 2.7
	 *		Toro(2009): (3.9), (8.58)
	 */

	// Initialize DB Parameters
	double Rg      = DB.Rg,
	       p_Total = DB.p_Total,
	       T_Total = DB.T_Total;

	// Standard datatypes
	unsigned int NnTotal;
	double       *WB_ptr[Nvar];
	const double *rhoL_ptr, *rhouL_ptr, *rhovL_ptr, *rhowL_ptr, *EL_ptr, *WL_ptr[Nvar], *n_ptr;

	// silence
	rhowL_ptr = NULL;
	WB[0] = XYZ[0];

	NnTotal = Nn*Nel;

	for (size_t var = 0; var < Nvar; var++) {
		WL_ptr[var] = &WL[var*NnTotal];
		WB_ptr[var] = &WB[var*NnTotal];
	}

	double zeros[NnTotal];

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
		unsigned int IndW = 0;
		double       rhoL, rhoL_inv, uL, vL, wL, EL, V2L, pL, cL, HL, n1, n2, n3, VnL, RL;

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
		cL  = sqrt(GAMMA*pL/rhoL);

		HL = (EL+pL)*rhoL_inv;

		n1 = *n_ptr++;
		n2 = *n_ptr++;
		if (d == 3)
			n3 = *n_ptr++;

		VnL = uL*n1+vL*n2+wL*n3;

		RL = VnL + 2.0/GM1*cL;

		// Solve for c
		double aQ, bQ, cQ, term1, term2, cM, cP, c, Vn, M, T, p, rho, u, v, w, E;

		aQ =  1.0 + 2.0/GM1;
		bQ = -2.0*RL;
		cQ =  0.5*GM1*(RL*RL - 2.0*HL);

		term1 = -bQ/(2.0*aQ);
		term2 = sqrt(bQ*bQ-4.0*aQ*cQ)/(2.0*aQ);

		cM = term1-term2;
		cP = term1+term2;

		// c = max(cM,cP)
		if (cM > cP)
			c = cM;
		else
			c = cP;

		Vn = RL - 2.0/GM1*c;

//		if (Vn > EPS)
//			printf("\nWarning: Velocity Outflow in boundary_Total_TP.\n");

		M = Vn/c;

		T = T_Total/(1+0.5*GM1*M*M);
		p = p_Total*pow(T/T_Total,GAMMA/GM1);

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

void boundary_SupersonicInflow(const unsigned int Nn, const unsigned int Nel, const double *const XYZ,
                               const double *const WL, double *const WB, const double *const nL, const unsigned int d,
                               const unsigned int Nvar)
{
	unsigned int NnTotal;
	double       *WB_ptr[Nvar];
	const double *X_ptr, *Y_ptr;

	// silence
	WB[0] = WL[0];
	WB[0] = nL[0];

	NnTotal = Nn*Nel;

	X_ptr = &XYZ[NnTotal*0];
	Y_ptr = &XYZ[NnTotal*1];

	for (size_t var = 0; var < Nvar; var++)
		WB_ptr[var] = &WB[(var)*NnTotal];

	for (size_t n = 0; n < NnTotal; n++) {
		unsigned int IndW = 0;
		double       X, Y, rhoR, uR, vR, wR, pR, V2R, ER;

		X = *X_ptr++;
		Y = *Y_ptr++;

		get_boundary_values(X,Y,&rhoR,&uR,&vR,&wR,&pR);

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

void boundary_SupersonicOutflow(const unsigned int Nn, const unsigned int Nel, const double *const XYZ,
                                const double *const WL, double *const WB, const double *const nL, const unsigned int d,
                                const unsigned int Nvar)
{
	unsigned int NnTotal;

	// silence
	WB[0] = XYZ[0];
	WB[0] = nL[0];
	NnTotal = d;

	NnTotal = Nn*Nel;

	for (size_t i = 0, iMax = NnTotal*Nvar; i < iMax; i++)
		WB[i] = WL[i];
}

void boundary_NoSlip_Dirichlet(struct S_BC *const BCdata)
{
	/*
	 *	Comments:
	 *		Following remark 11 in Nordstrom(2005), only four conditions are imposed for the Dirichlet boundary (No slip
	 *		with prescribed velocity).
	 *		This imposes the velocity and temperature on the boundary face by converting to the entropy variables,
	 *		setting the last d+1 values using the boundary conditions, then converting back to the conservative
	 *		variables.
	 *		Parsani(2014) discuss the correct method to impose the temperature boundary condition for the scheme to be
	 *		entropy stable (See Theorem 3.2 and eq. (56)). Investigate (ToBeModified).
	 *		The entropy variables of Barth(1998_Thesis, p. 16) are used.
	 *		See Hughes(1986) eq. (48)-(53) for conversion between entropy and conservative variables. Note the
	 *		additional GM1 factor and that rho*i == P/GM1.
	 *
	 *	References:
	 *		Nordstrom(2005)-Well-Posed_Boundary_Conditions_for_the_Navier-Stokes_Equations
	 *		Parsani(2014)-Entropy_Stable_Wall_Boundary_Conditions_for_the_Compressible_Navier-Stokes_Equations
	 *		Barth(1998_Thesis)-Simplified_Numerical_Methods_for_Gasdynamic_Systems_on_Triangulated_Domains
	 *		Hughes(1986)_A_New_Finite_Element_Formulation_for_Computational_Fluid_Dynamics_I
	 */

	unsigned int const d    = BCdata->d,
	                   Nvar = d+2,
	                   Nn   = BCdata->Nn,
	                   Nel  = BCdata->Nel;

	double const *const        XYZ    = BCdata->XYZ,
	             *const        WL     = BCdata->WL,
	             *const *const GradWL = BCdata->GradWL;

	double       *const        WB     = BCdata->WB,
	             *const *const GradWB = BCdata->GradWB;


	unsigned int const NnTotal = Nn*Nel;

	double const *X_ptr = &XYZ[NnTotal*0],
	             *Y_ptr = &XYZ[NnTotal*1];

	double const *WL_ptr[Nvar];
	double       *WB_ptr[Nvar];
	for (size_t var = 0; var < Nvar; var++) {
		WL_ptr[var] = &WL[var*NnTotal];
		WB_ptr[var] = &WB[var*NnTotal];
	}

	double zeros[NnTotal];
	for (size_t n = 0; n < NnTotal; n++)
		zeros[n] = 0.0;

	double const *rhoL_ptr  = WL_ptr[0],
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
		double const rhoL     = *rhoL_ptr++,
		             rhoL_inv = 1.0/rhoL,

		             uL   = (*rhouL_ptr++)*rhoL_inv,
		             vL   = (*rhovL_ptr++)*rhoL_inv,
		             wL   = (*rhowL_ptr++)*rhoL_inv,
		             EL   = *EL_ptr++,

		             V2L = uL*uL+vL*vL+wL*wL,
		             pL  = GM1*(EL-0.5*rhoL*V2L),
					 TL  = pL/(rhoL*DB.Rg);

		double uB = 0.0, vB = 0.0, wB = 0.0, TB = 0.0;
		if (strstr(DB.TestCase,"TaylorCouette")) {
			double const X  = X_ptr[n],
			             Y  = Y_ptr[n],
			             t  = atan2(Y,X),
			             Vt = DB.omega*DB.rIn;
if (0) {
			uB = -sin(t)*Vt;
			vB =  cos(t)*Vt;
			wB =  0.0;
			TB =  DB.TIn;
} else {
			uB = -uL + 2.0*(-sin(t)*Vt);
			vB = -vL + 2.0*( cos(t)*Vt);
			wB = -wL + 2.0*( 0.0);
			TB = -TL + 2.0*( DB.TIn);
}
		} else {
			EXIT_UNSUPPORTED;
		}

		// Compute boundary entropy variables
		double const sL = log(pL/pow(rhoL,GAMMA)),
		             rho_over_p = 1.0/(DB.Rg*TB); // Using the ideal gas law

		double V[NVAR3D];
		unsigned int IndV = 0;
		V[IndV++] =  (GAMMA+1.0-sL)/GM1-EL/pL;
		V[IndV++] =  rho_over_p*uB;
		V[IndV++] =  rho_over_p*vB;
		V[IndV++] =  rho_over_p*wB;
		V[IndV++] = -rho_over_p;

		// Convert to conservative variables
		double const V2 = V[1]*V[1]+V[2]*V[2]+V[3]*V[3], // V2 == (rho/p)^2*(u^2+v^2+w^2)
		             sB = GAMMA+GM1*(-V[0]+0.5*V2/V[4]),
		             pB = GM1*(pow(GM1/pow(-GM1*V[4],GAMMA),1.0/GM1)*exp(-sB/GM1));

		unsigned int IndW = 0;
		*WB_ptr[IndW++] = -pB*V[4];
		*WB_ptr[IndW++] =  pB*V[1];
		*WB_ptr[IndW++] =  pB*V[2];
		if (d == 3)
			*WB_ptr[IndW++] = pB*V[3];
		*WB_ptr[IndW++] = pB*(1.0/GM1-0.5*V2/V[4]);


		IndW = 0;
		double const rho = DB.rhoIn, p = DB.pIn;
		*WB_ptr[IndW++] = rho;
		*WB_ptr[IndW++] = rho*uB;
		*WB_ptr[IndW++] = rho*vB;
		if (d == 3)
			*WB_ptr[IndW++] = rho*wB;
		*WB_ptr[IndW++] = p/GM1+0.5*rho*(uB*uB+vB*vB+wB*wB);


		for (size_t var = 0; var < Nvar; var++)
			WB_ptr[var]++;
	}

	// Set QB == QL (if necessary)
	if (GradWL == NULL)
		return;

	for (size_t dim = 0; dim < d; dim++) {
		for (size_t n = 0; n < NnTotal*Nvar; n++) {
			GradWB[dim][n] = GradWL[dim][n];
		}
	}
}

void boundary_NoSlip_Adiabatic(struct S_BC *const BCdata)
{
	/*
	 *	Comments:
	 *		This imposes three velocity BCs on the solution and the final boundary condition on the solution gradient.
	 *		Parsani(2014) discusses the correct method to impose this last boundary condition for the scheme to be
	 *		entropy stable (See Theorem 3.2 and eq. (55)). Investigate (ToBeModified).
	 *
	 *	References:
	 *		Parsani(2014)-Entropy_Stable_Wall_Boundary_Conditions_for_the_Compressible_Navier-Stokes_Equations
	 */

	unsigned int const d    = BCdata->d,
	                   Nvar = d+2,
	                   Nn   = BCdata->Nn,
	                   Nel  = BCdata->Nel;

	double const *const        WL     = BCdata->WL,
	             *const *const GradWL = BCdata->GradWL;

	double       *const        WB     = BCdata->WB,
	             *const *const GradWB = BCdata->GradWB;

	unsigned int const NnTotal = Nn*Nel;

	double const *WL_ptr[Nvar];
	double       *WB_ptr[Nvar];

	for (size_t var = 0; var < Nvar; var++) {
		WL_ptr[var] = &WL[var*NnTotal];
		WB_ptr[var] = &WB[var*NnTotal];
	}

	double zeros[NnTotal];
	for (size_t n = 0; n < NnTotal; n++)
		zeros[n] = 0.0;

	double const *rhoL_ptr  = WL_ptr[0],
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
		double const rhoL     = *rhoL_ptr++,
		             rhoL_inv = 1.0/rhoL,

		             uL   = (*rhouL_ptr++)*rhoL_inv,
		             vL   = (*rhovL_ptr++)*rhoL_inv,
		             wL   = (*rhowL_ptr++)*rhoL_inv,
		             EL   = *EL_ptr++,

		             V2L = uL*uL+vL*vL+wL*wL,
		             pL  = GM1*(EL-0.5*rhoL*V2L);

		double uB = 0.0, vB = 0.0, wB = 0.0;
		if (strstr(DB.TestCase,"TaylorCouette")) {
			uB = -uL;
			vB = -vL;
			wB = -wL;
			; // Do nothing.
		} else {
			EXIT_UNSUPPORTED;
		}

		unsigned int IndW = 0;
		*WB_ptr[IndW++] = rhoL;
		*WB_ptr[IndW++] = rhoL*uB;
		*WB_ptr[IndW++] = rhoL*vB;
		if (d == 3)
			*WB_ptr[IndW++] = rhoL*wB;
		*WB_ptr[IndW++] = pL/GM1+0.5*rhoL*(uB*uB+vB*vB+wB*wB);

		for (size_t var = 0; var < Nvar; var++)
			WB_ptr[var]++;
	}

	if (GradWL == NULL)
		return;

	if (strstr(DB.TestCase,"TaylorCouette")) {
		// Set QB == QL but set the Energy equation component of the numerical flux to 0. Alternatively, set QB here
		// such that the computed numerical flux is zero (More expensive and more difficult, but more general).
		// (ToBeModified)
		for (size_t dim = 0; dim < d; dim++) {
			for (size_t n = 0; n < NnTotal*Nvar; n++) {
				GradWB[dim][n] = GradWL[dim][n];
			}
		}
	} else {
		EXIT_UNSUPPORTED;
	}
}
