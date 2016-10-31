// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "jacobian_boundary_conditions.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"

#include "variable_functions.h"

#include "array_print.h"

/*
 *	Purpose:
 *		Compute the jacobian of boundary conditions for given input parameters: XYZ, WL, WOut, nL
 *
 *	Comments:
 *		The jacobian_boundary_* functions compute the transpose of the standard dWdW matrices. (i.e. the terms are
 *		ordered in memory as: dWB0/dWL0 dWB1/dWL0 ... dWB(Neq)/dWL0 dWB0/dWL1 ...)
 *		If this is found to be slow after profiling, separate cases for different dimensions to avoid unnecessary
 *		operations (ToBeDeleted).
 *
 *	Notation:
 *
 *	References:
 */

void jacobian_boundary_Riemann(const unsigned int Nn, const unsigned int Nel, double *XYZ, double *WL, double *WOut,
                               double *dWdW, double *nL, const unsigned int d, const unsigned int Neq)
{
	/*
	 *	Jacobian Matrices [var * eq]
	 *
	 *	Supersonic Inlet
	 *	dWBdWL = [  0  0  0  0  0
	 *	            0  0  0  0  0
	 *	            0  0  0  0  0
	 *	            0  0  0  0  0
	 *	            0  0  0  0  0 ]
	 *
	 *	Supersonic Outlet
	 *	dWBdWL = [  1  0  0  0  0
	 *	            0  1  0  0  0
	 *	            0  0  1  0  0
	 *	            0  0  0  1  0
	 *	            0  0  0  0  1 ]
	 *
	 *	Subsonic Inlet/Outlet
	 *		See below.
	 *
	 */

	// Initialize DB Parameters
	char         *TestCase = DB.TestCase;
	double       rIn       = DB.rIn,
	             MIn       = DB.MIn,
	             rhoIn     = DB.rhoIn,
	             VIn       = DB.VIn;

	// Standard datatypes
	unsigned int i, iMax, n, eq, var, NnTotal, Nvar, InddWdW;
	double       *rhoL_ptr, *rhouL_ptr, *rhovL_ptr, *rhowL_ptr, *EL_ptr, *n_ptr, *X_ptr, *Y_ptr,
	             rhoL, rhoL_inv, uL, vL, wL, EL, V2L, pL, rhoR, uR, vR, wR, pR,
	             cL, RL, VnL, cR, RR, VnR, c, Vn,
	             X, Y, r, t, Vt, n1, n2, n3, *dWdW_ptr[Neq*Neq];

	// silence
	rhoL_ptr = WOut;
	rhovL_ptr = rhowL_ptr = NULL;
	n2 = n3 = 0.0;

	NnTotal = Nn*Nel;
	Nvar    = Neq;

	double zeros[NnTotal];

	for (i = 0; i < NnTotal; i++)
		zeros[i] = 0.0;

	rhoL_ptr  = &WL[        0];
	rhouL_ptr = &WL[NnTotal*1];
	EL_ptr    = &WL[NnTotal*(d+1)];

	n_ptr = nL;

	if (d == 3) {
		rhovL_ptr = &WL[NnTotal*2];
		rhowL_ptr = &WL[NnTotal*3];
	} else if (d == 2) {
		rhovL_ptr = &WL[NnTotal*2];
		rhowL_ptr = zeros;
	} else if (d == 1) {
		rhovL_ptr = zeros;
		rhowL_ptr = zeros;
	}

	X_ptr = &XYZ[NnTotal*0];
	Y_ptr = &XYZ[NnTotal*1];

	for (eq  = 0; eq  < Neq;  eq++)  {
	for (var = 0; var < Nvar; var++) {
		dWdW_ptr[eq*Nvar+var] = &dWdW[(eq*Nvar+var)*NnTotal];
	}}

	for (n = 0; n < NnTotal; n++) {
		InddWdW = 0;

		// Inner VOLUME
		rhoL     = *rhoL_ptr++;
		rhoL_inv = 1.0/rhoL;

		uL   = (*rhouL_ptr++)*rhoL_inv;
		vL   = (*rhovL_ptr++)*rhoL_inv;
		wL   = (*rhowL_ptr++)*rhoL_inv;
		EL   = *EL_ptr++;

		V2L = uL*uL+vL*vL+wL*wL;

		pL  = GM1*(EL-0.5*rhoL*V2L);

		n1 = *n_ptr++;
		if      (d == 3) { n2 = *n_ptr++; n3 = *n_ptr++; }
		else if (d == 2) { n2 = *n_ptr++; n3 = 0.0;      }
		else if (d == 1) { n2 = 0.0;      n3 = 0.0;      }

		VnL = uL*n1+vL*n2+wL*n3;

		X = *X_ptr++;
		Y = *Y_ptr++;

		// Outer VOLUME
		if (strstr(TestCase,"SupersonicVortex") ||
		    strstr(TestCase,"Test_linearization")) {
			// Use the exact solution for the Outer VOLUME
			r = sqrt(X*X+Y*Y);
			t = atan2(Y,X);

			rhoR = rhoIn*pow(1.0+0.5*GM1*MIn*MIn*(1.0-pow(rIn/r,2.0)),1.0/GM1);
			pR   = pow(rhoR,GAMMA)/GAMMA;

			Vt = -VIn/r;
			uR = -sin(t)*Vt;
			vR =  cos(t)*Vt;
			wR = 0.0;

			VnR = n1*uR+n2*vR; // wR == 0
		} else {
			printf("TestCase: %s\n",TestCase);
			printf("Error: Unsupported TestCase.\n"), EXIT_MSG;
		}

		cL = sqrt(GAMMA*pL/rhoL);
		cR = sqrt(GAMMA*pR/rhoR);

		// Riemann invariants
		RL = VnL + 2.0/GM1*cL;
		RR = VnR - 2.0/GM1*cR;

		Vn = 0.5*(RL+RR);
		c  = 0.25*GM1*(RL-RR);

		if (fabs(Vn) >= c) { // Supersonic
			if (Vn < 0.0) { // Inlet
//printf("j: Sup Inlet\n");
				for (var = 0; var < Nvar; var++) {
				for (eq = 0; eq < Neq; eq++) {
					*dWdW_ptr[InddWdW++] = 0.0;
				}}
			} else { // Outlet
//printf("j: Sup Outlet\n");
				for (var = 0; var < Nvar; var++) {
				for (eq = 0; eq < Neq; eq++) {
					if (var != eq)
						*dWdW_ptr[InddWdW++] = 0.0;
					else
						*dWdW_ptr[InddWdW++] = 1.0;
				}}
			}
		} else { // Subsonic
			double dcLdW, rho, u, v, w, V2, ut, vt, wt, un, vn, wn, cnst1, drhodW, dudW, dvdW, dwdW, dpdW,
			       drhoLdW[Nvar], duLdW[Nvar], dvLdW[Nvar], dwLdW[Nvar], dpLdW[Nvar],
			       dVnLdW[Nvar], dRLdW[Nvar], dcdW[Nvar];

			// silence
			un = vn = wn = 0.0;

			if (d == 3) {
				drhoLdW[0] = 1.0;     drhoLdW[1] = 0.0; drhoLdW[2] = 0.0; drhoLdW[3] = 0.0; drhoLdW[4] = 0.0;
				dpLdW[0]   = 0.5*V2L; dpLdW[1]   = -uL; dpLdW[2]   = -vL; dpLdW[3]   = -wL; dpLdW[4]   = 1.0;

				duLdW[0] = -uL*rhoL_inv; duLdW[1] = rhoL_inv; duLdW[2] = 0.0;      duLdW[3] = 0.0;      duLdW[4] = 0.0;
				dvLdW[0] = -vL*rhoL_inv; dvLdW[1] = 0.0;      dvLdW[2] = rhoL_inv; dvLdW[3] = 0.0;      dvLdW[4] = 0.0;
				dwLdW[0] = -wL*rhoL_inv; dwLdW[1] = 0.0;      dwLdW[2] = 0.0;      dwLdW[3] = rhoL_inv; dwLdW[4] = 0.0;

				for (var = 0; var < Nvar; var++) {
					dpLdW[var] *= GM1;
					dVnLdW[var] = duLdW[var]*n1 + dvLdW[var]*n2 + dwLdW[var]*n3;
				}

				un = Vn*n1;
				vn = Vn*n2;
				wn = Vn*n3;
			} else if (d == 2) {
				drhoLdW[0] = 1.0;     drhoLdW[1] = 0.0; drhoLdW[2] = 0.0; drhoLdW[3] = 0.0;
				dpLdW[0]   = 0.5*V2L; dpLdW[1]   = -uL; dpLdW[2]   = -vL; dpLdW[3]   = 1.0;

				duLdW[0] = -uL*rhoL_inv; duLdW[1] = rhoL_inv; duLdW[2] = 0.0;      duLdW[3] = 0.0;
				dvLdW[0] = -vL*rhoL_inv; dvLdW[1] = 0.0;      dvLdW[2] = rhoL_inv; dvLdW[3] = 0.0;

				for (var = 0; var < Nvar; var++) {
					dpLdW[var] *= GM1;
					dVnLdW[var] = duLdW[var]*n1 + dvLdW[var]*n2;
				}

				un = Vn*n1;
				vn = Vn*n2;
			} else if (d == 1) {
				drhoLdW[0] = 1.0;     drhoLdW[1] = 0.0; drhoLdW[2] = 0.0;
				dpLdW[0]   = 0.5*V2L; dpLdW[1]   = -uL; dpLdW[2]   = 1.0;

				duLdW[0] = -uL*rhoL_inv; duLdW[1] = rhoL_inv; duLdW[2] = 0.0;

				for (var = 0; var < Nvar; var++) {
					dpLdW[var] *= GM1;
					dVnLdW[var] = duLdW[var]*n1;
				}

				un = Vn*n1;
			}

			for (var = 0; var < Nvar; var++) {
				dcLdW      = 0.5*GAMMA/(cL*rhoL*rhoL)*(dpLdW[var]*rhoL-pL*drhoLdW[var]);
				dRLdW[var] = dVnLdW[var] + 2.0/GM1*dcLdW;
				dcdW[var]  = 0.25*GM1*dRLdW[var];
			}

			if (Vn < 0.0) { // Inlet
//printf("j: Sub Inlet\n");
				double sR;

				sR  = sqrt(pR/pow(rhoR,GAMMA));
				if (d == 3) {
					for (var = 0; var < Nvar; var++) {
						drhodW = pow(GAMMA,-1.0/GM1)*2.0/(GM1*sR)*pow(c/sR,-GM3/GM1)*dcdW[var];

						rho = pow(c*c/(GAMMA*sR*sR),1.0/GM1);

						ut = uR - VnR*n1;
						vt = vR - VnR*n2;
						wt = wR - VnR*n3;

						u = un + ut;
						v = vn + vt;
						w = wn + wt;
						V2 = u*u+v*v+w*w;

						dudW = 0.5*dRLdW[var]*n1;
						dvdW = 0.5*dRLdW[var]*n2;
						dwdW = 0.5*dRLdW[var]*n3;
						dpdW = (2.0*c*dcdW[var]*rho+c*c*drhodW)/GAMMA;

						*dWdW_ptr[InddWdW++] = drhodW;
						*dWdW_ptr[InddWdW++] = drhodW*u + rho*dudW;
						*dWdW_ptr[InddWdW++] = drhodW*v + rho*dvdW;
						*dWdW_ptr[InddWdW++] = drhodW*w + rho*dwdW;
						*dWdW_ptr[InddWdW++] = dpdW/GM1 + 0.5*(drhodW*V2+2.0*rho*(u*dudW+v*dvdW+w*dwdW));
					}
				} else if (d == 2) {
					for (var = 0; var < Nvar; var++) {
						drhodW = pow(GAMMA,-1.0/GM1)*2.0/(GM1*sR)*pow(c/sR,-GM3/GM1)*dcdW[var];

						rho = pow(c*c/(GAMMA*sR*sR),1.0/GM1);

						ut = uR - VnR*n1;
						vt = vR - VnR*n2;

						u = un + ut;
						v = vn + vt;
						V2 = u*u+v*v;

						dudW = 0.5*dRLdW[var]*n1;
						dvdW = 0.5*dRLdW[var]*n2;
						dpdW = (2.0*c*dcdW[var]*rho+c*c*drhodW)/GAMMA;

						*dWdW_ptr[InddWdW++] = drhodW;
						*dWdW_ptr[InddWdW++] = drhodW*u + rho*dudW;
						*dWdW_ptr[InddWdW++] = drhodW*v + rho*dvdW;
						*dWdW_ptr[InddWdW++] = dpdW/GM1 + 0.5*(drhodW*V2+2.0*rho*(u*dudW+v*dvdW));
					}
				} else if (d == 1) {
					for (var = 0; var < Nvar; var++) {
						drhodW = pow(GAMMA,-1.0/GM1)*2.0/(GM1*sR)*pow(c/sR,-GM3/GM1)*dcdW[var];

						rho = pow(c*c/(GAMMA*sR*sR),1.0/GM1);

						ut = uR - VnR*n1;
						u  = un + ut;
						V2 = u*u;

						dudW = 0.5*dRLdW[var]*n1;
						dpdW = (2.0*c*dcdW[var]*rho+c*c*drhodW)/GAMMA;

						*dWdW_ptr[InddWdW++] = drhodW;
						*dWdW_ptr[InddWdW++] = drhodW*u + rho*dudW;
						*dWdW_ptr[InddWdW++] = dpdW/GM1 + 0.5*(drhodW*V2+2.0*rho*(u*dudW));
					}
				}
			} else { // Outlet
//printf("j: Sub Outlet\n");
				double sL, dsLdW;

				sL  = sqrt(pL/pow(rhoL,GAMMA));
				if (d == 3) {
					for (var = 0; var < Nvar; var++) {
						dsLdW = 0.5*sqrt(pow(rhoL,GAMMA)/pL)/pow(rhoL,2.0*GAMMA)*
								(dpLdW[var]*pow(rhoL,GAMMA)-GAMMA*pL*pow(rhoL,GM1)*drhoLdW[var]);

						drhodW = pow(GAMMA,-1.0/GM1)*2.0/GM1*pow(c,-GM3/GM1)*
								 pow(sL,-(GAMMA+1.0)/GM1)*(dcdW[var]*sL-c*dsLdW);

						rho = pow(c*c/(GAMMA*sL*sL),1.0/GM1);

						ut = uL - VnL*n1;
						vt = vL - VnL*n2;
						wt = wL - VnL*n3;

						u = un + ut;
						v = vn + vt;
						w = wn + wt;
						V2 = u*u+v*v+w*w;

						cnst1 = 0.5*dRLdW[var]-dVnLdW[var];
						dudW = duLdW[var]+n1*cnst1;
						dvdW = dvLdW[var]+n2*cnst1;
						dwdW = dwLdW[var]+n3*cnst1;
						dpdW = (2.0*c*dcdW[var]*rho+c*c*drhodW)/GAMMA;

						*dWdW_ptr[InddWdW++] = drhodW;
						*dWdW_ptr[InddWdW++] = drhodW*u + rho*dudW;
						*dWdW_ptr[InddWdW++] = drhodW*v + rho*dvdW;
						*dWdW_ptr[InddWdW++] = drhodW*w + rho*dwdW;
						*dWdW_ptr[InddWdW++] = dpdW/GM1 + 0.5*(drhodW*V2+2.0*rho*(u*dudW+v*dvdW+w*dwdW));
					}
				} else if (d == 2) {
					for (var = 0; var < Nvar; var++) {
						dsLdW = 0.5*sqrt(pow(rhoL,GAMMA)/pL)/pow(rhoL,2.0*GAMMA)*
								(dpLdW[var]*pow(rhoL,GAMMA)-GAMMA*pL*pow(rhoL,GM1)*drhoLdW[var]);

						drhodW = pow(GAMMA,-1.0/GM1)*2.0/GM1*pow(c,-GM3/GM1)*
								 pow(sL,-(GAMMA+1.0)/GM1)*(dcdW[var]*sL-c*dsLdW);

						rho = pow(c*c/(GAMMA*sL*sL),1.0/GM1);

						ut = uL - VnL*n1;
						vt = vL - VnL*n2;

						u = un + ut;
						v = vn + vt;
						V2 = u*u+v*v;

						cnst1 = 0.5*dRLdW[var]-dVnLdW[var];
						dudW = duLdW[var]+n1*cnst1;
						dvdW = dvLdW[var]+n2*cnst1;
						dpdW = (2.0*c*dcdW[var]*rho+c*c*drhodW)/GAMMA;

						*dWdW_ptr[InddWdW++] = drhodW;
						*dWdW_ptr[InddWdW++] = drhodW*u + rho*dudW;
						*dWdW_ptr[InddWdW++] = drhodW*v + rho*dvdW;
						*dWdW_ptr[InddWdW++] = dpdW/GM1 + 0.5*(drhodW*V2+2.0*rho*(u*dudW+v*dvdW));
					}
				} else if (d == 1) {
					for (var = 0; var < Nvar; var++) {
						dsLdW = 0.5*sqrt(pow(rhoL,GAMMA)/pL)/pow(rhoL,2.0*GAMMA)*
								(dpLdW[var]*pow(rhoL,GAMMA)-GAMMA*pL*pow(rhoL,GM1)*drhoLdW[var]);

						drhodW = pow(GAMMA,-1.0/GM1)*2.0/GM1*pow(c,-GM3/GM1)*
								 pow(sL,-(GAMMA+1.0)/GM1)*(dcdW[var]*sL-c*dsLdW);

						rho = pow(c*c/(GAMMA*sL*sL),1.0/GM1);

						ut = uL - VnL*n1;
						u = un + ut;
						V2 = u*u;

						cnst1 = 0.5*dRLdW[var]-dVnLdW[var];
						dudW = duLdW[var]+n1*cnst1;
						dpdW = (2.0*c*dcdW[var]*rho+c*c*drhodW)/GAMMA;

						*dWdW_ptr[InddWdW++] = drhodW;
						*dWdW_ptr[InddWdW++] = drhodW*u + rho*dudW;
						*dWdW_ptr[InddWdW++] = dpdW/GM1 + 0.5*(drhodW*V2+2.0*rho*(u*dudW));
					}
				}
			}
		}
		for (i = 0, iMax = Neq*Nvar; i < iMax; i++)
			dWdW_ptr[i]++;
	}
}

void jacobian_boundary_SlipWall(const unsigned int Nn, const unsigned int Nel, double *WL, double *dWdW, double *nL,
                                const unsigned int d, const unsigned int Neq)
{
	/*
	 *	Jacobian Matrix [var * eq]
	 *
	 *	dWBdWL = [  1  0          0         0         0
	 *	            0  1-2*n1*n1   -2*n2*n1  -2*n3*n1 0
	 *	            0   -2*n1*n2  1-2*n2*n2  -2*n3*n2 0
	 *	            0   -2*n1*n3   -2*n2*n3 1-2*n3*n3 0
	 *	            0  0          0         0         1 ]
	 */

	// Standard datatypes
	unsigned int i, iMax, n, eq, var, NnTotal, Nvar, InddWdW;
	double       n1, n2, n3, *dWdW_ptr[Neq*Neq], *n_ptr;

	// silence
	dWdW_ptr[0] = WL;

	NnTotal = Nn*Nel;
	Nvar    = Neq;

	for (eq  = 0; eq  < Neq;  eq++)  {
	for (var = 0; var < Nvar; var++) {
		dWdW_ptr[eq*Nvar+var] = &dWdW[(eq*Nvar+var)*NnTotal];
	}}
	n_ptr = nL;

	if (d == 3) {
		for (n = 0; n < NnTotal; n++) {
			n1 = *n_ptr++;
			n2 = *n_ptr++;
			n3 = *n_ptr++;

			InddWdW = 0;

			// *** var 1 ***
			*dWdW_ptr[InddWdW++] = 1.0;
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 0.0;

			// *** var 2 ***
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 1.0 - 2.0*n1*n1;
			*dWdW_ptr[InddWdW++] =     - 2.0*n2*n1;
			*dWdW_ptr[InddWdW++] =     - 2.0*n3*n1;
			*dWdW_ptr[InddWdW++] = 0.0;

			// *** var 3 ***
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] =     - 2.0*n1*n2;
			*dWdW_ptr[InddWdW++] = 1.0 - 2.0*n2*n2;
			*dWdW_ptr[InddWdW++] =     - 2.0*n3*n2;
			*dWdW_ptr[InddWdW++] = 0.0;

			// *** var 4 ***
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] =     - 2.0*n1*n3;
			*dWdW_ptr[InddWdW++] =     - 2.0*n2*n3;
			*dWdW_ptr[InddWdW++] = 1.0 - 2.0*n3*n3;
			*dWdW_ptr[InddWdW++] = 0.0;

			// *** var 5 ***
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 1.0;

			for (i = 0, iMax = Neq*Nvar; i < iMax; i++)
				dWdW_ptr[i]++;
		}
	} else if (d == 2) {
		for (n = 0; n < NnTotal; n++) {
			n1 = *n_ptr++;
			n2 = *n_ptr++;

			InddWdW = 0;

			// *** var 1 ***
			*dWdW_ptr[InddWdW++] = 1.0;
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 0.0;

			// *** var 2 ***
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 1.0 - 2.0*n1*n1;
			*dWdW_ptr[InddWdW++] =     - 2.0*n2*n1;
			*dWdW_ptr[InddWdW++] = 0.0;

			// *** var 3 ***
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] =     - 2.0*n1*n2;
			*dWdW_ptr[InddWdW++] = 1.0 - 2.0*n2*n2;
			*dWdW_ptr[InddWdW++] = 0.0;

			// *** var 4 ***
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 1.0;

			for (i = 0, iMax = Neq*Nvar; i < iMax; i++)
				dWdW_ptr[i]++;
		}
	} else if (d == 1) {
		for (n = 0; n < NnTotal; n++) {
			n1 = *n_ptr++;

			InddWdW = 0;

			// *** var 1 ***
			*dWdW_ptr[InddWdW++] = 1.0;
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 0.0;

			// *** var 2 ***
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 1.0 - 2.0*n1*n1;
			*dWdW_ptr[InddWdW++] = 0.0;

			// *** var 3 ***
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 0.0;
			*dWdW_ptr[InddWdW++] = 1.0;

			for (i = 0, iMax = Neq*Nvar; i < iMax; i++)
				dWdW_ptr[i]++;
		}
	}
}
