// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "jacobian_fluxes_inviscid.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Parameters.h"

/*
 *	Purpose:
 *		Compute inviscid flux jacobians from input W in conservative form.
 *
 *	Comments:
 *		It is assumed that inputs: W, n and outputs: dFdW (Add others: ToBeModified) are vectorized (i.e. the memory
 *		ordering is by equation and not by element).
 *
 *	Notation:
 *
 *	References:
 *		Toro(2009)-Riemann_Solvers_and_Numerical_Methods_for_Fluid_Dynamics (Ch. 3.2)
 *		   Note: There is a typo in the 3D Jacobian matrix (eq. 3.79), compare with eq. 3.70.
 */

void jacobian_flux_inviscid(const unsigned int Nn, const unsigned int Nel, double *W, double *dFdW,
                            const unsigned int d, const unsigned int Neq)
{
	/*
	 *	Jacobian Matrices [ eq * var]
	 *
	 *	dF1dW = [  0                1          0          0          0
	 *	          -u^2+GM1/2*V^2   -GM3*u     -GM1*v     -GM1*w      GM1
	 *	          -u*v              v          u          0          0
	 *	          -u*w              w          0          u          0
	 *	           u*(GM1/2*V^2-H)  H-GM1*u^2 -GM1*u*v   -GM1*u*w    G*u ]
	 *
	 *	dF2dW = [  0                0          1          0          0
	 *	          -u*v              v          u          0          0
	 *	          -v^2+GM1/2*V^2   -GM1*u     -GM3*v     -GM1*w      GM1
	 *	          -v*w              0          w          v          0
	 *	           v*(GM1/2*V^2-H) -GM1*u*v    H-GM1*v^2 -GM1*v*w    G*v ]
	 *
	 *	dF1dW = [  0                0          0          1          0
	 *	          -u*w              w          0          u          0
	 *	          -v*w              0          w          v          0
	 *	          -w^2+GM1/2*V^2   -GM1*u     -GM1*v     -GM3*w      GM1
	 *	           w*(GM1/2*V^2-H) -GM1*u*w   -GM1*v*w    H-GM1*w^2  G*w ]
	 *
	 */

	// Standard datatypes
	unsigned int i, n, eq, var, dim, iMax, Nvar, NnTotal, InddFdW;
	double *rho_ptr, *rhou_ptr, *rhov_ptr, *rhow_ptr, *E_ptr,
	       rho, u, v, w, u2, uv, uw, v2, vw, w2, V2, E, p, H, alpha, beta, *dFdW_ptr[DMAX*Neq*Neq];

	Nvar    = Neq;
	NnTotal = Nn*Nel;

	rho_ptr  = &W[NnTotal*0];
	rhou_ptr = &W[NnTotal*1];
	E_ptr    = &W[NnTotal*(d+1)];

	for (eq  = 0; eq  < Neq;  eq++)  {
	for (var = 0; var < Nvar; var++) {
	for (dim = 0; dim < d;    dim++) {
		dFdW_ptr[(eq*Nvar+var)*DMAX+dim] = &dFdW[((eq*Nvar+var)*d+dim)*NnTotal];
	}}}

	if (d == 3) {
		rhov_ptr = &W[NnTotal*2];
		rhow_ptr = &W[NnTotal*3];

		for (n = 0; n < NnTotal; n++) {
			rho = *rho_ptr;
			u   = (*rhou_ptr)/rho;
			v   = (*rhov_ptr)/rho;
			w   = (*rhow_ptr)/rho;
			E   = *E_ptr;

			u2 = u*u;
			uv = u*v;
			uw = u*w;
			v2 = v*v;
			vw = v*w;
			w2 = w*w;

			V2 = u2+v2+w2;
			p  = GM1*(E-0.5*rho*V2);
			H  = (E+p)/rho;

			alpha = 0.5*GM1*V2;
			beta  = alpha-H;

			InddFdW = 0;
			// *** eq 1 ***
			// var 1
			*dFdW_ptr[InddFdW++] =  0.0;
			*dFdW_ptr[InddFdW++] =  0.0;
			*dFdW_ptr[InddFdW++] =  0.0;

			// var 2
			*dFdW_ptr[InddFdW++] =  1.0;
			*dFdW_ptr[InddFdW++] =  0.0;
			*dFdW_ptr[InddFdW++] =  0.0;

			// var 3
			*dFdW_ptr[InddFdW++] =  0.0;
			*dFdW_ptr[InddFdW++] =  1.0;
			*dFdW_ptr[InddFdW++] =  0.0;

			// var 4
			*dFdW_ptr[InddFdW++] =  0.0;
			*dFdW_ptr[InddFdW++] =  0.0;
			*dFdW_ptr[InddFdW++] =  1.0;

			// var 5
			*dFdW_ptr[InddFdW++] =  0.0;
			*dFdW_ptr[InddFdW++] =  0.0;
			*dFdW_ptr[InddFdW++] =  0.0;

			// *** eq 2 ***
			// var 1
			*dFdW_ptr[InddFdW++] = -u2+alpha;
			*dFdW_ptr[InddFdW++] = -uv;
			*dFdW_ptr[InddFdW++] = -uw;

			// var 2
			*dFdW_ptr[InddFdW++] = -GM3*u;
			*dFdW_ptr[InddFdW++] =  v;
			*dFdW_ptr[InddFdW++] =  w;

			// var 3
			*dFdW_ptr[InddFdW++] = -GM1*v;
			*dFdW_ptr[InddFdW++] =  u;
			*dFdW_ptr[InddFdW++] =  0.0;

			// var 4
			*dFdW_ptr[InddFdW++] = -GM1*w;
			*dFdW_ptr[InddFdW++] =  0.0;
			*dFdW_ptr[InddFdW++] =  u;

			// var 5
			*dFdW_ptr[InddFdW++] =  GM1;
			*dFdW_ptr[InddFdW++] =  0.0;
			*dFdW_ptr[InddFdW++] =  0.0;

			// *** eq 3 ***
			// var 1
			*dFdW_ptr[InddFdW++] = -uv;
			*dFdW_ptr[InddFdW++] = -v2+alpha;
			*dFdW_ptr[InddFdW++] = -vw;

			// var 2
			*dFdW_ptr[InddFdW++] =  v;
			*dFdW_ptr[InddFdW++] = -GM1*u;
			*dFdW_ptr[InddFdW++] =  0.0;

			// var 3
			*dFdW_ptr[InddFdW++] =  u;
			*dFdW_ptr[InddFdW++] = -GM3*v;
			*dFdW_ptr[InddFdW++] =  w;

			// var 4
			*dFdW_ptr[InddFdW++] =  0.0;
			*dFdW_ptr[InddFdW++] = -GM1*w;
			*dFdW_ptr[InddFdW++] =  v;

			// var 5
			*dFdW_ptr[InddFdW++] =  0.0;
			*dFdW_ptr[InddFdW++] =  GM1;
			*dFdW_ptr[InddFdW++] =  0.0;

			// *** eq 4 ***
			// var 1
			*dFdW_ptr[InddFdW++] = -uw;
			*dFdW_ptr[InddFdW++] = -vw;
			*dFdW_ptr[InddFdW++] = -w2+alpha;

			// var 2
			*dFdW_ptr[InddFdW++] =  w;
			*dFdW_ptr[InddFdW++] =  0.0;
			*dFdW_ptr[InddFdW++] = -GM1*u;

			// var 3
			*dFdW_ptr[InddFdW++] =  0.0;
			*dFdW_ptr[InddFdW++] =  w;
			*dFdW_ptr[InddFdW++] = -GM1*v;

			// var 4
			*dFdW_ptr[InddFdW++] =  u;
			*dFdW_ptr[InddFdW++] =  v;
			*dFdW_ptr[InddFdW++] = -GM3*w;

			// var 5
			*dFdW_ptr[InddFdW++] =  0.0;
			*dFdW_ptr[InddFdW++] =  0.0;
			*dFdW_ptr[InddFdW++] =  GM1;

			// *** eq 5 ***
			// var 1
			*dFdW_ptr[InddFdW++] =  u*beta;
			*dFdW_ptr[InddFdW++] =  v*beta;
			*dFdW_ptr[InddFdW++] =  w*beta;

			// var 2
			*dFdW_ptr[InddFdW++] =  H-GM1*u2;
			*dFdW_ptr[InddFdW++] = -GM1*uv;
			*dFdW_ptr[InddFdW++] = -GM1*uw;

			// var 3
			*dFdW_ptr[InddFdW++] = -GM1*uv;
			*dFdW_ptr[InddFdW++] =  H-GM1*v2;
			*dFdW_ptr[InddFdW++] = -GM1*vw;

			// var 4
			*dFdW_ptr[InddFdW++] = -GM1*uw;
			*dFdW_ptr[InddFdW++] = -GM1*vw;
			*dFdW_ptr[InddFdW++] =  H-GM1*w2;

			// var 5
			*dFdW_ptr[InddFdW++] =  GAMMA*u;
			*dFdW_ptr[InddFdW++] =  GAMMA*v;
			*dFdW_ptr[InddFdW++] =  GAMMA*w;

			rho_ptr++; rhou_ptr++; rhov_ptr++; rhow_ptr++; E_ptr++;
			for (i = 0, iMax = Neq*Nvar*DMAX; i < iMax; i++)
				dFdW_ptr[i]++;
		}
	} else if (d == 2) {
		rhov_ptr = &W[NnTotal*2];

		for (n = 0; n < NnTotal; n++) {
			rho = *rho_ptr;
			u   = (*rhou_ptr)/rho;
			v   = (*rhov_ptr)/rho;
			E   = *E_ptr;

			u2 = u*u;
			uv = u*v;
			v2 = v*v;

			V2 = u2+v2;
			p  = GM1*(E-0.5*rho*V2);
			H  = (E+p)/rho;

			alpha = 0.5*GM1*V2;
			beta  = alpha-H;

			InddFdW = 0;
			// *** eq 1 ***
			// var 1
			*dFdW_ptr[InddFdW++] =  0.0;
			*dFdW_ptr[InddFdW++] =  0.0;
			InddFdW += 1;

			// var 2
			*dFdW_ptr[InddFdW++] =  1.0;
			*dFdW_ptr[InddFdW++] =  0.0;
			InddFdW += 1;

			// var 3
			*dFdW_ptr[InddFdW++] =  0.0;
			*dFdW_ptr[InddFdW++] =  1.0;
			InddFdW += 1;

			// var 4
			*dFdW_ptr[InddFdW++] =  0.0;
			*dFdW_ptr[InddFdW++] =  0.0;
			InddFdW += 1;

			// *** eq 2 ***
			// var 1
			*dFdW_ptr[InddFdW++] = -u2+alpha;
			*dFdW_ptr[InddFdW++] = -uv;
			InddFdW += 1;

			// var 2
			*dFdW_ptr[InddFdW++] = -GM3*u;
			*dFdW_ptr[InddFdW++] =  v;
			InddFdW += 1;

			// var 3
			*dFdW_ptr[InddFdW++] = -GM1*v;
			*dFdW_ptr[InddFdW++] =  u;
			InddFdW += 1;

			// var 4
			*dFdW_ptr[InddFdW++] =  GM1;
			*dFdW_ptr[InddFdW++] =  0.0;
			InddFdW += 1;

			// *** eq 3 ***
			// var 1
			*dFdW_ptr[InddFdW++] = -uv;
			*dFdW_ptr[InddFdW++] = -v2+alpha;
			InddFdW += 1;

			// var 2
			*dFdW_ptr[InddFdW++] =  v;
			*dFdW_ptr[InddFdW++] = -GM1*u;
			InddFdW += 1;

			// var 3
			*dFdW_ptr[InddFdW++] =  u;
			*dFdW_ptr[InddFdW++] = -GM3*v;
			InddFdW += 1;

			// var 4
			*dFdW_ptr[InddFdW++] =  0.0;
			*dFdW_ptr[InddFdW++] =  GM1;
			InddFdW += 1;

			// *** eq 4 ***
			// var 1
			*dFdW_ptr[InddFdW++] =  u*beta;
			*dFdW_ptr[InddFdW++] =  v*beta;
			InddFdW += 1;

			// var 2
			*dFdW_ptr[InddFdW++] =  H-GM1*u2;
			*dFdW_ptr[InddFdW++] = -GM1*uv;
			InddFdW += 1;

			// var 3
			*dFdW_ptr[InddFdW++] = -GM1*uv;
			*dFdW_ptr[InddFdW++] =  H-GM1*v2;
			InddFdW += 1;

			// var 4
			*dFdW_ptr[InddFdW++] =  GAMMA*u;
			*dFdW_ptr[InddFdW++] =  GAMMA*v;

			rho_ptr++; rhou_ptr++; rhov_ptr++; E_ptr++;
			for (i = 0, iMax = Neq*Nvar*DMAX; i < iMax; i++)
				dFdW_ptr[i]++;
		}
	} else if (d == 1) {
		for (n = 0; n < NnTotal; n++) {
			rho = *rho_ptr;
			u   = (*rhou_ptr)/rho;
			E   = *E_ptr;

			u2 = u*u;

			V2 = u2;
			p  = GM1*(E-0.5*rho*V2);
			H  = (E+p)/rho;

			alpha = 0.5*GM1*V2;
			beta  = alpha-H;

			InddFdW = 0;
			// *** eq 1 ***
			// var 1
			*dFdW_ptr[InddFdW++] =  0.0;
			InddFdW += 2;

			// var 2
			*dFdW_ptr[InddFdW++] =  1.0;
			InddFdW += 2;

			// var 3
			*dFdW_ptr[InddFdW++] =  0.0;
			InddFdW += 2;

			// *** eq 2 ***
			// var 1
			*dFdW_ptr[InddFdW++] = -u2+alpha;
			InddFdW += 2;

			// var 2
			*dFdW_ptr[InddFdW++] = -GM3*u;
			InddFdW += 2;

			// var 3
			*dFdW_ptr[InddFdW++] =  GM1;
			InddFdW += 2;

			// *** eq 3 ***
			// var 1
			*dFdW_ptr[InddFdW++] =  u*beta;
			InddFdW += 2;

			// var 2
			*dFdW_ptr[InddFdW++] =  H-GM1*u2;
			InddFdW += 2;

			// var 3
			*dFdW_ptr[InddFdW++] =  GAMMA*u;

			rho_ptr++; rhou_ptr++; E_ptr++;
			for (i = 0, iMax = Neq*Nvar*DMAX; i < iMax; i++)
				dFdW_ptr[i]++;
		}
	}
}

void jacobian_flux_LF(const unsigned int Nn, const unsigned int Nel, double *WL, double *WR, double *dnFdW,
                      double *nL, const unsigned int d, const unsigned int Neq, const char side)
{
	// Standard datatypes
	char         sideMaxV;
	unsigned int i, n, eq, var, iMax, Nvar, NnTotal, InddnFdW;
	double       *rhoL_ptr, *rhouL_ptr, *rhovL_ptr, *rhowL_ptr, *EL_ptr,
	             *rhoR_ptr, *rhouR_ptr, *rhovR_ptr, *rhowR_ptr, *ER_ptr,
	             rhoL, rhouL, rhovL, rhowL, EL, rhoL_inv, uL, vL, wL, V2L, VL, pL, cL,
	             rhoR, rhouR, rhovR, rhowR, ER, rhoR_inv, uR, vR, wR, V2R, VR, pR, cR,
	             maxlL, maxlR, maxV,
	             n1, n2, n3, *n_ptr, *dnFdW_ptr[Neq*Neq];

	Nvar    = Neq;
	NnTotal = Nn*Nel;

	rhoL_ptr  = &WL[NnTotal*0];
	rhouL_ptr = &WL[NnTotal*1];
	EL_ptr    = &WL[NnTotal*(d+1)];

	rhoR_ptr  = &WR[NnTotal*0];
	rhouR_ptr = &WR[NnTotal*1];
	ER_ptr    = &WR[NnTotal*(d+1)];

	n_ptr = nL;

	for (eq  = 0; eq  < Neq;  eq++)  {
	for (var = 0; var < Nvar; var++) {
		dnFdW_ptr[eq*Nvar+var] = &dnFdW[(eq*Nvar+var)*NnTotal];
	}}

	if (d == 3) {
		rhovL_ptr = &WL[NnTotal*2];
		rhowL_ptr = &WL[NnTotal*3];

		rhovR_ptr = &WR[NnTotal*2];
		rhowR_ptr = &WR[NnTotal*3];

		for (n = 0; n < NnTotal; n++) {
			InddnFdW = 0;

			// Inner VOLUME
			rhoL  = *rhoL_ptr++;
			rhouL = *rhouL_ptr++;
			rhovL = *rhovL_ptr++;
			rhowL = *rhowL_ptr++;
			EL    = *EL_ptr++;

			rhoL_inv = 1.0/rhoL;
			uL = rhouL*rhoL_inv;
			vL = rhovL*rhoL_inv;
			wL = rhowL*rhoL_inv;

			V2L = uL*uL+vL*vL+wL*wL;
			VL  = sqrt(V2L);

			pL  = GM1*(EL-0.5*rhoL*V2L);
			cL  = sqrt(GAMMA*pL/rhoL);

			// Outer VOLUME
			rhoR  = *rhoR_ptr++;
			rhouR = *rhouR_ptr++;
			rhovR = *rhovR_ptr++;
			rhowR = *rhowR_ptr++;
			ER    = *ER_ptr++;

			rhoR_inv = 1.0/rhoR;
			uR = rhouR*rhoR_inv;
			vR = rhovR*rhoR_inv;
			wR = rhowR*rhoR_inv;

			V2R = uR*uR+vR*vR+wR*wR;
			VR  = sqrt(V2R);

			pR  = GM1*(ER-0.5*rhoR*V2R);
			cR  = sqrt(GAMMA*pR/rhoR);


			maxlL = VL+cL;
			maxlR = VR+cR;

			if (maxlL > maxlR) {
				sideMaxV = 'L';
				maxV = maxlL;
			} else {
				sideMaxV = 'R';
				maxV = maxlR;
			}

			n1 = *n_ptr++;
			n2 = *n_ptr++;
			n3 = *n_ptr++;

			double dFdWn[Nvar*Nvar*d], *dF1dW_ptr, *dF2dW_ptr, *dF3dW_ptr;
			if (side == 'L') {
				// Flux term
				double WLn[Nvar];

				WLn[0] = rhoL; WLn[1] = rhouL; WLn[2] = rhovL; WLn[3] = rhowL; WLn[4] = EL;
				jacobian_flux_inviscid(1,1,WLn,dFdWn,d,Neq);

				dF1dW_ptr = dFdWn;
				dF2dW_ptr = dF1dW_ptr+1;
				dF3dW_ptr = dF2dW_ptr+1;

				InddnFdW = 0;
				for (eq = 0; eq < Neq; eq++) {
				for (var = 0; var < Nvar; var++) {
					*dnFdW_ptr[InddnFdW++] = 0.5*(n1*(*dF1dW_ptr)+n2*(*dF2dW_ptr)+n3*(*dF3dW_ptr));

					dF1dW_ptr += d;
					dF2dW_ptr += d;
					dF3dW_ptr += d;
				}}

				// Dissipation term
				InddnFdW = 0;
				for (eq = 0; eq < Neq; eq++) {
				for (var = 0; var < Nvar; var++) {
					if (var == eq)
						*dnFdW_ptr[InddnFdW] += 0.5*maxV;
					InddnFdW++;
				}}
				if (sideMaxV == 'L') {
					double WRn[Nvar], dudW[Nvar], dvdW[Nvar], dwdW[Nvar], drhodW[Nvar], dpdW[Nvar],
					       dVdW, dcdW, dmaxVdW[Nvar];

					WRn[0] = rhoR; WRn[1] = rhouR; WRn[2] = rhovR; WRn[3] = rhowR; WRn[4] = ER;

					dudW[0] = -uL*rhoL_inv; dudW[1] = rhoL_inv; dudW[2] = 0.0;      dudW[3] = 0.0;      dudW[4] = 0.0;
					dvdW[0] = -vL*rhoL_inv; dvdW[1] = 0.0;      dvdW[2] = rhoL_inv; dvdW[3] = 0.0;      dvdW[4] = 0.0;
					dwdW[0] = -wL*rhoL_inv; dwdW[1] = 0.0;      dwdW[2] = 0.0;      dwdW[3] = rhoL_inv; dwdW[4] = 0.0;

					drhodW[0] = 1.0;     drhodW[1] = 0.0; drhodW[2] = 0.0; drhodW[3] = 0.0; drhodW[4] = 0.0;
					dpdW[0]   = 0.5*V2L; dpdW[1]   = -uL; dpdW[2]   = -vL; dpdW[3]   = -wL; dpdW[4]   = 1.0;

					for (var = 0; var < Nvar; var++) {
						dpdW[var] *= GM1;
						dVdW = 1/VL*(uL*dudW[var]+vL*dvdW[var]+wL*dwdW[var]);
						dcdW = GAMMA/(2.0*cL*rhoL*rhoL)*(dpdW[var]*rhoL-pL*drhodW[var]);
						dmaxVdW[var] = dVdW+dcdW;
					}

					InddnFdW = 0;
					for (eq = 0; eq < Neq; eq++) {
					for (var = 0; var < Nvar; var++) {
						*dnFdW_ptr[InddnFdW++] += 0.5*dmaxVdW[var]*(WLn[eq]-WRn[eq]);
					}}
				}
			} else {
				// Flux term
				double WRn[Nvar];

				WRn[0] = rhoR; WRn[1] = rhouR; WRn[2] = rhovR; WRn[3] = rhowR; WRn[4] = ER;
				jacobian_flux_inviscid(1,1,WRn,dFdWn,d,Neq);

				dF1dW_ptr = dFdWn;
				dF2dW_ptr = dF1dW_ptr+1;
				dF3dW_ptr = dF2dW_ptr+1;

				InddnFdW = 0;
				for (eq = 0; eq < Neq; eq++) {
				for (var = 0; var < Nvar; var++) {
					*dnFdW_ptr[InddnFdW++] = 0.5*(n1*(*dF1dW_ptr)+n2*(*dF2dW_ptr)+n3*(*dF3dW_ptr));

					dF1dW_ptr += d;
					dF2dW_ptr += d;
					dF3dW_ptr += d;
				}}

				// Dissipation term
				InddnFdW = 0;
				for (eq = 0; eq < Neq; eq++) {
				for (var = 0; var < Nvar; var++) {
					if (var == eq)
						*dnFdW_ptr[InddnFdW] -= 0.5*maxV;
					InddnFdW++;
				}}
				if (sideMaxV == 'R') {
					double WLn[Nvar], dudW[Nvar], dvdW[Nvar], dwdW[Nvar], drhodW[Nvar], dpdW[Nvar],
					       dVdW, dcdW, dmaxVdW[Nvar];

					WLn[0] = rhoL; WLn[1] = rhouL; WLn[2] = rhovL; WLn[3] = rhowL; WLn[4] = EL;

					dudW[0] = -uR*rhoR_inv; dudW[1] = rhoR_inv; dudW[2] = 0.0;      dudW[3] = 0.0;      dudW[4] = 0.0;
					dvdW[0] = -vR*rhoR_inv; dvdW[1] = 0.0;      dvdW[2] = rhoR_inv; dvdW[3] = 0.0;      dvdW[4] = 0.0;
					dwdW[0] = -wR*rhoR_inv; dwdW[1] = 0.0;      dwdW[2] = 0.0;      dwdW[3] = rhoR_inv; dwdW[4] = 0.0;

					drhodW[0] = 1.0;     drhodW[1] = 0.0; drhodW[2] = 0.0; drhodW[3] = 0.0; drhodW[4] = 0.0;
					dpdW[0]   = 0.5*V2R; dpdW[1]   = -uR; dpdW[2]   = -vR; dpdW[3]   = -wR; dpdW[4]   = 1.0;

					for (var = 0; var < Nvar; var++) {
						dpdW[var] *= GM1;
						dVdW = 1/VR*(uR*dudW[var]+vR*dvdW[var]+wR*dwdW[var]);
						dcdW = GAMMA/(2.0*cR*rhoR*rhoR)*(dpdW[var]*rhoR-pR*drhodW[var]);
						dmaxVdW[var] = dVdW+dcdW;
					}

					InddnFdW = 0;
					for (eq = 0; eq < Neq; eq++) {
					for (var = 0; var < Nvar; var++) {
						*dnFdW_ptr[InddnFdW++] += 0.5*dmaxVdW[var]*(WLn[eq]-WRn[eq]);
					}}
				}
			}

			for (i = 0, iMax = Neq*Nvar; i < iMax; i++)
				dnFdW_ptr[i]++;
		}
	} else if (d == 2) {
		rhovL_ptr = &WL[NnTotal*2];
		rhovR_ptr = &WR[NnTotal*2];

		for (n = 0; n < NnTotal; n++) {
			InddnFdW = 0;

			// Inner VOLUME
			rhoL  = *rhoL_ptr++;
			rhouL = *rhouL_ptr++;
			rhovL = *rhovL_ptr++;
			EL    = *EL_ptr++;

			rhoL_inv = 1.0/rhoL;
			uL = rhouL*rhoL_inv;
			vL = rhovL*rhoL_inv;

			V2L = uL*uL+vL*vL;
			VL  = sqrt(V2L);

			pL  = GM1*(EL-0.5*rhoL*V2L);
			cL  = sqrt(GAMMA*pL/rhoL);

			// Outer VOLUME
			rhoR  = *rhoR_ptr++;
			rhouR = *rhouR_ptr++;
			rhovR = *rhovR_ptr++;
			ER    = *ER_ptr++;

			rhoR_inv = 1.0/rhoR;
			uR = rhouR*rhoR_inv;
			vR = rhovR*rhoR_inv;

			V2R = uR*uR+vR*vR;
			VR  = sqrt(V2R);

			pR  = GM1*(ER-0.5*rhoR*V2R);
			cR  = sqrt(GAMMA*pR/rhoR);


			maxlL = VL+cL;
			maxlR = VR+cR;

			if (maxlL > maxlR) {
				sideMaxV = 'L';
				maxV = maxlL;
			} else {
				sideMaxV = 'R';
				maxV = maxlR;
			}

			n1 = *n_ptr++;
			n2 = *n_ptr++;

			double dFdWn[Nvar*Nvar*d], *dF1dW_ptr, *dF2dW_ptr;
			if (side == 'L') {
				// Flux term
				double WLn[Nvar];

				WLn[0] = rhoL; WLn[1] = rhouL; WLn[2] = rhovL; WLn[3] = EL;
				jacobian_flux_inviscid(1,1,WLn,dFdWn,d,Neq);

				dF1dW_ptr = dFdWn;
				dF2dW_ptr = dF1dW_ptr+1;

				InddnFdW = 0;
				for (eq = 0; eq < Neq; eq++) {
				for (var = 0; var < Nvar; var++) {
					*dnFdW_ptr[InddnFdW++] = 0.5*(n1*(*dF1dW_ptr)+n2*(*dF2dW_ptr));

					dF1dW_ptr += d;
					dF2dW_ptr += d;
				}}

				// Dissipation term
				InddnFdW = 0;
				for (eq = 0; eq < Neq; eq++) {
				for (var = 0; var < Nvar; var++) {
					if (var == eq)
						*dnFdW_ptr[InddnFdW] += 0.5*maxV;
					InddnFdW++;
				}}
				if (sideMaxV == 'L') {
					double WRn[Nvar], dudW[Nvar], dvdW[Nvar], drhodW[Nvar], dpdW[Nvar],
					       dVdW, dcdW, dmaxVdW[Nvar];

					WRn[0] = rhoR; WRn[1] = rhouR; WRn[2] = rhovR; WRn[3] = ER;

					dudW[0] = -uL*rhoL_inv; dudW[1] = rhoL_inv; dudW[2] = 0.0;      dudW[3] = 0.0;
					dvdW[0] = -vL*rhoL_inv; dvdW[1] = 0.0;      dvdW[2] = rhoL_inv; dvdW[3] = 0.0;

					drhodW[0] = 1.0;     drhodW[1] = 0.0; drhodW[2] = 0.0; drhodW[3] = 0.0;
					dpdW[0]   = 0.5*V2L; dpdW[1]   = -uL; dpdW[2]   = -vL; dpdW[3]   = 1.0;

					for (var = 0; var < Nvar; var++) {
						dpdW[var] *= GM1;
						dVdW = 1/VL*(uL*dudW[var]+vL*dvdW[var]);
						dcdW = GAMMA/(2.0*cL*rhoL*rhoL)*(dpdW[var]*rhoL-pL*drhodW[var]);
						dmaxVdW[var] = dVdW+dcdW;
					}

					InddnFdW = 0;
					for (eq = 0; eq < Neq; eq++) {
					for (var = 0; var < Nvar; var++) {
						*dnFdW_ptr[InddnFdW++] += 0.5*dmaxVdW[var]*(WLn[eq]-WRn[eq]);
					}}
				}
			} else {
				// Flux term
				double WRn[Nvar];

				WRn[0] = rhoR; WRn[1] = rhouR; WRn[2] = rhovR; WRn[3] = ER;
				jacobian_flux_inviscid(1,1,WRn,dFdWn,d,Neq);

				dF1dW_ptr = dFdWn;
				dF2dW_ptr = dF1dW_ptr+1;

				InddnFdW = 0;
				for (eq = 0; eq < Neq; eq++) {
				for (var = 0; var < Nvar; var++) {
					*dnFdW_ptr[InddnFdW++] = 0.5*(n1*(*dF1dW_ptr)+n2*(*dF2dW_ptr));

					dF1dW_ptr += d;
					dF2dW_ptr += d;
				}}

				// Dissipation term
				InddnFdW = 0;
				for (eq = 0; eq < Neq; eq++) {
				for (var = 0; var < Nvar; var++) {
					if (var == eq)
						*dnFdW_ptr[InddnFdW] -= 0.5*maxV;
					InddnFdW++;
				}}
				if (sideMaxV == 'R') {
					double WLn[Nvar], dudW[Nvar], dvdW[Nvar], drhodW[Nvar], dpdW[Nvar],
					       dVdW, dcdW, dmaxVdW[Nvar];

					WLn[0] = rhoL; WLn[1] = rhouL; WLn[2] = rhovL; WLn[3] = EL;

					dudW[0] = -uR*rhoR_inv; dudW[1] = rhoR_inv; dudW[2] = 0.0;      dudW[3] = 0.0;
					dvdW[0] = -vR*rhoR_inv; dvdW[1] = 0.0;      dvdW[2] = rhoR_inv; dvdW[3] = 0.0;

					drhodW[0] = 1.0;     drhodW[1] = 0.0; drhodW[2] = 0.0; drhodW[3] = 0.0;
					dpdW[0]   = 0.5*V2R; dpdW[1]   = -uR; dpdW[2]   = -vR; dpdW[3]   = 1.0;

					for (var = 0; var < Nvar; var++) {
						dpdW[var] *= GM1;
						dVdW = 1/VR*(uR*dudW[var]+vR*dvdW[var]);
						dcdW = GAMMA/(2.0*cR*rhoR*rhoR)*(dpdW[var]*rhoR-pR*drhodW[var]);
						dmaxVdW[var] = dVdW+dcdW;
					}

					InddnFdW = 0;
					for (eq = 0; eq < Neq; eq++) {
					for (var = 0; var < Nvar; var++) {
						*dnFdW_ptr[InddnFdW++] += 0.5*dmaxVdW[var]*(WLn[eq]-WRn[eq]);
					}}
				}
			}

			for (i = 0, iMax = Neq*Nvar; i < iMax; i++)
				dnFdW_ptr[i]++;
		}
	} else if (d == 1) {
		for (n = 0; n < NnTotal; n++) {
			InddnFdW = 0;

			// Inner VOLUME
			rhoL  = *rhoL_ptr++;
			rhouL = *rhouL_ptr++;
			EL    = *EL_ptr++;

			rhoL_inv = 1.0/rhoL;
			uL = rhouL*rhoL_inv;

			V2L = uL*uL;
			VL  = sqrt(V2L);

			pL  = GM1*(EL-0.5*rhoL*V2L);
			cL  = sqrt(GAMMA*pL/rhoL);

			// Outer VOLUME
			rhoR  = *rhoR_ptr++;
			rhouR = *rhouR_ptr++;
			ER    = *ER_ptr++;

			rhoR_inv = 1.0/rhoR;
			uR = rhouR*rhoR_inv;

			V2R = uR*uR;
			VR  = sqrt(V2R);

			pR  = GM1*(ER-0.5*rhoR*V2R);
			cR  = sqrt(GAMMA*pR/rhoR);


			maxlL = VL+cL;
			maxlR = VR+cR;

			if (maxlL > maxlR) {
				sideMaxV = 'L';
				maxV = maxlL;
			} else {
				sideMaxV = 'R';
				maxV = maxlR;
			}

			n1 = *n_ptr++;

			double dFdWn[Nvar*Nvar*d], *dF1dW_ptr;
			if (side == 'L') {
				// Flux term
				double WLn[Nvar];

				WLn[0] = rhoL; WLn[1] = rhouL; WLn[2] = EL;
				jacobian_flux_inviscid(1,1,WLn,dFdWn,d,Neq);

				dF1dW_ptr = dFdWn;

				InddnFdW = 0;
				for (eq = 0; eq < Neq; eq++) {
				for (var = 0; var < Nvar; var++) {
					*dnFdW_ptr[InddnFdW++] = 0.5*(n1*(*dF1dW_ptr));

					dF1dW_ptr += d;
				}}

				// Dissipation term
				InddnFdW = 0;
				for (eq = 0; eq < Neq; eq++) {
				for (var = 0; var < Nvar; var++) {
					if (var == eq)
						*dnFdW_ptr[InddnFdW] += 0.5*maxV;
					InddnFdW++;
				}}
				if (sideMaxV == 'L') {
					double WRn[Nvar], dudW[Nvar], drhodW[Nvar], dpdW[Nvar],
					       dVdW, dcdW, dmaxVdW[Nvar];

					WRn[0] = rhoR; WRn[1] = rhouR; WRn[2] = ER;

					dudW[0] = -uL*rhoL_inv; dudW[1] = rhoL_inv; dudW[2] = 0.0;

					drhodW[0] = 1.0;     drhodW[1] = 0.0; drhodW[2] = 0.0;
					dpdW[0]   = 0.5*V2L; dpdW[1]   = -uL; dpdW[2]   = 1.0;

					for (var = 0; var < Nvar; var++) {
						dpdW[var] *= GM1;
						dVdW = 1/VL*(uL*dudW[var]);
						dcdW = GAMMA/(2.0*cL*rhoL*rhoL)*(dpdW[var]*rhoL-pL*drhodW[var]);
						dmaxVdW[var] = dVdW+dcdW;
					}

					InddnFdW = 0;
					for (eq = 0; eq < Neq; eq++) {
					for (var = 0; var < Nvar; var++) {
						*dnFdW_ptr[InddnFdW++] += 0.5*dmaxVdW[var]*(WLn[eq]-WRn[eq]);
					}}
				}
			} else {
				// Flux term
				double WRn[Nvar];

				WRn[0] = rhoR; WRn[1] = rhouR; WRn[2] = ER;
				jacobian_flux_inviscid(1,1,WRn,dFdWn,d,Neq);

				dF1dW_ptr = dFdWn;

				InddnFdW = 0;
				for (eq = 0; eq < Neq; eq++) {
				for (var = 0; var < Nvar; var++) {
					*dnFdW_ptr[InddnFdW++] = 0.5*(n1*(*dF1dW_ptr));

					dF1dW_ptr += d;
				}}

				// Dissipation term
				InddnFdW = 0;
				for (eq = 0; eq < Neq; eq++) {
				for (var = 0; var < Nvar; var++) {
					if (var == eq)
						*dnFdW_ptr[InddnFdW] -= 0.5*maxV;
					InddnFdW++;
				}}
				if (sideMaxV == 'R') {
					double WLn[Nvar], dudW[Nvar], drhodW[Nvar], dpdW[Nvar],
					       dVdW, dcdW, dmaxVdW[Nvar];

					WLn[0] = rhoL; WLn[1] = rhouL; WLn[2] = EL;

					dudW[0] = -uR*rhoR_inv; dudW[1] = rhoR_inv; dudW[2] = 0.0;

					drhodW[0] = 1.0;     drhodW[1] = 0.0; drhodW[2] = 0.0;
					dpdW[0]   = 0.5*V2R; dpdW[1]   = -uR; dpdW[2]   = 1.0;

					for (var = 0; var < Nvar; var++) {
						dpdW[var] *= GM1;
						dVdW = 1/VR*(uR*dudW[var]);
						dcdW = GAMMA/(2.0*cR*rhoR*rhoR)*(dpdW[var]*rhoR-pR*drhodW[var]);
						dmaxVdW[var] = dVdW+dcdW;
					}

					InddnFdW = 0;
					for (eq = 0; eq < Neq; eq++) {
					for (var = 0; var < Nvar; var++) {
						*dnFdW_ptr[InddnFdW++] += 0.5*dmaxVdW[var]*(WLn[eq]-WRn[eq]);
					}}
				}
			}

			for (i = 0, iMax = Neq*Nvar; i < iMax; i++)
				dnFdW_ptr[i]++;
		}
	}
}

void jacobian_flux_Roe(const unsigned int Nn, const unsigned int Nel, double *WL, double *WR, double *dnFdW,
                       double *nL, const unsigned int d, const unsigned int Neq, const char side)
{
	// Standard datatypes
	unsigned int i, n, eq, var, iMax, Nvar, NnTotal, InddnFdW;
	double       *rhoL_ptr, *rhouL_ptr, *rhovL_ptr, *rhowL_ptr, *EL_ptr,
	             *rhoR_ptr, *rhouR_ptr, *rhovR_ptr, *rhowR_ptr, *ER_ptr,
	             rhoL, rhouL, rhovL, rhowL, EL, rhoL_inv, uL, vL, wL, V2L, VL, VnL, pL, cL,
	             rhoR, rhouR, rhovR, rhowR, ER, rhoR_inv, uR, vR, wR, V2R, VR, VnR, pR, cR,
	             n1, n2, n3, *n_ptr, *dnFdW_ptr[Neq*Neq];


	Nvar    = Neq;
	NnTotal = Nn*Nel;

	rhoL_ptr  = &WL[NnTotal*0];
	rhouL_ptr = &WL[NnTotal*1];
	EL_ptr    = &WL[NnTotal*(d+1)];

	rhoR_ptr  = &WR[NnTotal*0];
	rhouR_ptr = &WR[NnTotal*1];
	ER_ptr    = &WR[NnTotal*(d+1)];

	n_ptr = nL;

	for (eq  = 0; eq  < Neq;  eq++)  {
	for (var = 0; var < Nvar; var++) {
		dnFdW_ptr[eq*Nvar+var] = &dnFdW[(eq*Nvar+var)*NnTotal];
	}}

	if (d == 3) {
		rhovL_ptr = &WL[NnTotal*2];
		rhowL_ptr = &WL[NnTotal*3];

		rhovR_ptr = &WR[NnTotal*2];
		rhowR_ptr = &WR[NnTotal*3];

		for (n = 0; n < NnTotal; n++) {
			n1 = *n_ptr++;
			n2 = *n_ptr++;
			n3 = *n_ptr++;

			// Inner VOLUME
			rhoL  = *rhoL_ptr++;
			rhouL = *rhouL_ptr++;
			rhovL = *rhovL_ptr++;
			rhowL = *rhowL_ptr++;
			EL    = *EL_ptr++;

			rhoL_inv = 1.0/rhoL;
			uL = rhouL*rhoL_inv;
			vL = rhovL*rhoL_inv;
			wL = rhowL*rhoL_inv;

			V2L = uL*uL+vL*vL+wL*wL;
			VL  = sqrt(V2L);
			VnL = n1*uL+n2*vL+n3*wL;

			pL  = GM1*(EL-0.5*rhoL*V2L);
			cL  = sqrt(GAMMA*pL/rhoL);

			// Outer VOLUME
			rhoR  = *rhoR_ptr++;
			rhouR = *rhouR_ptr++;
			rhovR = *rhovR_ptr++;
			rhowR = *rhowR_ptr++;
			ER    = *ER_ptr++;

			rhoR_inv = 1.0/rhoR;
			uR = rhouR*rhoR_inv;
			vR = rhovR*rhoR_inv;
			wR = rhowR*rhoR_inv;

			V2R = uR*uR+vR*vR+wR*wR;
			VR  = sqrt(V2R);
			VnR = n1*uR+n2*vR+n3*wR;

			pR  = GM1*(ER-0.5*rhoR*V2R);
			cR  = sqrt(GAMMA*pR/rhoR);



			double dudW[Nvar], dvdW[Nvar], dwdW[Nvar], drhodW[Nvar], dEdW[Nvar], dpdW[Nvar],
				   rhoVn, dVndW, drhoVndW,
				   dnF1dW[Nvar], dnF2dW[Nvar], dnF3dW[Nvar], dnF4dW[Nvar], dnF5dW[Nvar];

			if (side == 'L') {
				// Flux term
				rhoVn = rhoL*VnL;

				dudW[0] = -uL*rhoL_inv; dudW[1] = rhoL_inv; dudW[2] = 0.0;      dudW[3] = 0.0;      dudW[4] = 0.0;
				dvdW[0] = -vL*rhoL_inv; dvdW[1] = 0.0;      dvdW[2] = rhoL_inv; dvdW[3] = 0.0;      dvdW[4] = 0.0;
				dwdW[0] = -wL*rhoL_inv; dwdW[1] = 0.0;      dwdW[2] = 0.0;      dwdW[3] = rhoL_inv; dwdW[4] = 0.0;

				drhodW[0] = 1.0;     drhodW[1] = 0.0; drhodW[2] = 0.0; drhodW[3] = 0.0; drhodW[4] = 0.0;
				dEdW[0]   = 0.0;     dEdW[1]   = 0.0; dEdW[2]   = 0.0; dEdW[3]   = 0.0; dEdW[4]   = 1.0;
				dpdW[0]   = 0.5*V2L; dpdW[1]   = -uL; dpdW[2]   = -vL; dpdW[3]   = -wL; dpdW[4]   = 1.0;

				for (var = 0; var < Nvar; var++) {
					dpdW[var] *= GM1;

					dVndW    = n1*dudW[var]+n2*dvdW[var]+n3*dwdW[var];
					drhoVndW = drhodW[var]*VnL + rhoL*dVndW;

					dnF1dW[var] = drhoVndW;
					dnF2dW[var] = drhoVndW*uL + rhoVn*dudW[var] + n1*dpdW[var];
					dnF3dW[var] = drhoVndW*vL + rhoVn*dvdW[var] + n2*dpdW[var];
					dnF4dW[var] = drhoVndW*wL + rhoVn*dwdW[var] + n3*dpdW[var];
					dnF5dW[var] = dVndW*(EL+pL) * VnL*(dEdW[var]+dpdW[var]);
				}

				InddnFdW = 0;
				for (var = 0; var < Nvar; var++) *dnFdW_ptr[InddnFdW++] = 0.5*dnF1dW[var];
				for (var = 0; var < Nvar; var++) *dnFdW_ptr[InddnFdW++] = 0.5*dnF2dW[var];
				for (var = 0; var < Nvar; var++) *dnFdW_ptr[InddnFdW++] = 0.5*dnF3dW[var];
				for (var = 0; var < Nvar; var++) *dnFdW_ptr[InddnFdW++] = 0.5*dnF4dW[var];
				for (var = 0; var < Nvar; var++) *dnFdW_ptr[InddnFdW++] = 0.5*dnF5dW[var];
			} else {
				pL = cR = VnR = VR = cL = VL;
			}

			for (i = 0, iMax = Neq*Nvar; i < iMax; i++)
				dnFdW_ptr[i]++;
		}
	} else if (d == 2) {
	} else if (d == 1) {
	}

}
