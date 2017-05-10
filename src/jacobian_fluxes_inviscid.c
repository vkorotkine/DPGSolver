// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "jacobian_fluxes_inviscid.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Parameters.h"

#include "fluxes_structs.h"

/*
 *	Purpose:
 *		Compute inviscid flux jacobians from input W in conservative form. Also returns the flux if applicable.
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

void jacobian_flux_inviscid(struct S_FLUX *const FLUXDATA)
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

	unsigned int const d   = FLUXDATA->d,
	                   Neq = d+2,
	                   Nn  = FLUXDATA->Nn,
	                   Nel = FLUXDATA->Nel;

	double const *const W    = FLUXDATA->W;
	double       *const F    = FLUXDATA->F,
	             *const dFdW = FLUXDATA->dFdW;

	// Standard datatypes
	unsigned int i, n, eq, var, dim, iMax, Nvar, NnTotal, InddFdW;
	double       rho, u, v, w, u2, uv, uw, v2, vw, w2, V2, E, p, H, alpha, beta, *dFdW_ptr[DMAX*Neq*Neq];
	const double *rho_ptr, *rhou_ptr, *rhov_ptr, *rhow_ptr, *E_ptr;

	Nvar    = Neq;
	NnTotal = Nn*Nel;

	rho_ptr  = &W[NnTotal*0];
	rhou_ptr = &W[NnTotal*1];
	E_ptr    = &W[NnTotal*(d+1)];

	double *F_ptr[DMAX*Neq];
	if (F != NULL) {
		for (eq  = 0; eq  < Neq;  eq++)  {
		for (dim = 0; dim < d;    dim++) {
			F_ptr[eq*DMAX+dim] = &F[(eq*d+dim)*NnTotal];
		}}
	}

	for (eq  = 0; eq  < Neq;  eq++)  {
	for (var = 0; var < Nvar; var++) {
	for (dim = 0; dim < d;    dim++) {
		dFdW_ptr[(eq*Nvar+var)*DMAX+dim] = &dFdW[((eq*Nvar+var)*d+dim)*NnTotal];
	}}}

	if (d == 3) {
		rhov_ptr = &W[NnTotal*2];
		rhow_ptr = &W[NnTotal*3];

		for (n = 0; n < NnTotal; n++) {
			rho  = *rho_ptr;
			double const rhou = *rhou_ptr,
			             rhov = *rhov_ptr,
			             rhow = *rhow_ptr;
			u   = rhou/rho;
			v   = rhov/rho;
			w   = rhow/rho;
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

			if (F != NULL) {
				size_t IndF = 0;
				// eq 1
				*F_ptr[IndF++] = rhou;
				*F_ptr[IndF++] = rhov;
				*F_ptr[IndF++] = rhow;

				// eq 2
				*F_ptr[IndF++] = rhou*u + p;
				*F_ptr[IndF++] = rhou*v;
				*F_ptr[IndF++] = rhou*w;

				// eq 3
				*F_ptr[IndF++] = rhov*u;
				*F_ptr[IndF++] = rhov*v + p;
				*F_ptr[IndF++] = rhov*w;

				// eq 4
				*F_ptr[IndF++] = rhow*u;
				*F_ptr[IndF++] = rhow*v;
				*F_ptr[IndF++] = rhow*w + p;

				// eq 5
				*F_ptr[IndF++] = (E+p)*u;
				*F_ptr[IndF++] = (E+p)*v;
				*F_ptr[IndF++] = (E+p)*w;

				for (i = 0, iMax = Neq*DMAX; i < iMax; i++)
					F_ptr[i]++;
			}

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
			double const rhou = *rhou_ptr,
			             rhov = *rhov_ptr;
			u   = rhou/rho;
			v   = rhov/rho;
			E   = *E_ptr;

			u2 = u*u;
			uv = u*v;
			v2 = v*v;

			V2 = u2+v2;
			p  = GM1*(E-0.5*rho*V2);
			H  = (E+p)/rho;

			alpha = 0.5*GM1*V2;
			beta  = alpha-H;

			if (F != NULL) {
				size_t IndF = 0;
				// eq 1
				*F_ptr[IndF++] = rhou;
				*F_ptr[IndF++] = rhov;
				IndF += 1;

				// eq 2
				*F_ptr[IndF++] = rhou*u + p;
				*F_ptr[IndF++] = rhou*v;
				IndF += 1;

				// eq 3
				*F_ptr[IndF++] = rhov*u;
				*F_ptr[IndF++] = rhov*v + p;
				IndF += 1;

				// eq 4
				*F_ptr[IndF++] = (E+p)*u;
				*F_ptr[IndF++] = (E+p)*v;

				for (i = 0, iMax = Neq*DMAX; i < iMax; i++)
					F_ptr[i]++;
			}

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
			double const rhou = *rhou_ptr;
			u   = rhou/rho;
			E   = *E_ptr;

			u2 = u*u;

			V2 = u2;
			p  = GM1*(E-0.5*rho*V2);
			H  = (E+p)/rho;

			alpha = 0.5*GM1*V2;
			beta  = alpha-H;

			if (F != NULL) {
				size_t IndF = 0;
				// eq 1
				*F_ptr[IndF++] = rhou;
				IndF += 2;

				// eq 2
				*F_ptr[IndF++] = rhou*u + p;
				IndF += 2;

				// eq 3
				*F_ptr[IndF++] = (E+p)*u;
				IndF += 2;

				for (i = 0, iMax = Neq*DMAX; i < iMax; i++)
					F_ptr[i]++;
			}

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

void jacobian_flux_LF(const unsigned int Nn, const unsigned int Nel, const double *WL, const double *WR, double *dnFdW,
                      const double *nL, const unsigned int d, const unsigned int Neq, const char side)
{
	// Standard datatypes
	char         sideMaxV;
	unsigned int i, n, eq, var, iMax, Nvar, NnTotal, InddnFdW;
	double       rhoL, rhouL, rhovL, rhowL, EL, rhoL_inv, uL, vL, wL, V2L, VL, pL, cL,
	             rhoR, rhouR, rhovR, rhowR, ER, rhoR_inv, uR, vR, wR, V2R, VR, pR, cR,
	             maxlL, maxlR, maxV, n1, n2, n3, *dnFdW_ptr[Neq*Neq];
	const double *rhoL_ptr, *rhouL_ptr, *rhovL_ptr, *rhowL_ptr, *EL_ptr,
	             *rhoR_ptr, *rhouR_ptr, *rhovR_ptr, *rhowR_ptr, *ER_ptr, *n_ptr;

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

	struct S_FLUX *const FLUXDATA = malloc(sizeof *FLUXDATA); // free
	FLUXDATA->d   = d;
	FLUXDATA->Nn  = 1;
	FLUXDATA->Nel = 1;
	FLUXDATA->F   = NULL;

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
				FLUXDATA->W    = WLn;
				FLUXDATA->dFdW = dFdWn;
				jacobian_flux_inviscid(FLUXDATA);

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
				FLUXDATA->W    = WRn;
				FLUXDATA->dFdW = dFdWn;
				jacobian_flux_inviscid(FLUXDATA);

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
				FLUXDATA->W    = WLn;
				FLUXDATA->dFdW = dFdWn;
				jacobian_flux_inviscid(FLUXDATA);

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
				FLUXDATA->W    = WRn;
				FLUXDATA->dFdW = dFdWn;
				jacobian_flux_inviscid(FLUXDATA);

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
				FLUXDATA->W    = WLn;
				FLUXDATA->dFdW = dFdWn;
				jacobian_flux_inviscid(FLUXDATA);

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
				FLUXDATA->W    = WRn;
				FLUXDATA->dFdW = dFdWn;
				jacobian_flux_inviscid(FLUXDATA);

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
	free(FLUXDATA);
}

void jacobian_flux_Roe(const unsigned int Nn, const unsigned int Nel, const double *WL, const double *WR, double *dnFdW,
                       const double *nL, const unsigned int d, const unsigned int Neq, const char side)
{
	// Standard datatypes
	unsigned int i, n, eq, var, iMax, Nvar, NnTotal, InddnFdW, case_l1, case_l5;
	double       rhoL, rhouL, rhovL, rhowL, EL, rhoL_inv, rhoLr_inv, uL, vL, wL, V2L, VnL, pL, HL,
	             rhoR, rhouR, rhovR, rhowR, ER, rhoR_inv, rhoRr_inv, uR, vR, wR, V2R, VnR, pR, HR,
	             n1, n2, n3, *dnFdW_ptr[Neq*Neq],
	             r, rP1, rho, u, v, w, H, Vn, V2, c2, c, Den, mult, drho, drhou, drhov, drhow, dE, dp, dVn,
	             sign_l1, sign_l234, sign_l5, l1, l234, l5, l1L, l5R, lc1, lc2, disInter1, disInter2;
	double       duLdW[Neq], dvLdW[Neq], dwLdW[Neq], drhoLdW[Neq], dVnLdW[Neq], dEdW[Neq], dpdW[Neq],
	             duRdW[Neq], dvRdW[Neq], dwRdW[Neq], drhoRdW[Neq], dVnRdW[Neq],
	             rhoVn, dVndW, drhoVndW, dcdW, drhoudW, drhovdW, drhowdW,
	             dudW[Neq], dvdW[Neq], dwdW[Neq], dHdW[Neq], drhodW[Neq],
	             dl1dW, dl234dW, dl5dW, dlc1dW, dlc2dW, ddisInter1dW, ddisInter2dW,
	             dnF1dW[Neq], dnF2dW[Neq], dnF3dW[Neq], dnF4dW[Neq], dnF5dW[Neq],
	             ddis1dW[Neq], ddis2dW[Neq], ddis3dW[Neq], ddis4dW[Neq], ddis5dW[Neq];
	const double *rhoL_ptr, *rhouL_ptr, *rhovL_ptr, *rhowL_ptr, *EL_ptr,
	             *rhoR_ptr, *rhouR_ptr, *rhovR_ptr, *rhowR_ptr, *ER_ptr, *n_ptr;


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

			rhoL_inv  = 1.0/rhoL;
			rhoLr_inv = sqrt(rhoL_inv);

			uL = rhouL*rhoL_inv;
			vL = rhovL*rhoL_inv;
			wL = rhowL*rhoL_inv;

			V2L = uL*uL+vL*vL+wL*wL;
			VnL = n1*uL+n2*vL+n3*wL;

			pL  = GM1*(EL-0.5*rhoL*V2L);
			HL  = (EL+pL)*rhoL_inv;

			// Outer VOLUME
			rhoR  = *rhoR_ptr++;
			rhouR = *rhouR_ptr++;
			rhovR = *rhovR_ptr++;
			rhowR = *rhowR_ptr++;
			ER    = *ER_ptr++;

			rhoR_inv  = 1.0/rhoR;
			rhoRr_inv = sqrt(rhoR_inv);

			uR = rhouR*rhoR_inv;
			vR = rhovR*rhoR_inv;
			wR = rhowR*rhoR_inv;

			V2R = uR*uR+vR*vR+wR*wR;
			VnR = n1*uR+n2*vR+n3*wR;

			pR  = GM1*(ER-0.5*rhoR*V2R);
			HR  = (ER+pR)*rhoR_inv;

			// Roe-averaged states
			r = sqrt(rhoR/rhoL);
			rP1 = r+1;

			rho = r*rhoL;
			u   = (r*uR+uL)/rP1;
			v   = (r*vR+vL)/rP1;
			w   = (r*wR+wL)/rP1;
			H   = (r*HR+HL)/rP1;
			Vn  = n1*u+n2*v+n3*w;
			V2  = u*u+v*v+w*w;
			c2  = GM1*(H-0.5*V2);
			c   = sqrt(c2);

			Den  = sqrt(rhoL) + sqrt(rhoR);

			// Compute eigenvalues (with entropy fix)
			l1L = VnL-c;
			l1  = Vn-c;

			sign_l1 = 1.0;
			if (fabs(l1L) < fabs(l1)) {
				case_l1 = 0;
				if (l1L < 0.0)
					sign_l1 = -1.0;
				l1 = sign_l1*l1L;
			} else {
				case_l1 = 1;
				if (l1 < 0.0)
					sign_l1 = -1.0;
				l1 = sign_l1*l1;
			}

			l5R = VnR+c;
			l5  = Vn+c;

			sign_l5 = 1.0;
			if (fabs(l5R) > fabs(l5)) {
				case_l5 = 0;
				if (l5R < 0.0)
					sign_l5 = -1.0;
				l5 = sign_l5*l5R;
			} else {
				case_l5 = 1;
				if (l5 < 0.0)
					sign_l5 = -1.0;
				l5 = sign_l5*l5;
			}

			sign_l234 = 1.0;
			if (Vn < 0.0)
				sign_l234 = -1.0;
			l234 = sign_l234*(Vn);


			drho  = rhoR-rhoL;
			drhou = rhoR*uR-rhoL*uL;
			drhov = rhoR*vR-rhoL*vL;
			drhow = rhoR*wR-rhoL*wL;
			dE    = ER-EL;
			dp    = pR-pL;
			dVn   = VnR-VnL;

			lc1 = 0.5*(l5+l1) - l234;
			lc2 = 0.5*(l5-l1);

			disInter1 = lc1*dp/(c*c) + lc2*rho*dVn/c;
			disInter2 = lc1*rho*dVn  + lc2*dp/c;

			if (side == 'L') {
				// Flux term
				rhoVn = rhoL*VnL;

				duLdW[0] = -uL*rhoL_inv; duLdW[1] = rhoL_inv; duLdW[2] = 0.0;      duLdW[3] = 0.0;      duLdW[4] = 0.0;
				dvLdW[0] = -vL*rhoL_inv; dvLdW[1] = 0.0;      dvLdW[2] = rhoL_inv; dvLdW[3] = 0.0;      dvLdW[4] = 0.0;
				dwLdW[0] = -wL*rhoL_inv; dwLdW[1] = 0.0;      dwLdW[2] = 0.0;      dwLdW[3] = rhoL_inv; dwLdW[4] = 0.0;

				drhoLdW[0] = 1.0;     drhoLdW[1] = 0.0; drhoLdW[2] = 0.0; drhoLdW[3] = 0.0; drhoLdW[4] = 0.0;
				dEdW[0]    = 0.0;     dEdW[1]    = 0.0; dEdW[2]    = 0.0; dEdW[3]    = 0.0; dEdW[4]    = 1.0;
				dpdW[0]    = 0.5*V2L; dpdW[1]    = -uL; dpdW[2]    = -vL; dpdW[3]    = -wL; dpdW[4]    = 1.0;

				for (var = 0; var < Nvar; var++) {
					dpdW[var] *= GM1;

					dVnLdW[var] = n1*duLdW[var]+n2*dvLdW[var]+n3*dwLdW[var];
					drhoVndW    = drhoLdW[var]*VnL + rhoL*dVnLdW[var];

					dnF1dW[var] = drhoVndW;
					dnF2dW[var] = drhoVndW*uL + rhoVn*duLdW[var] + n1*dpdW[var];
					dnF3dW[var] = drhoVndW*vL + rhoVn*dvLdW[var] + n2*dpdW[var];
					dnF4dW[var] = drhoVndW*wL + rhoVn*dwLdW[var] + n3*dpdW[var];
					dnF5dW[var] = dVnLdW[var]*(EL+pL) + VnL*(dEdW[var]+dpdW[var]);
				}

				// Dissipation term
				mult      = rhoLr_inv/Den;
				drhodW[0] = 0.5*rhoR/rho;     drhodW[1] = 0.0;  drhodW[2] = 0.0;  drhodW[3] = 0.0;  drhodW[4] = 0.0;
				dudW[0]   = -0.5*(uL+u)*mult; dudW[1]   = mult; dudW[2]   = 0.0;  dudW[3]   = 0.0;  dudW[4]   = 0.0;
				dvdW[0]   = -0.5*(vL+v)*mult; dvdW[1]   = 0.0;  dvdW[2]   = mult; dvdW[3]   = 0.0;  dvdW[4]   = 0.0;
				dwdW[0]   = -0.5*(wL+w)*mult; dwdW[1]   = 0.0;  dwdW[2]   = 0.0;  dwdW[3]   = mult; dwdW[4]   = 0.0;

				dHdW[0] = -0.5*(HL+H-GM1*V2L)*mult;
				dHdW[1] = -GM1*uL*mult;
				dHdW[2] = -GM1*vL*mult;
				dHdW[3] = -GM1*wL*mult;
				dHdW[4] = GAMMA*mult;

				for (var = 0; var < Nvar; var++) {
					dcdW = 0.5*GM1/c*(dHdW[var]-(u*dudW[var]+v*dvdW[var]+w*dwdW[var]));
					dVndW = n1*dudW[var]+n2*dvdW[var]+n3*dwdW[var];

					if (case_l1)
						dl1dW = sign_l1*(dVndW-dcdW);
					else
						dl1dW = sign_l1*(dVnLdW[var]-dcdW);

					if (case_l5)
						dl5dW = sign_l5*(dVndW+dcdW);
					else
						dl5dW = sign_l5*(dcdW);

					dl234dW = sign_l234*dVndW;

					dlc1dW = 0.5*(dl5dW+dl1dW) - dl234dW;
					dlc2dW = 0.5*(dl5dW-dl1dW);

					ddisInter1dW = (dlc1dW*dp-lc1*dpdW[var])/c2-(2.0*lc1*dp*dcdW)/(c*c2) +
					               (dlc2dW*rho*dVn+lc2*drhodW[var]*dVn-lc2*rho*dVnLdW[var])/c-(lc2*rho*dVn*dcdW)/c2;
					ddisInter2dW = dlc1dW*rho*dVn+lc1*drhodW[var]*dVn-lc1*rho*dVnLdW[var] +
					               (dlc2dW*dp-lc2*dpdW[var])/c-(lc2*dp*dcdW)/c2;

					drhoudW = drhoLdW[var]*uL+rhoL*duLdW[var];
					drhovdW = drhoLdW[var]*vL+rhoL*dvLdW[var];
					drhowdW = drhoLdW[var]*wL+rhoL*dwLdW[var];

					ddis1dW[var] = dl234dW*drho -l234*drhoLdW[var] + ddisInter1dW;
					ddis2dW[var] = dl234dW*drhou-l234*drhoudW      + ddisInter1dW*u  + disInter1*dudW[var] + ddisInter2dW*n1;
					ddis3dW[var] = dl234dW*drhov-l234*drhovdW      + ddisInter1dW*v  + disInter1*dvdW[var] + ddisInter2dW*n2;
					ddis4dW[var] = dl234dW*drhow-l234*drhowdW      + ddisInter1dW*w  + disInter1*dwdW[var] + ddisInter2dW*n3;
					ddis5dW[var] = dl234dW*dE   -l234*dEdW[var]    + ddisInter1dW*H  + disInter1*dHdW[var] +
					                                                 ddisInter2dW*Vn + disInter2*dVndW;
				}
			} else {
				// Flux term
				rhoVn = rhoR*VnR;

				duRdW[0] = -uR*rhoR_inv; duRdW[1] = rhoR_inv; duRdW[2] = 0.0;      duRdW[3] = 0.0;      duRdW[4] = 0.0;
				dvRdW[0] = -vR*rhoR_inv; dvRdW[1] = 0.0;      dvRdW[2] = rhoR_inv; dvRdW[3] = 0.0;      dvRdW[4] = 0.0;
				dwRdW[0] = -wR*rhoR_inv; dwRdW[1] = 0.0;      dwRdW[2] = 0.0;      dwRdW[3] = rhoR_inv; dwRdW[4] = 0.0;

				drhoRdW[0] = 1.0;     drhoRdW[1] = 0.0; drhoRdW[2] = 0.0; drhoRdW[3] = 0.0; drhoRdW[4] = 0.0;
				dEdW[0]    = 0.0;     dEdW[1]    = 0.0; dEdW[2]    = 0.0; dEdW[3]    = 0.0; dEdW[4]    = 1.0;
				dpdW[0]    = 0.5*V2R; dpdW[1]    = -uR; dpdW[2]    = -vR; dpdW[3]    = -wR; dpdW[4]    = 1.0;

				for (var = 0; var < Nvar; var++) {
					dpdW[var] *= GM1;

					dVnRdW[var] = n1*duRdW[var]+n2*dvRdW[var]+n3*dwRdW[var];
					drhoVndW    = drhoRdW[var]*VnR + rhoR*dVnRdW[var];

					dnF1dW[var] = drhoVndW;
					dnF2dW[var] = drhoVndW*uR + rhoVn*duRdW[var] + n1*dpdW[var];
					dnF3dW[var] = drhoVndW*vR + rhoVn*dvRdW[var] + n2*dpdW[var];
					dnF4dW[var] = drhoVndW*wR + rhoVn*dwRdW[var] + n3*dpdW[var];
					dnF5dW[var] = dVnRdW[var]*(ER+pR) + VnR*(dEdW[var]+dpdW[var]);
				}

				// Dissipation term
				mult      = rhoRr_inv/Den;
				drhodW[0] = 0.5*rhoL/rho;     drhodW[1] = 0.0;  drhodW[2] = 0.0;  drhodW[3] = 0.0;  drhodW[4] = 0.0;
				dudW[0]   = -0.5*(uR+u)*mult; dudW[1]   = mult; dudW[2]   = 0.0;  dudW[3]   = 0.0;  dudW[4]   = 0.0;
				dvdW[0]   = -0.5*(vR+v)*mult; dvdW[1]   = 0.0;  dvdW[2]   = mult; dvdW[3]   = 0.0;  dvdW[4]   = 0.0;
				dwdW[0]   = -0.5*(wR+w)*mult; dwdW[1]   = 0.0;  dwdW[2]   = 0.0;  dwdW[3]   = mult; dwdW[4]   = 0.0;

				dHdW[0] = -0.5*(HR+H-GM1*V2R)*mult;
				dHdW[1] = -GM1*uR*mult;
				dHdW[2] = -GM1*vR*mult;
				dHdW[3] = -GM1*wR*mult;
				dHdW[4] = GAMMA*mult;

				for (var = 0; var < Nvar; var++) {
					dcdW = 0.5*GM1/c*(dHdW[var]-(u*dudW[var]+v*dvdW[var]+w*dwdW[var]));
					dVndW = n1*dudW[var]+n2*dvdW[var]+n3*dwdW[var];

					if (case_l1)
						dl1dW = sign_l1*(dVndW-dcdW);
					else
						dl1dW = sign_l1*(-dcdW);

					if (case_l5)
						dl5dW = sign_l5*(dVndW+dcdW);
					else
						dl5dW = sign_l5*(dVnRdW[var]+dcdW);

					dl234dW = sign_l234*dVndW;

					dlc1dW = 0.5*(dl5dW+dl1dW) - dl234dW;
					dlc2dW = 0.5*(dl5dW-dl1dW);

					ddisInter1dW = (dlc1dW*dp+lc1*dpdW[var])/c2-(2.0*lc1*dp*dcdW)/(c*c2) +
					               (dlc2dW*rho*dVn+lc2*drhodW[var]*dVn+lc2*rho*dVnRdW[var])/c-(lc2*rho*dVn*dcdW)/c2;
					ddisInter2dW = dlc1dW*rho*dVn+lc1*drhodW[var]*dVn+lc1*rho*dVnRdW[var] +
					               (dlc2dW*dp+lc2*dpdW[var])/c-(lc2*dp*dcdW)/c2;

					drhoudW = drhoRdW[var]*uR+rhoR*duRdW[var];
					drhovdW = drhoRdW[var]*vR+rhoR*dvRdW[var];
					drhowdW = drhoRdW[var]*wR+rhoR*dwRdW[var];

					ddis1dW[var] = dl234dW*drho +l234*drhoRdW[var] + ddisInter1dW;
					ddis2dW[var] = dl234dW*drhou+l234*drhoudW      + ddisInter1dW*u  + disInter1*dudW[var] + ddisInter2dW*n1;
					ddis3dW[var] = dl234dW*drhov+l234*drhovdW      + ddisInter1dW*v  + disInter1*dvdW[var] + ddisInter2dW*n2;
					ddis4dW[var] = dl234dW*drhow+l234*drhowdW      + ddisInter1dW*w  + disInter1*dwdW[var] + ddisInter2dW*n3;
					ddis5dW[var] = dl234dW*dE   +l234*dEdW[var]    + ddisInter1dW*H  + disInter1*dHdW[var] +
					                                                 ddisInter2dW*Vn + disInter2*dVndW;
				}
			}

			InddnFdW = 0;
			for (var = 0; var < Nvar; var++) *dnFdW_ptr[InddnFdW++] = 0.5*(dnF1dW[var]-ddis1dW[var]);
			for (var = 0; var < Nvar; var++) *dnFdW_ptr[InddnFdW++] = 0.5*(dnF2dW[var]-ddis2dW[var]);
			for (var = 0; var < Nvar; var++) *dnFdW_ptr[InddnFdW++] = 0.5*(dnF3dW[var]-ddis3dW[var]);
			for (var = 0; var < Nvar; var++) *dnFdW_ptr[InddnFdW++] = 0.5*(dnF4dW[var]-ddis4dW[var]);
			for (var = 0; var < Nvar; var++) *dnFdW_ptr[InddnFdW++] = 0.5*(dnF5dW[var]-ddis5dW[var]);

			for (i = 0, iMax = Neq*Nvar; i < iMax; i++)
				dnFdW_ptr[i]++;
		}
	} else if (d == 2) {
		rhovL_ptr = &WL[NnTotal*2];
		rhovR_ptr = &WR[NnTotal*2];

		for (n = 0; n < NnTotal; n++) {
			n1 = *n_ptr++;
			n2 = *n_ptr++;

			// Inner VOLUME
			rhoL  = *rhoL_ptr++;
			rhouL = *rhouL_ptr++;
			rhovL = *rhovL_ptr++;
			EL    = *EL_ptr++;

			rhoL_inv  = 1.0/rhoL;
			rhoLr_inv = sqrt(rhoL_inv);

			uL = rhouL*rhoL_inv;
			vL = rhovL*rhoL_inv;

			V2L = uL*uL+vL*vL;
			VnL = n1*uL+n2*vL;

			pL  = GM1*(EL-0.5*rhoL*V2L);
			HL  = (EL+pL)*rhoL_inv;

			// Outer VOLUME
			rhoR  = *rhoR_ptr++;
			rhouR = *rhouR_ptr++;
			rhovR = *rhovR_ptr++;
			ER    = *ER_ptr++;

			rhoR_inv  = 1.0/rhoR;
			rhoRr_inv = sqrt(rhoR_inv);

			uR = rhouR*rhoR_inv;
			vR = rhovR*rhoR_inv;

			V2R = uR*uR+vR*vR;
			VnR = n1*uR+n2*vR;

			pR  = GM1*(ER-0.5*rhoR*V2R);
			HR  = (ER+pR)*rhoR_inv;

			// Roe-averaged states
			r = sqrt(rhoR/rhoL);
			rP1 = r+1;

			rho = r*rhoL;
			u   = (r*uR+uL)/rP1;
			v   = (r*vR+vL)/rP1;
			H   = (r*HR+HL)/rP1;
			Vn  = n1*u+n2*v;
			V2  = u*u+v*v;
			c2  = GM1*(H-0.5*V2);
			c   = sqrt(c2);

			Den  = sqrt(rhoL) + sqrt(rhoR);

			// Compute eigenvalues (with entropy fix)
			l1L = VnL-c;
			l1  = Vn-c;

			sign_l1 = 1.0;
			if (fabs(l1L) < fabs(l1)) {
				case_l1 = 0;
				if (l1L < 0.0)
					sign_l1 = -1.0;
				l1 = sign_l1*l1L;
			} else {
				case_l1 = 1;
				if (l1 < 0.0)
					sign_l1 = -1.0;
				l1 = sign_l1*l1;
			}

			l5R = VnR+c;
			l5  = Vn+c;

			sign_l5 = 1.0;
			if (fabs(l5R) > fabs(l5)) {
				case_l5 = 0;
				if (l5R < 0.0)
					sign_l5 = -1.0;
				l5 = sign_l5*l5R;
			} else {
				case_l5 = 1;
				if (l5 < 0.0)
					sign_l5 = -1.0;
				l5 = sign_l5*l5;
			}

			sign_l234 = 1.0;
			if (Vn < 0.0)
				sign_l234 = -1.0;
			l234 = sign_l234*(Vn);


			drho  = rhoR-rhoL;
			drhou = rhoR*uR-rhoL*uL;
			drhov = rhoR*vR-rhoL*vL;
			dE    = ER-EL;
			dp    = pR-pL;
			dVn   = VnR-VnL;

			lc1 = 0.5*(l5+l1) - l234;
			lc2 = 0.5*(l5-l1);

			disInter1 = lc1*dp/(c*c) + lc2*rho*dVn/c;
			disInter2 = lc1*rho*dVn  + lc2*dp/c;

			if (side == 'L') {
				// Flux term
				rhoVn = rhoL*VnL;

				duLdW[0] = -uL*rhoL_inv; duLdW[1] = rhoL_inv; duLdW[2] = 0.0;      duLdW[3] = 0.0;
				dvLdW[0] = -vL*rhoL_inv; dvLdW[1] = 0.0;      dvLdW[2] = rhoL_inv; dvLdW[3] = 0.0;

				drhoLdW[0] = 1.0;     drhoLdW[1] = 0.0; drhoLdW[2] = 0.0; drhoLdW[3] = 0.0;
				dEdW[0]    = 0.0;     dEdW[1]    = 0.0; dEdW[2]    = 0.0; dEdW[3]    = 1.0;
				dpdW[0]    = 0.5*V2L; dpdW[1]    = -uL; dpdW[2]    = -vL; dpdW[3]    = 1.0;

				for (var = 0; var < Nvar; var++) {
					dpdW[var] *= GM1;

					dVnLdW[var] = n1*duLdW[var]+n2*dvLdW[var];
					drhoVndW    = drhoLdW[var]*VnL + rhoL*dVnLdW[var];

					dnF1dW[var] = drhoVndW;
					dnF2dW[var] = drhoVndW*uL + rhoVn*duLdW[var] + n1*dpdW[var];
					dnF3dW[var] = drhoVndW*vL + rhoVn*dvLdW[var] + n2*dpdW[var];
					dnF5dW[var] = dVnLdW[var]*(EL+pL) + VnL*(dEdW[var]+dpdW[var]);
				}

				// Dissipation term
				mult      = rhoLr_inv/Den;
				drhodW[0] = 0.5*rhoR/rho;     drhodW[1] = 0.0;  drhodW[2] = 0.0;  drhodW[3] = 0.0;
				dudW[0]   = -0.5*(uL+u)*mult; dudW[1]   = mult; dudW[2]   = 0.0;  dudW[3]   = 0.0;
				dvdW[0]   = -0.5*(vL+v)*mult; dvdW[1]   = 0.0;  dvdW[2]   = mult; dvdW[3]   = 0.0;

				dHdW[0] = -0.5*(HL+H-GM1*V2L)*mult;
				dHdW[1] = -GM1*uL*mult;
				dHdW[2] = -GM1*vL*mult;
				dHdW[3] = GAMMA*mult;

				for (var = 0; var < Nvar; var++) {
					dcdW = 0.5*GM1/c*(dHdW[var]-(u*dudW[var]+v*dvdW[var]));
					dVndW = n1*dudW[var]+n2*dvdW[var];

					if (case_l1)
						dl1dW = sign_l1*(dVndW-dcdW);
					else
						dl1dW = sign_l1*(dVnLdW[var]-dcdW);

					if (case_l5)
						dl5dW = sign_l5*(dVndW+dcdW);
					else
						dl5dW = sign_l5*(dcdW);

					dl234dW = sign_l234*dVndW;

					dlc1dW = 0.5*(dl5dW+dl1dW) - dl234dW;
					dlc2dW = 0.5*(dl5dW-dl1dW);

					ddisInter1dW = (dlc1dW*dp-lc1*dpdW[var])/c2-(2.0*lc1*dp*dcdW)/(c*c2) +
					               (dlc2dW*rho*dVn+lc2*drhodW[var]*dVn-lc2*rho*dVnLdW[var])/c-(lc2*rho*dVn*dcdW)/c2;
					ddisInter2dW = dlc1dW*rho*dVn+lc1*drhodW[var]*dVn-lc1*rho*dVnLdW[var] +
					               (dlc2dW*dp-lc2*dpdW[var])/c-(lc2*dp*dcdW)/c2;

					drhoudW = drhoLdW[var]*uL+rhoL*duLdW[var];
					drhovdW = drhoLdW[var]*vL+rhoL*dvLdW[var];

					ddis1dW[var] = dl234dW*drho -l234*drhoLdW[var] + ddisInter1dW;
					ddis2dW[var] = dl234dW*drhou-l234*drhoudW      + ddisInter1dW*u  + disInter1*dudW[var] + ddisInter2dW*n1;
					ddis3dW[var] = dl234dW*drhov-l234*drhovdW      + ddisInter1dW*v  + disInter1*dvdW[var] + ddisInter2dW*n2;
					ddis5dW[var] = dl234dW*dE   -l234*dEdW[var]    + ddisInter1dW*H  + disInter1*dHdW[var] +
					                                                 ddisInter2dW*Vn + disInter2*dVndW;
				}
			} else {
				// Flux term
				rhoVn = rhoR*VnR;

				duRdW[0] = -uR*rhoR_inv; duRdW[1] = rhoR_inv; duRdW[2] = 0.0;      duRdW[3] = 0.0;
				dvRdW[0] = -vR*rhoR_inv; dvRdW[1] = 0.0;      dvRdW[2] = rhoR_inv; dvRdW[3] = 0.0;

				drhoRdW[0] = 1.0;     drhoRdW[1] = 0.0; drhoRdW[2] = 0.0; drhoRdW[3] = 0.0;
				dEdW[0]    = 0.0;     dEdW[1]    = 0.0; dEdW[2]    = 0.0; dEdW[3]    = 1.0;
				dpdW[0]    = 0.5*V2R; dpdW[1]    = -uR; dpdW[2]    = -vR; dpdW[3]    = 1.0;

				for (var = 0; var < Nvar; var++) {
					dpdW[var] *= GM1;

					dVnRdW[var] = n1*duRdW[var]+n2*dvRdW[var];
					drhoVndW    = drhoRdW[var]*VnR + rhoR*dVnRdW[var];

					dnF1dW[var] = drhoVndW;
					dnF2dW[var] = drhoVndW*uR + rhoVn*duRdW[var] + n1*dpdW[var];
					dnF3dW[var] = drhoVndW*vR + rhoVn*dvRdW[var] + n2*dpdW[var];
					dnF5dW[var] = dVnRdW[var]*(ER+pR) + VnR*(dEdW[var]+dpdW[var]);
				}

				// Dissipation term
				mult      = rhoRr_inv/Den;
				drhodW[0] = 0.5*rhoL/rho;     drhodW[1] = 0.0;  drhodW[2] = 0.0;  drhodW[3] = 0.0;
				dudW[0]   = -0.5*(uR+u)*mult; dudW[1]   = mult; dudW[2]   = 0.0;  dudW[3]   = 0.0;
				dvdW[0]   = -0.5*(vR+v)*mult; dvdW[1]   = 0.0;  dvdW[2]   = mult; dvdW[3]   = 0.0;

				dHdW[0] = -0.5*(HR+H-GM1*V2R)*mult;
				dHdW[1] = -GM1*uR*mult;
				dHdW[2] = -GM1*vR*mult;
				dHdW[3] = GAMMA*mult;

				for (var = 0; var < Nvar; var++) {
					dcdW = 0.5*GM1/c*(dHdW[var]-(u*dudW[var]+v*dvdW[var]));
					dVndW = n1*dudW[var]+n2*dvdW[var];

					if (case_l1)
						dl1dW = sign_l1*(dVndW-dcdW);
					else
						dl1dW = sign_l1*(-dcdW);

					if (case_l5)
						dl5dW = sign_l5*(dVndW+dcdW);
					else
						dl5dW = sign_l5*(dVnRdW[var]+dcdW);

					dl234dW = sign_l234*dVndW;

					dlc1dW = 0.5*(dl5dW+dl1dW) - dl234dW;
					dlc2dW = 0.5*(dl5dW-dl1dW);

					ddisInter1dW = (dlc1dW*dp+lc1*dpdW[var])/c2-(2.0*lc1*dp*dcdW)/(c*c2) +
					               (dlc2dW*rho*dVn+lc2*drhodW[var]*dVn+lc2*rho*dVnRdW[var])/c-(lc2*rho*dVn*dcdW)/c2;
					ddisInter2dW = dlc1dW*rho*dVn+lc1*drhodW[var]*dVn+lc1*rho*dVnRdW[var] +
					               (dlc2dW*dp+lc2*dpdW[var])/c-(lc2*dp*dcdW)/c2;

					drhoudW = drhoRdW[var]*uR+rhoR*duRdW[var];
					drhovdW = drhoRdW[var]*vR+rhoR*dvRdW[var];

					ddis1dW[var] = dl234dW*drho +l234*drhoRdW[var] + ddisInter1dW;
					ddis2dW[var] = dl234dW*drhou+l234*drhoudW      + ddisInter1dW*u  + disInter1*dudW[var] + ddisInter2dW*n1;
					ddis3dW[var] = dl234dW*drhov+l234*drhovdW      + ddisInter1dW*v  + disInter1*dvdW[var] + ddisInter2dW*n2;
					ddis5dW[var] = dl234dW*dE   +l234*dEdW[var]    + ddisInter1dW*H  + disInter1*dHdW[var] +
					                                                 ddisInter2dW*Vn + disInter2*dVndW;
				}
			}

			InddnFdW = 0;
			for (var = 0; var < Nvar; var++) *dnFdW_ptr[InddnFdW++] = 0.5*(dnF1dW[var]-ddis1dW[var]);
			for (var = 0; var < Nvar; var++) *dnFdW_ptr[InddnFdW++] = 0.5*(dnF2dW[var]-ddis2dW[var]);
			for (var = 0; var < Nvar; var++) *dnFdW_ptr[InddnFdW++] = 0.5*(dnF3dW[var]-ddis3dW[var]);
			for (var = 0; var < Nvar; var++) *dnFdW_ptr[InddnFdW++] = 0.5*(dnF5dW[var]-ddis5dW[var]);

			for (i = 0, iMax = Neq*Nvar; i < iMax; i++)
				dnFdW_ptr[i]++;
		}
	} else if (d == 1) {
		for (n = 0; n < NnTotal; n++) {
			n1 = *n_ptr++;

			// Inner VOLUME
			rhoL  = *rhoL_ptr++;
			rhouL = *rhouL_ptr++;
			EL    = *EL_ptr++;

			rhoL_inv  = 1.0/rhoL;
			rhoLr_inv = sqrt(rhoL_inv);

			uL = rhouL*rhoL_inv;

			V2L = uL*uL;
			VnL = n1*uL;

			pL  = GM1*(EL-0.5*rhoL*V2L);
			HL  = (EL+pL)*rhoL_inv;

			// Outer VOLUME
			rhoR  = *rhoR_ptr++;
			rhouR = *rhouR_ptr++;
			ER    = *ER_ptr++;

			rhoR_inv  = 1.0/rhoR;
			rhoRr_inv = sqrt(rhoR_inv);

			uR = rhouR*rhoR_inv;

			V2R = uR*uR;
			VnR = n1*uR;

			pR  = GM1*(ER-0.5*rhoR*V2R);
			HR  = (ER+pR)*rhoR_inv;

			// Roe-averaged states
			r = sqrt(rhoR/rhoL);
			rP1 = r+1;

			rho = r*rhoL;
			u   = (r*uR+uL)/rP1;
			H   = (r*HR+HL)/rP1;
			Vn  = n1*u;
			V2  = u*u;
			c2  = GM1*(H-0.5*V2);
			c   = sqrt(c2);

			Den  = sqrt(rhoL) + sqrt(rhoR);

			// Compute eigenvalues (with entropy fix)
			l1L = VnL-c;
			l1  = Vn-c;

			sign_l1 = 1.0;
			if (fabs(l1L) < fabs(l1)) {
				case_l1 = 0;
				if (l1L < 0.0)
					sign_l1 = -1.0;
				l1 = sign_l1*l1L;
			} else {
				case_l1 = 1;
				if (l1 < 0.0)
					sign_l1 = -1.0;
				l1 = sign_l1*l1;
			}

			l5R = VnR+c;
			l5  = Vn+c;

			sign_l5 = 1.0;
			if (fabs(l5R) > fabs(l5)) {
				case_l5 = 0;
				if (l5R < 0.0)
					sign_l5 = -1.0;
				l5 = sign_l5*l5R;
			} else {
				case_l5 = 1;
				if (l5 < 0.0)
					sign_l5 = -1.0;
				l5 = sign_l5*l5;
			}

			sign_l234 = 1.0;
			if (Vn < 0.0)
				sign_l234 = -1.0;
			l234 = sign_l234*(Vn);


			drho  = rhoR-rhoL;
			drhou = rhoR*uR-rhoL*uL;
			dE    = ER-EL;
			dp    = pR-pL;
			dVn   = VnR-VnL;

			lc1 = 0.5*(l5+l1) - l234;
			lc2 = 0.5*(l5-l1);

			disInter1 = lc1*dp/(c*c) + lc2*rho*dVn/c;
			disInter2 = lc1*rho*dVn  + lc2*dp/c;

			if (side == 'L') {
				// Flux term
				rhoVn = rhoL*VnL;

				duLdW[0] = -uL*rhoL_inv; duLdW[1] = rhoL_inv; duLdW[2] = 0.0;

				drhoLdW[0] = 1.0;     drhoLdW[1] = 0.0; drhoLdW[2] = 0.0;
				dEdW[0]    = 0.0;     dEdW[1]    = 0.0; dEdW[2]    = 1.0;
				dpdW[0]    = 0.5*V2L; dpdW[1]    = -uL; dpdW[2]    = 1.0;

				for (var = 0; var < Nvar; var++) {
					dpdW[var] *= GM1;

					dVnLdW[var] = n1*duLdW[var];
					drhoVndW    = drhoLdW[var]*VnL + rhoL*dVnLdW[var];

					dnF1dW[var] = drhoVndW;
					dnF2dW[var] = drhoVndW*uL + rhoVn*duLdW[var] + n1*dpdW[var];
					dnF5dW[var] = dVnLdW[var]*(EL+pL) + VnL*(dEdW[var]+dpdW[var]);
				}

				// Dissipation term
				mult      = rhoLr_inv/Den;
				drhodW[0] = 0.5*rhoR/rho;     drhodW[1] = 0.0;  drhodW[2] = 0.0;
				dudW[0]   = -0.5*(uL+u)*mult; dudW[1]   = mult; dudW[2]   = 0.0;

				dHdW[0] = -0.5*(HL+H-GM1*V2L)*mult;
				dHdW[1] = -GM1*uL*mult;
				dHdW[2] = GAMMA*mult;

				for (var = 0; var < Nvar; var++) {
					dcdW = 0.5*GM1/c*(dHdW[var]-(u*dudW[var]));
					dVndW = n1*dudW[var];

					if (case_l1)
						dl1dW = sign_l1*(dVndW-dcdW);
					else
						dl1dW = sign_l1*(dVnLdW[var]-dcdW);

					if (case_l5)
						dl5dW = sign_l5*(dVndW+dcdW);
					else
						dl5dW = sign_l5*(dcdW);

					dl234dW = sign_l234*dVndW;

					dlc1dW = 0.5*(dl5dW+dl1dW) - dl234dW;
					dlc2dW = 0.5*(dl5dW-dl1dW);

					ddisInter1dW = (dlc1dW*dp-lc1*dpdW[var])/c2-(2.0*lc1*dp*dcdW)/(c*c2) +
					               (dlc2dW*rho*dVn+lc2*drhodW[var]*dVn-lc2*rho*dVnLdW[var])/c-(lc2*rho*dVn*dcdW)/c2;
					ddisInter2dW = dlc1dW*rho*dVn+lc1*drhodW[var]*dVn-lc1*rho*dVnLdW[var] +
					               (dlc2dW*dp-lc2*dpdW[var])/c-(lc2*dp*dcdW)/c2;

					drhoudW = drhoLdW[var]*uL+rhoL*duLdW[var];

					ddis1dW[var] = dl234dW*drho -l234*drhoLdW[var] + ddisInter1dW;
					ddis2dW[var] = dl234dW*drhou-l234*drhoudW      + ddisInter1dW*u  + disInter1*dudW[var] + ddisInter2dW*n1;
					ddis5dW[var] = dl234dW*dE   -l234*dEdW[var]    + ddisInter1dW*H  + disInter1*dHdW[var] +
					                                                 ddisInter2dW*Vn + disInter2*dVndW;
				}
			} else {
				// Flux term
				rhoVn = rhoR*VnR;

				duRdW[0] = -uR*rhoR_inv; duRdW[1] = rhoR_inv; duRdW[2] = 0.0;

				drhoRdW[0] = 1.0;     drhoRdW[1] = 0.0; drhoRdW[2] = 0.0;
				dEdW[0]    = 0.0;     dEdW[1]    = 0.0; dEdW[2]    = 1.0;
				dpdW[0]    = 0.5*V2R; dpdW[1]    = -uR; dpdW[2]    = 1.0;

				for (var = 0; var < Nvar; var++) {
					dpdW[var] *= GM1;

					dVnRdW[var] = n1*duRdW[var];
					drhoVndW    = drhoRdW[var]*VnR + rhoR*dVnRdW[var];

					dnF1dW[var] = drhoVndW;
					dnF2dW[var] = drhoVndW*uR + rhoVn*duRdW[var] + n1*dpdW[var];
					dnF5dW[var] = dVnRdW[var]*(ER+pR) + VnR*(dEdW[var]+dpdW[var]);
				}

				// Dissipation term
				mult      = rhoRr_inv/Den;
				drhodW[0] = 0.5*rhoL/rho;     drhodW[1] = 0.0;  drhodW[2] = 0.0;
				dudW[0]   = -0.5*(uR+u)*mult; dudW[1]   = mult; dudW[2]   = 0.0;

				dHdW[0] = -0.5*(HR+H-GM1*V2R)*mult;
				dHdW[1] = -GM1*uR*mult;
				dHdW[2] = GAMMA*mult;

				for (var = 0; var < Nvar; var++) {
					dcdW = 0.5*GM1/c*(dHdW[var]-(u*dudW[var]));
					dVndW = n1*dudW[var];

					if (case_l1)
						dl1dW = sign_l1*(dVndW-dcdW);
					else
						dl1dW = sign_l1*(-dcdW);

					if (case_l5)
						dl5dW = sign_l5*(dVndW+dcdW);
					else
						dl5dW = sign_l5*(dVnRdW[var]+dcdW);

					dl234dW = sign_l234*dVndW;

					dlc1dW = 0.5*(dl5dW+dl1dW) - dl234dW;
					dlc2dW = 0.5*(dl5dW-dl1dW);

					ddisInter1dW = (dlc1dW*dp+lc1*dpdW[var])/c2-(2.0*lc1*dp*dcdW)/(c*c2) +
					               (dlc2dW*rho*dVn+lc2*drhodW[var]*dVn+lc2*rho*dVnRdW[var])/c-(lc2*rho*dVn*dcdW)/c2;
					ddisInter2dW = dlc1dW*rho*dVn+lc1*drhodW[var]*dVn+lc1*rho*dVnRdW[var] +
					               (dlc2dW*dp+lc2*dpdW[var])/c-(lc2*dp*dcdW)/c2;

					drhoudW = drhoRdW[var]*uR+rhoR*duRdW[var];

					ddis1dW[var] = dl234dW*drho +l234*drhoRdW[var] + ddisInter1dW;
					ddis2dW[var] = dl234dW*drhou+l234*drhoudW      + ddisInter1dW*u  + disInter1*dudW[var] + ddisInter2dW*n1;
					ddis5dW[var] = dl234dW*dE   +l234*dEdW[var]    + ddisInter1dW*H  + disInter1*dHdW[var] +
					                                                 ddisInter2dW*Vn + disInter2*dVndW;
				}
			}

			InddnFdW = 0;
			for (var = 0; var < Nvar; var++) *dnFdW_ptr[InddnFdW++] = 0.5*(dnF1dW[var]-ddis1dW[var]);
			for (var = 0; var < Nvar; var++) *dnFdW_ptr[InddnFdW++] = 0.5*(dnF2dW[var]-ddis2dW[var]);
			for (var = 0; var < Nvar; var++) *dnFdW_ptr[InddnFdW++] = 0.5*(dnF5dW[var]-ddis5dW[var]);

			for (i = 0, iMax = Neq*Nvar; i < iMax; i++)
				dnFdW_ptr[i]++;
		}
	}

}
