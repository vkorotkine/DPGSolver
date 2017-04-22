// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "fluxes_inviscid_c.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

#include "Parameters.h"
#include "Macros.h"
#include "Test.h"

#include "variable_functions_c.h"

#include "array_print.h"

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

void flux_inviscid_c(const unsigned int Nn, const unsigned int Nel, const double complex *const W,
                     double complex *const F, const unsigned int d, const unsigned int Neq)
{
	// Standard datatypes
	unsigned int   i, n, eq, dim, iMax, NnTotal, IndF;
	double complex rho, rhou, rhov, rhow, E, u, v, w, p, *F_ptr[DMAX*Neq];
	const double complex *rho_ptr, *rhou_ptr, *rhov_ptr, *rhow_ptr, *E_ptr;

	NnTotal = Nn*Nel;

	rho_ptr  = &W[NnTotal*0];
	rhou_ptr = &W[NnTotal*1];
	E_ptr    = &W[NnTotal*(d+1)];

	for (eq  = 0; eq  < Neq;  eq++)  {
	for (dim = 0; dim < d;    dim++) {
		F_ptr[eq*DMAX+dim] = &F[(eq*d+dim)*NnTotal];
	}}

	if (d == 3) {
		rhov_ptr = &W[NnTotal*2];
		rhow_ptr = &W[NnTotal*3];

		for (n = 0; n < NnTotal; n++) {
			rho  = *rho_ptr;
			rhou = *rhou_ptr;
			rhov = *rhov_ptr;
			rhow = *rhow_ptr;
			E    = *E_ptr;

			u   = rhou/rho;
			v   = rhov/rho;
			w   = rhow/rho;

			p = GM1*(E-0.5*rho*(u*u+v*v+w*w));

			IndF = 0;
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

			rho_ptr++; rhou_ptr++; rhov_ptr++; rhow_ptr++; E_ptr++;
			for (i = 0, iMax = Neq*DMAX; i < iMax; i++)
				F_ptr[i]++;
		}
	} else if (d == 2) {
		rhov_ptr = &W[NnTotal*2];

		for (n = 0; n < NnTotal; n++) {
			rho  = *rho_ptr;
			rhou = *rhou_ptr;
			rhov = *rhov_ptr;
			E    = *E_ptr;

			u   = rhou/rho;
			v   = rhov/rho;

			p = GM1*(E-0.5*rho*(u*u+v*v));

			IndF = 0;
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

			rho_ptr++; rhou_ptr++; rhov_ptr++; E_ptr++;
			for (i = 0, iMax = Neq*DMAX; i < iMax; i++)
				F_ptr[i]++;
		}
	} else if (d == 1) {
		for (n = 0; n < NnTotal; n++) {
			rho  = *rho_ptr;
			rhou = *rhou_ptr;
			E    = *E_ptr;

			u   = rhou/rho;

			p = GM1*(E-0.5*rho*(u*u));

			IndF = 0;
			// eq 1
			*F_ptr[IndF++] = rhou;
			IndF += 2;

			// eq 2
			*F_ptr[IndF++] = rhou*u + p;
			IndF += 2;

			// eq 3
			*F_ptr[IndF++] = (E+p)*u;
			IndF += 2;

			rho_ptr++; rhou_ptr++; E_ptr++;
			for (i = 0, iMax = Neq*DMAX; i < iMax; i++)
				F_ptr[i]++;
		}
	}
}

void flux_LF_c(const unsigned int Nn, const unsigned int Nel, const double complex *const WL,
               const double complex *const WR, double complex *const nFluxNum, const double *const nL,
               const unsigned int d, const unsigned int Neq)
{
	TestDB.EnteredInviscidFlux[0]++;

	// Standard datatypes
	unsigned int   i, iMax, jMax, NnTotal;
	double complex *rhoL, *uL, *vL, *wL, *pL, *UL, *rhoR, *uR, *vR, *wR, *pR, *UR, *FL, *FR, *maxV,
	               *maxV_ptr, *nFluxNum_ptr,
	               *FxL_ptr, *FyL_ptr, *FzL_ptr, *FxR_ptr, *FyR_ptr, *FzR_ptr;
	const double   *nx_ptr, *ny_ptr, *nz_ptr;
	const double complex *WL_ptr, *WR_ptr;

	NnTotal = Nn*Nel;

	UL   = malloc(NnTotal*Neq   * sizeof *UL);   // free
	UR   = malloc(NnTotal*Neq   * sizeof *UR);   // free
	FL   = malloc(NnTotal*Neq*d * sizeof *FL);   // free
	FR   = malloc(NnTotal*Neq*d * sizeof *FR);   // free
	maxV = malloc(NnTotal       * sizeof *maxV); // free

	convert_variables_c(WL,UL,d,d,Nn,Nel,'c','p');
	convert_variables_c(WR,UR,d,d,Nn,Nel,'c','p');

	flux_inviscid_c(Nn,Nel,WL,FL,d,Neq);
	flux_inviscid_c(Nn,Nel,WR,FR,d,Neq);

	rhoL = &UL[NnTotal*0];
	uL   = &UL[NnTotal*1];
	pL   = &UL[NnTotal*(d+1)];

	rhoR = &UR[NnTotal*0];
	uR   = &UR[NnTotal*1];
	pR   = &UR[NnTotal*(d+1)];

	if (d == 3) {
		vL = &UL[NnTotal*2];
		wL = &UL[NnTotal*3];

		vR = &UR[NnTotal*2];
		wR = &UR[NnTotal*3];

		// Compute wave speed
		maxV_ptr = maxV;
		for (iMax = NnTotal; iMax--; ) {

			if (creal(csqrt((*uL)*(*uL)+(*vL)*(*vL)+(*wL)*(*wL)) + csqrt(GAMMA*(*pL)/(*rhoL))) >
			    creal(csqrt((*uR)*(*uR)+(*vR)*(*vR)+(*wR)*(*wR)) + csqrt(GAMMA*(*pR)/(*rhoR)))) {
					*maxV_ptr = csqrt((*uL)*(*uL)+(*vL)*(*vL)+(*wL)*(*wL)) + csqrt(GAMMA*(*pL)/(*rhoL));
					TestDB.EnteredLF[0]++;
			} else {
				*maxV_ptr = csqrt((*uR)*(*uR)+(*vR)*(*vR)+(*wR)*(*wR)) + csqrt(GAMMA*(*pR)/(*rhoR));
				TestDB.EnteredLF[1]++;
			}

			maxV_ptr++;
			rhoL++; uL++; vL++; wL++; pL++;
			rhoR++; uR++; vR++; wR++; pR++;
		}

		// Compute n (dot) FluxNum
		nFluxNum_ptr = nFluxNum;
		WL_ptr       = WL;
		WR_ptr       = WR;
		for (i = 0; i < Neq; i++) {
			maxV_ptr = maxV;

			nx_ptr = &nL[0];
			ny_ptr = &nL[1];
			nz_ptr = &nL[2];

			FxL_ptr = &FL[NnTotal*(i*d+0)];
			FyL_ptr = &FL[NnTotal*(i*d+1)];
			FzL_ptr = &FL[NnTotal*(i*d+2)];
			FxR_ptr = &FR[NnTotal*(i*d+0)];
			FyR_ptr = &FR[NnTotal*(i*d+1)];
			FzR_ptr = &FR[NnTotal*(i*d+2)];

			for (jMax = NnTotal; jMax--; ) {
				*nFluxNum_ptr = 0.5 * ( (*nx_ptr)*((*FxL_ptr)+(*FxR_ptr))
				                       +(*ny_ptr)*((*FyL_ptr)+(*FyR_ptr))
				                       +(*nz_ptr)*((*FzL_ptr)+(*FzR_ptr)) + (*maxV_ptr)*((*WL_ptr)-(*WR_ptr)));

				nFluxNum_ptr++;
				nx_ptr += d; ny_ptr += d; nz_ptr += d;
				FxL_ptr++; FxR_ptr++;
				FyL_ptr++; FyR_ptr++;
				FzL_ptr++; FzR_ptr++;

				maxV_ptr++;
				WL_ptr++; WR_ptr++;
			}
		}
	} else if (d == 2) {
		vL = &UL[NnTotal*2];

		vR = &UR[NnTotal*2];

		// Compute wave speed
		maxV_ptr = maxV;
		for (iMax = NnTotal; iMax--; ) {
			if (creal(csqrt((*uL)*(*uL)+(*vL)*(*vL)) + csqrt(GAMMA*(*pL)/(*rhoL))) >
			    creal(csqrt((*uR)*(*uR)+(*vR)*(*vR)) + csqrt(GAMMA*(*pR)/(*rhoR)))) {
					*maxV_ptr = csqrt((*uL)*(*uL)+(*vL)*(*vL)) + csqrt(GAMMA*(*pL)/(*rhoL));
					TestDB.EnteredLF[0]++;
			} else {
				*maxV_ptr = csqrt((*uR)*(*uR)+(*vR)*(*vR)) + csqrt(GAMMA*(*pR)/(*rhoR));
				TestDB.EnteredLF[1]++;
			}

			maxV_ptr++;
			rhoL++; uL++; vL++; pL++;
			rhoR++; uR++; vR++; pR++;
		}

		// Compute n (dot) FluxNum
		nFluxNum_ptr = nFluxNum;
		WL_ptr       = WL;
		WR_ptr       = WR;
		for (i = 0; i < Neq; i++) {
			maxV_ptr = maxV;

			nx_ptr = &nL[0];
			ny_ptr = &nL[1];

			FxL_ptr = &FL[NnTotal*(i*d+0)];
			FyL_ptr = &FL[NnTotal*(i*d+1)];
			FxR_ptr = &FR[NnTotal*(i*d+0)];
			FyR_ptr = &FR[NnTotal*(i*d+1)];

			for (jMax = NnTotal; jMax--; ) {
				*nFluxNum_ptr = 0.5 * ( (*nx_ptr)*((*FxL_ptr)+(*FxR_ptr))
				                       +(*ny_ptr)*((*FyL_ptr)+(*FyR_ptr)) + (*maxV_ptr)*((*WL_ptr)-(*WR_ptr)));

				nFluxNum_ptr++;
				nx_ptr += d; ny_ptr += d;
				FxL_ptr++; FxR_ptr++;
				FyL_ptr++; FyR_ptr++;

				maxV_ptr++;
				WL_ptr++; WR_ptr++;
			}
		}
	} else if (d == 1) {
		// Compute wave speed
		maxV_ptr = maxV;
		for (iMax = NnTotal; iMax--; ) {
			if (creal(csqrt((*uL)*(*uL)) + csqrt(GAMMA*(*pL)/(*rhoL))) >
			    creal(csqrt((*uR)*(*uR)) + csqrt(GAMMA*(*pR)/(*rhoR)))) {
					*maxV_ptr = csqrt((*uL)*(*uL)) + csqrt(GAMMA*(*pL)/(*rhoL));
					TestDB.EnteredLF[0]++;
			} else {
				*maxV_ptr = csqrt((*uR)*(*uR)) + csqrt(GAMMA*(*pR)/(*rhoR));
				TestDB.EnteredLF[1]++;
			}

			maxV_ptr++;
			rhoL++; uL++; pL++;
			rhoR++; uR++; pR++;
		}

		// Compute n (dot) FluxNum
		nFluxNum_ptr = nFluxNum;
		WL_ptr       = WL;
		WR_ptr       = WR;
		for (i = 0; i < Neq; i++) {
			maxV_ptr = maxV;

			nx_ptr = &nL[0];

			FxL_ptr = &FL[NnTotal*(i*d+0)];
			FxR_ptr = &FR[NnTotal*(i*d+0)];

			for (jMax = NnTotal; jMax--; ) {
				*nFluxNum_ptr = 0.5 * ( (*nx_ptr)*((*FxL_ptr)+(*FxR_ptr)) + (*maxV_ptr)*((*WL_ptr)-(*WR_ptr)));

				nFluxNum_ptr++;
				nx_ptr += d;
				FxL_ptr++; FxR_ptr++;

				maxV_ptr++;
				WL_ptr++; WR_ptr++;
			}
		}
	}

	free(UL);
	free(UR);
	free(FL);
	free(FR);
	free(maxV);
}

void flux_Roe_c(const unsigned int Nn, const unsigned int Nel, const double complex *const WL,
                const double complex *const WR, double complex *const nFluxNum, const double *const nL,
                const unsigned int d, const unsigned int Neq)
{
	TestDB.EnteredInviscidFlux[1]++;

	// Standard datatypes
	unsigned int   iMax, NnTotal;
	double         sign_l1, sign_l234, sign_l5;
	double complex r, rP1, rho, u, v, w, H, Vn, V2, c, l1, l234, l5, l1L, l5R,
	               VnL, rhoVnL, VnR, rhoVnR, pLR, drho, drhou, drhov, drhow, dE, dp, dVn, lc1, lc2, disInter1, disInter2,
	               rhoL, uL, vL, wL, pL, EL, rhoR, uR, vR, wR, pR, ER,
	               *nFluxNum_ptr1, *nFluxNum_ptr2, *nFluxNum_ptr3, *nFluxNum_ptr4, *nFluxNum_ptr5,
	               dis1, dis2, dis3, dis4, dis5, nF1, nF2, nF3, nF4, nF5;
	const double   *nx, *ny, *nz;
	const double complex *W1L, *W2L, *W3L, *W4L, *W5L, *W1R, *W2R, *W3R, *W4R, *W5R;

	// silence
	iMax = Neq;


	NnTotal = Nn*Nel;

	if (d == 3) {
		nx = &nL[0];
		ny = &nL[1];
		nz = &nL[2];

		W1L = &WL[NnTotal*0];
		W2L = &WL[NnTotal*1];
		W3L = &WL[NnTotal*2];
		W4L = &WL[NnTotal*3];
		W5L = &WL[NnTotal*(d+1)];

		W1R = &WR[NnTotal*0];
		W2R = &WR[NnTotal*1];
		W3R = &WR[NnTotal*2];
		W4R = &WR[NnTotal*3];
		W5R = &WR[NnTotal*(d+1)];

		nFluxNum_ptr1 = &nFluxNum[NnTotal*0];
		nFluxNum_ptr2 = &nFluxNum[NnTotal*1];
		nFluxNum_ptr3 = &nFluxNum[NnTotal*2];
		nFluxNum_ptr4 = &nFluxNum[NnTotal*3];
		nFluxNum_ptr5 = &nFluxNum[NnTotal*(d+1)];

		for (iMax = NnTotal; iMax--; ) {
			// Initialize left and right states at the current node
			rhoL = (*W1L++);
			uL   = (*W2L++)/rhoL;
			vL   = (*W3L++)/rhoL;
			wL   = (*W4L++)/rhoL;
			EL   = (*W5L++);
			pL   = GM1*(EL-0.5*rhoL*(uL*uL+vL*vL+wL*wL));

			rhoR = (*W1R++);
			uR   = (*W2R++)/rhoR;
			vR   = (*W3R++)/rhoR;
			wR   = (*W4R++)/rhoR;
			ER   = (*W5R++);
			pR   = GM1*(ER-0.5*rhoR*(uR*uR+vR*vR+wR*wR));

			// Compute Roe-averaged states
			r = csqrt(rhoR/rhoL);
			rP1 = r+1;

			rho = r*rhoL;
			u   = (r*uR+uL)/rP1;
			v   = (r*vR+vL)/rP1;
			w   = (r*wR+wL)/rP1;
			H   = (r*(ER+pR)/rhoR+(EL+pL)/rhoL)/rP1;
			Vn  = (*nx)*u+(*ny)*v+(*nz)*w;
			V2  = u*u+v*v+w*w;
			c   = csqrt(GM1*(H-0.5*V2));

			// Compute eigenvalues (with entropy fix)
			VnL = (*nx)*uL+(*ny)*vL+(*nz)*wL;
			VnR = (*nx)*uR+(*ny)*vR+(*nz)*wR;

			l1L = VnL-c;
			l1  = Vn-c;

			sign_l1 = 1.0;
			if (cabs(l1L) < cabs(l1)) {
				TestDB.EnteredRoe[0]++;
				if (creal(l1L) < 0.0)
					sign_l1 = -1.0;
				l1 = sign_l1*l1L;
			} else {
				TestDB.EnteredRoe[1]++;
				if (creal(l1) < 0.0)
					sign_l1 = -1.0;
				l1 = sign_l1*l1;
			}

			l5R = VnR+c;
			l5  = Vn+c;

			sign_l5 = 1.0;
			if (cabs(l5R) > cabs(l5)) {
				TestDB.EnteredRoe[2]++;
				if (creal(l5R) < 0.0)
					sign_l5 = -1.0;
				l5 = sign_l5*l5R;
			} else {
				TestDB.EnteredRoe[3]++;
				if (creal(l5) < 0.0)
					sign_l5 = -1.0;
				l5 = sign_l5*l5;
			}

			sign_l234 = 1.0;
			if (creal(Vn) < 0.0)
				sign_l234 = -1.0;
			l234 = sign_l234*Vn;


			// Compute combined eigenvalues, eigenvectors and linearized wave strengths
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

			dis1 = l234*drho  + disInter1;
			dis2 = l234*drhou + disInter1*u + disInter2*(*nx);
			dis3 = l234*drhov + disInter1*v + disInter2*(*ny);
			dis4 = l234*drhow + disInter1*w + disInter2*(*nz);
			dis5 = l234*dE    + disInter1*H + disInter2*(Vn);

			// Compute contribution of normal flux components (multiplied by 0.5 below)
			rhoVnL = rhoL*VnL;
			rhoVnR = rhoR*VnR;
			pLR    = pL + pR;

			nF1 = rhoVnL      + rhoVnR;
			nF2 = rhoVnL*uL   + rhoVnR*uR  + (*nx)*pLR;
			nF3 = rhoVnL*vL   + rhoVnR*vR  + (*ny)*pLR;
			nF4 = rhoVnL*wL   + rhoVnR*wR  + (*nz)*pLR;
			nF5 = VnL*(EL+pL) + VnR*(ER+pR);

			// Assemble components
			*nFluxNum_ptr1++ = 0.5*(nF1 - dis1);
			*nFluxNum_ptr2++ = 0.5*(nF2 - dis2);
			*nFluxNum_ptr3++ = 0.5*(nF3 - dis3);
			*nFluxNum_ptr4++ = 0.5*(nF4 - dis4);
			*nFluxNum_ptr5++ = 0.5*(nF5 - dis5);

			nx += d; ny += d; nz += d;
		}
	} else if (d == 2) {
		nx = &nL[0];
		ny = &nL[1];

		W1L = &WL[NnTotal*0];
		W2L = &WL[NnTotal*1];
		W3L = &WL[NnTotal*2];
		W5L = &WL[NnTotal*(d+1)];

		W1R = &WR[NnTotal*0];
		W2R = &WR[NnTotal*1];
		W3R = &WR[NnTotal*2];
		W5R = &WR[NnTotal*(d+1)];

		nFluxNum_ptr1 = &nFluxNum[NnTotal*0];
		nFluxNum_ptr2 = &nFluxNum[NnTotal*1];
		nFluxNum_ptr3 = &nFluxNum[NnTotal*2];
		nFluxNum_ptr5 = &nFluxNum[NnTotal*(d+1)];

		for (iMax = NnTotal; iMax--; ) {
			// Initialize left and right states at the current node
			rhoL = (*W1L++);
			uL   = (*W2L++)/rhoL;
			vL   = (*W3L++)/rhoL;
			EL   = (*W5L++);
			pL   = GM1*(EL-0.5*rhoL*(uL*uL+vL*vL));

			rhoR = (*W1R++);
			uR   = (*W2R++)/rhoR;
			vR   = (*W3R++)/rhoR;
			ER   = (*W5R++);
			pR   = GM1*(ER-0.5*rhoR*(uR*uR+vR*vR));

			// Compute Roe-averaged states
			r = csqrt(rhoR/rhoL);
			rP1 = r+1;

			rho = r*rhoL;
			u   = (r*uR+uL)/rP1;
			v   = (r*vR+vL)/rP1;
			H   = (r*(ER+pR)/rhoR+(EL+pL)/rhoL)/rP1;
			Vn  = (*nx)*u+(*ny)*v;
			V2  = u*u+v*v;
			c   = csqrt(GM1*(H-0.5*V2));

			// Compute eigenvalues (with entropy fix)
			VnL = (*nx)*uL+(*ny)*vL;
			VnR = (*nx)*uR+(*ny)*vR;

			l1L = VnL-c;
			l1  = Vn-c;

			sign_l1 = 1.0;
			if (cabs(l1L) < cabs(l1)) {
				TestDB.EnteredRoe[0]++;
				if (creal(l1L) < 0.0)
					sign_l1 = -1.0;
				l1 = sign_l1*l1L;
			} else {
				TestDB.EnteredRoe[1]++;
				if (creal(l1) < 0.0)
					sign_l1 = -1.0;
				l1 = sign_l1*l1;
			}

			l5R = VnR+c;
			l5  = Vn+c;

			sign_l5 = 1.0;
			if (cabs(l5R) > cabs(l5)) {
				TestDB.EnteredRoe[2]++;
				if (creal(l5R) < 0.0)
					sign_l5 = -1.0;
				l5 = sign_l5*l5R;
			} else {
				TestDB.EnteredRoe[3]++;
				if (creal(l5) < 0.0)
					sign_l5 = -1.0;
				l5 = sign_l5*l5;
			}

			sign_l234 = 1.0;
			if (creal(Vn) < 0.0)
				sign_l234 = -1.0;
			l234 = sign_l234*Vn;


			// Compute combined eigenvalues, eigenvectors and linearized wave strengths
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

			dis1 = l234*drho  + disInter1;
			dis2 = l234*drhou + disInter1*u + disInter2*(*nx);
			dis3 = l234*drhov + disInter1*v + disInter2*(*ny);
			dis5 = l234*dE    + disInter1*H + disInter2*(Vn);

			// Compute contribution of normal flux components (multiplied by 0.5 below)
			rhoVnL = rhoL*VnL;
			rhoVnR = rhoR*VnR;
			pLR    = pL + pR;

			nF1 = rhoVnL      + rhoVnR;
			nF2 = rhoVnL*uL   + rhoVnR*uR  + (*nx)*pLR;
			nF3 = rhoVnL*vL   + rhoVnR*vR  + (*ny)*pLR;
			nF5 = VnL*(EL+pL) + VnR*(ER+pR);

			// Assemble components
			*nFluxNum_ptr1++ = 0.5*(nF1 - dis1);
			*nFluxNum_ptr2++ = 0.5*(nF2 - dis2);
			*nFluxNum_ptr3++ = 0.5*(nF3 - dis3);
			*nFluxNum_ptr5++ = 0.5*(nF5 - dis5);

			nx += d; ny += d;
		}
	} else if (d == 1) {
		nx = &nL[0];

		W1L = &WL[NnTotal*0];
		W2L = &WL[NnTotal*1];
		W5L = &WL[NnTotal*(d+1)];

		W1R = &WR[NnTotal*0];
		W2R = &WR[NnTotal*1];
		W5R = &WR[NnTotal*(d+1)];

		nFluxNum_ptr1 = &nFluxNum[NnTotal*0];
		nFluxNum_ptr2 = &nFluxNum[NnTotal*1];
		nFluxNum_ptr5 = &nFluxNum[NnTotal*(d+1)];

		for (iMax = NnTotal; iMax--; ) {
			// Initialize left and right states at the current node
			rhoL = (*W1L++);
			uL   = (*W2L++)/rhoL;
			EL   = (*W5L++);
			pL   = GM1*(EL-0.5*rhoL*(uL*uL));

			rhoR = (*W1R++);
			uR   = (*W2R++)/rhoR;
			ER   = (*W5R++);
			pR   = GM1*(ER-0.5*rhoR*(uR*uR));

			// Compute Roe-averaged states
			r = csqrt(rhoR/rhoL);
			rP1 = r+1;

			rho = r*rhoL;
			u   = (r*uR+uL)/rP1;
			H   = (r*(ER+pR)/rhoR+(EL+pL)/rhoL)/rP1;
			Vn  = (*nx)*u;
			V2  = u*u;
			c   = csqrt(GM1*(H-0.5*V2));

			// Compute eigenvalues (with entropy fix)
			VnL = (*nx)*uL;
			VnR = (*nx)*uR;

			l1L = VnL-c;
			l1  = Vn-c;

			sign_l1 = 1.0;
			if (cabs(l1L) < cabs(l1)) {
				TestDB.EnteredRoe[0]++;
				if (creal(l1L) < 0.0)
					sign_l1 = -1.0;
				l1 = sign_l1*l1L;
			} else {
				TestDB.EnteredRoe[1]++;
				if (creal(l1) < 0.0)
					sign_l1 = -1.0;
				l1 = sign_l1*l1;
			}

			l5R = VnR+c;
			l5  = Vn+c;

			sign_l5 = 1.0;
			if (cabs(l5R) > cabs(l5)) {
				TestDB.EnteredRoe[2]++;
				if (creal(l5R) < 0.0)
					sign_l5 = -1.0;
				l5 = sign_l5*l5R;
			} else {
				TestDB.EnteredRoe[3]++;
				if (creal(l5) < 0.0)
					sign_l5 = -1.0;
				l5 = sign_l5*l5;
			}

			sign_l234 = 1.0;
			if (creal(Vn) < 0.0)
				sign_l234 = -1.0;
			l234 = sign_l234*Vn;


			// Compute combined eigenvalues, eigenvectors and linearized wave strengths
			drho  = rhoR-rhoL;
			drhou = rhoR*uR-rhoL*uL;
			dE    = ER-EL;
			dp    = pR-pL;
			dVn   = VnR-VnL;

			lc1 = 0.5*(l5+l1) - l234;
			lc2 = 0.5*(l5-l1);

			disInter1 = lc1*dp/(c*c) + lc2*rho*dVn/c;
			disInter2 = lc1*rho*dVn  + lc2*dp/c;

			dis1 = l234*drho  + disInter1;
			dis2 = l234*drhou + disInter1*u + disInter2*(*nx);
			dis5 = l234*dE    + disInter1*H + disInter2*(Vn);

			// Compute contribution of normal flux components (multiplied by 0.5 below)
			rhoVnL = rhoL*VnL;
			rhoVnR = rhoR*VnR;
			pLR    = pL + pR;

			nF1 = rhoVnL      + rhoVnR;
			nF2 = rhoVnL*uL   + rhoVnR*uR  + (*nx)*pLR;
			nF5 = VnL*(EL+pL) + VnR*(ER+pR);

			// Assemble components
			*nFluxNum_ptr1++ = 0.5*(nF1 - dis1);
			*nFluxNum_ptr2++ = 0.5*(nF2 - dis2);
			*nFluxNum_ptr5++ = 0.5*(nF5 - dis5);

			nx += d;
		}
	}
}
