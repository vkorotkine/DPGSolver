// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Compute inviscid fluxes from input W in conservative form.
 *
 *	Comments:
 *		The memory allocations in flux_inviscid can be removed as in jacobian_flux_inviscid (ToBeDeleted).
 *		It is assumed that inputs: W, n and outputs: F, nFluxNum are vectorized (i.e. the memory ordering is by equation
 *		and not by element).
 *		Numerical flux functions were intentionally written with if statements for each dimension in order to avoid
 *		redundant calculations if d < 3.
 *		Try using BLAS calls for dot products and check if there is a speed-up. (ToBeDeleted)
 *		The Roe flux function is very difficult to parallelize due to the large amount of required memory.
 *		There is some additional discussion in Obward(2015) regarding incorrect dissipation of the ROE scheme at low
 *		mach numbers. The proposed method to correct this requires evaluation of a shock indicator on all FACETs
 *		adjacent to the inner VOLUME of the FACET on which the numerical flux is being computed (which seems quite
 *		expensive to do). In the same paper, there is a reference to Thornber(2008) which motivates a local normal
 *		velocity modification to correct the overly dissipative numerical flux. However, this is done only in the
 *		context of a finite volume scheme with MUSCL reconstruction. Perhaps look into whether this problem is important
 *		for high-order methods and add a fix similar to Thornber(2008). (ToBeModified)
 *
 *	Notation:
 *
 *	References:
 *		flux_LF : Based off of Hesthaven's NodalDG function.(https://github.com/tcew/nodal-dg/.../CFD2D/EulerLF2D)
 *		flux_ROE: Toro(2009)-Riemann_Solvers_and_Numerical_Methods_for_Fluid_Dynamics (Ch. 11.3)
 *		        : Obwald(2015)-L2Roe:_A_Low_Dissipation_Version_of_Roe's_Approximate_Riemann_Solver_for_Low_Mach_Numbers
 *				: Thornber(2008)-An_improved_reconstruction_method_for_compressible_flows_with_low_Mach_number_features
 */

void flux_inviscid(const unsigned int Nn, const unsigned int Nel, double *W, double *F, const unsigned int d,
                   const unsigned int Neq)
{
	// Standard datatypes
	unsigned int i, n, eq, dim, iMax, NnTotal, IndF;
	double       *rho_ptr, *rhou_ptr, *rhov_ptr, *rhow_ptr, *E_ptr,
	             rho, rhou, rhov, rhow, E, u, v, w, p, *F_ptr[DMAX*Neq];

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

void flux_LF(const unsigned int Nn, const unsigned int Nel, double *WL, double *WR, double *nFluxNum,
             double *nL, const unsigned int d, const unsigned int Neq)
{
	// Standard datatypes
	unsigned int i, iMax, jMax, NnTotal;
	double       *rhoL, *uL, *vL, *wL, *pL, *UL, *rhoR, *uR, *vR, *wR, *pR, *UR, *FL, *FR, *maxV,
	             *maxV_ptr, *WL_ptr, *WR_ptr, *nx_ptr, *ny_ptr, *nz_ptr, *nFluxNum_ptr,
	             *FxL_ptr, *FyL_ptr, *FzL_ptr, *FxR_ptr, *FyR_ptr, *FzR_ptr;

	NnTotal = Nn*Nel;

	UL   = malloc(NnTotal*Neq   * sizeof *UL);   // free
	UR   = malloc(NnTotal*Neq   * sizeof *UR);   // free
	FL   = malloc(NnTotal*Neq*d * sizeof *FL);   // free
	FR   = malloc(NnTotal*Neq*d * sizeof *FR);   // free
	maxV = malloc(NnTotal       * sizeof *maxV); // free

	convert_variables(WL,UL,d,d,Nn,Nel,'c','p');
	convert_variables(WR,UR,d,d,Nn,Nel,'c','p');

	flux_inviscid(Nn,Nel,WL,FL,d,Neq);
	flux_inviscid(Nn,Nel,WR,FR,d,Neq);

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
			*maxV_ptr = max(sqrt((*uL)*(*uL)+(*vL)*(*vL)+(*wL)*(*wL)) + sqrt(GAMMA*(*pL)/(*rhoL)),
			                sqrt((*uR)*(*uR)+(*vR)*(*vR)+(*wR)*(*wR)) + sqrt(GAMMA*(*pR)/(*rhoR)));

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
			*maxV_ptr = max(sqrt((*uL)*(*uL)+(*vL)*(*vL)) + sqrt(GAMMA*(*pL)/(*rhoL)),
			                sqrt((*uR)*(*uR)+(*vR)*(*vR)) + sqrt(GAMMA*(*pR)/(*rhoR)));

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
			*maxV_ptr = max(sqrt((*uL)*(*uL)) + sqrt(GAMMA*(*pL)/(*rhoL)),
			                sqrt((*uR)*(*uR)) + sqrt(GAMMA*(*pR)/(*rhoR)));

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

void flux_ROE(const unsigned int Nn, const unsigned int Nel, double *WL, double *WR, double *nFluxNum,
              double *nL, const unsigned int d, const unsigned int Neq)
{
	/*
	 *	Comments:
	 *		This is the Roe-Pike version of the scheme which is different from the original Roe scheme in that the wave
	 *		numbers are linearized for faster computation.
	 */
	// Standard datatypes
	unsigned int iMax, NnTotal;
	double       eps, r, rP1, rho, u, v, w, H, Vn, V2, c, l1, l234, l5,
	             VnL, rhoVnL, VnR, rhoVnR, pLR, drho, drhou, drhov, drhow, dE, dp, dVn, lc1, lc2, disInter1, disInter2,
	             *W1L, *W2L, *W3L, *W4L, *W5L, *W1R, *W2R, *W3R, *W4R, *W5R,
	             rhoL, uL, vL, wL, pL, EL, cL, rhoR, uR, vR, wR, pR, ER, cR, dl1, dl5,
	             *nx, *ny, *nz,
	             *nFluxNum_ptr1, *nFluxNum_ptr2, *nFluxNum_ptr3, *nFluxNum_ptr4, *nFluxNum_ptr5,
	             dis1, dis2, dis3, dis4, dis5, nF1, nF2, nF3, nF4, nF5;

	eps = 1e-4;

	// silence
	iMax = Neq;
	r    = eps;

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
			r = sqrt(rhoR/rhoL);
			rP1 = r+1;

			rho = r*rhoL;
			u   = (r*uR+uL)/rP1;
			v   = (r*vR+vL)/rP1;
			w   = (r*wR+wL)/rP1;
			H   = (r*(ER+pR)/rhoR+(EL+pL)/rhoL)/rP1;
			Vn  = (*nx)*u+(*ny)*v+(*nz)*w;
			V2  = u*u+v*v+w*w;
			c   = sqrt(GM1*(H-0.5*V2));

			// Compute eigenvalues (with entropy fix if required)
			cL  = sqrt(GAMMA*pL/rhoL);
			VnL = (*nx)*uL+(*ny)*vL+(*nz)*wL;

			cR  = sqrt(GAMMA*pR/rhoR);
			VnR = (*nx)*uR+(*ny)*vR+(*nz)*wR;

			l1   = fabs(Vn-c);
			l234 = fabs(Vn);
			l5   = fabs(Vn+c);

			dl1 = max(fabs(VnR-cR)-fabs(VnL-cL),0.0);
			dl5 = max(fabs(VnR+cR)-fabs(VnL+cL),0.0);

			if (l1 < 2*dl1)
				l1 = (l1*l1)/(4*dl1)+dl1;
			if (l5 < 2*dl5)
				l5 = (l5*l5)/(4*dl5)+dl5;

/*
			if (l1 < eps)
				l1 = (l1*l1+eps*eps)/(2*eps);
			if (l5 < eps)
				l5 = (l5*l5+eps*eps)/(2*eps);
*/

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
			r = sqrt(rhoR/rhoL);
			rP1 = r+1;

			rho = r*rhoL;
			u   = (r*uR+uL)/rP1;
			v   = (r*vR+vL)/rP1;
			H   = (r*(ER+pR)/rhoR+(EL+pL)/rhoL)/rP1;
			Vn  = (*nx)*u+(*ny)*v;
			V2  = u*u+v*v;
			c   = sqrt(GM1*(H-0.5*V2));

			// Compute eigenvalues (with entropy fix if required)
			cL  = sqrt(GAMMA*pL/rhoL);
			VnL = (*nx)*uL+(*ny)*vL;

			cR  = sqrt(GAMMA*pR/rhoR);
			VnR = (*nx)*uR+(*ny)*vR;

			l1   = fabs(Vn-c);
			l234 = fabs(Vn);
			l5   = fabs(Vn+c);

			dl1 = max(fabs(VnR-cR)-fabs(VnL-cL),0.0);
			dl5 = max(fabs(VnR+cR)-fabs(VnL+cL),0.0);

			if (l1 < 2*dl1)
				l1 = (l1*l1)/(4*dl1)+dl1;
			if (l5 < 2*dl5)
				l5 = (l5*l5)/(4*dl5)+dl5;

			// Compute combined eigenvalues, eigenvectors and linearized wave strengths
			VnL = (*nx)*uL+(*ny)*vL;
			VnR = (*nx)*uR+(*ny)*vR;

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
			r = sqrt(rhoR/rhoL);
			rP1 = r+1;

			rho = r*rhoL;
			u   = (r*uR+uL)/rP1;
			H   = (r*(ER+pR)/rhoR+(EL+pL)/rhoL)/rP1;
			Vn  = (*nx)*u;
			V2  = u*u;
			c   = sqrt(GM1*(H-0.5*V2));

			// Compute eigenvalues (with entropy fix if required)
			cL  = sqrt(GAMMA*pL/rhoL);
			VnL = (*nx)*uL;

			cR  = sqrt(GAMMA*pR/rhoR);
			VnR = (*nx)*uR;

			l1   = fabs(Vn-c);
			l234 = fabs(Vn);
			l5   = fabs(Vn+c);

			dl1 = max(fabs(VnR-cR)-fabs(VnL-cL),0.0);
			dl5 = max(fabs(VnR+cR)-fabs(VnL+cL),0.0);

			if (l1 < 2*dl1)
				l1 = (l1*l1)/(4*dl1)+dl1;
			if (l5 < 2*dl5)
				l5 = (l5*l5)/(4*dl5)+dl5;

			// Compute combined eigenvalues, eigenvectors and linearized wave strengths
			VnL = (*nx)*uL;
			VnR = (*nx)*uR;

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
