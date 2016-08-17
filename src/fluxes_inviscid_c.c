// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "fluxes_inviscid_c.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

#include "Parameters.h"
#include "Macros.h"
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

void flux_inviscid_c(const unsigned int Nn, const unsigned int Nel, double complex *W, double complex *F,
                     const unsigned int d, const unsigned int Neq)
{
	// Standard datatypes
	unsigned int   i, n, eq, dim, iMax, NnTotal, IndF;
	double complex *rho_ptr, *rhou_ptr, *rhov_ptr, *rhow_ptr, *E_ptr,
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

void flux_LF_c(const unsigned int Nn, const unsigned int Nel, double complex *WL, double complex *WR,
               double complex *nFluxNum, double *nL, const unsigned int d, const unsigned int Neq)
{
	// Standard datatypes
	unsigned int   i, iMax, jMax, NnTotal;
	double         *nx_ptr, *ny_ptr, *nz_ptr;
	double complex *rhoL, *uL, *vL, *wL, *pL, *UL, *rhoR, *uR, *vR, *wR, *pR, *UR, *FL, *FR, *maxV,
	               *maxV_ptr, *WL_ptr, *WR_ptr, *nFluxNum_ptr,
	               *FxL_ptr, *FyL_ptr, *FzL_ptr, *FxR_ptr, *FyR_ptr, *FzR_ptr;

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
