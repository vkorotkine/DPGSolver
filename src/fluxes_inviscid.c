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
 *		It is assumed that inputs: W, n and outputs: F, nFluxNum are vectorized (i.e. the memory ordering is by equation
 *		and not by element).
 *		Try using BLAS calls for dot products and check if there is a speed-up. (ToBeDeleted)
 *
 *	Notation:
 *
 *	References:
 */

void flux_inviscid(const unsigned int Nn, const unsigned int Nel, double *W, double *F, const unsigned int d,
                   const unsigned int Neq)
{
	// Standard datatypes
	unsigned int i, j, iMax, NnTotal;
	double *rho, *rhou, *rhov, *rhow, *E, *u, *v, *w, *p, *Fptr[15];

	NnTotal = Nn*Nel;

	rho  = &W[NnTotal*0];
	rhou = &W[NnTotal*1];
	E    = &W[NnTotal*(d+1)];

	for (i = 0; i < Neq; i++) {
	for (j = 0; j < d;   j++) {
		Fptr[i*3+j] = &F[(i*d+j)*NnTotal];
	}}

	if (d == 3) {
		rhov = &W[NnTotal*2];
		rhow = &W[NnTotal*3];

		u = malloc(Nn*Nel * sizeof *u); // free
		v = malloc(Nn*Nel * sizeof *v); // free
		w = malloc(Nn*Nel * sizeof *w); // free
		p = malloc(Nn*Nel * sizeof *p); // free

		for (i = 0; i < NnTotal; i++) {
			u[i] = rhou[i]/rho[i];
			v[i] = rhov[i]/rho[i];
			w[i] = rhow[i]/rho[i];
			p[i] = GM1*(E[i]-0.5*rho[i]*(u[i]*u[i]+v[i]*v[i]+w[i]*w[i]));
		}

		for (i = 0; i < NnTotal; i++) {
			// eq 1
			Fptr[0][i] = rhou[i];
			Fptr[1][i] = rhov[i];
			Fptr[2][i] = rhow[i];

			// eq 2
			Fptr[3][i] = rhou[i]*u[i] + p[i];
			Fptr[4][i] = rhou[i]*v[i];
			Fptr[5][i] = rhou[i]*w[i];

			// eq 3
			Fptr[6][i] = rhov[i]*u[i];
			Fptr[7][i] = rhov[i]*v[i] + p[i];
			Fptr[8][i] = rhov[i]*w[i];

			// eq 4
			Fptr[9][i]  = rhow[i]*u[i];
			Fptr[10][i] = rhow[i]*v[i];
			Fptr[11][i] = rhow[i]*w[i] + p[i];

			// eq 5
			Fptr[12][i] = (E[i]+p[i])*u[i];
			Fptr[13][i] = (E[i]+p[i])*v[i];
			Fptr[14][i] = (E[i]+p[i])*w[i];
		}
		free(u);
		free(v);
		free(w);
		free(p);
	} else if (d == 2) {
		rhov = &W[NnTotal*2];

		u = malloc(Nn*Nel * sizeof *u); // free
		v = malloc(Nn*Nel * sizeof *v); // free
		p = malloc(Nn*Nel * sizeof *p); // free

		for (i = 0; i < NnTotal; i++) {
			u[i] = rhou[i]/rho[i];
			v[i] = rhov[i]/rho[i];
			p[i] = GM1*(E[i]-0.5*rho[i]*(u[i]*u[i]+v[i]*v[i]));
		}

		for (i = 0; i < NnTotal; i++) {
			// eq 1
			Fptr[0][i] = rhou[i];
			Fptr[1][i] = rhov[i];

			// eq 2
			Fptr[3][i] = rhou[i]*u[i] + p[i];
			Fptr[4][i] = rhou[i]*v[i];

			// eq 3
			Fptr[6][i] = rhov[i]*u[i];
			Fptr[7][i] = rhov[i]*v[i] + p[i];

			// eq 4
			Fptr[9][i]  = (E[i]+p[i])*u[i];
			Fptr[10][i] = (E[i]+p[i])*v[i];
		}
		free(u);
		free(v);
		free(p);
	} else if (d == 1) {
		u = malloc(Nn*Nel * sizeof *u); // free
		p = malloc(Nn*Nel * sizeof *p); // free

		for (i = 0, iMax = Nn*Nel; i < iMax; i++) {
			u[i] = rhou[i]/rho[i];
			p[i] = GM1*(E[i]-0.5*rho[i]*(u[i]*u[i]));
		}

		for (i = 0; i < NnTotal; i++) {
			// eq 1
			Fptr[0][i] = rhou[i];

			// eq 2
			Fptr[3][i] = rhou[i]*u[i] + p[i];

			// eq 3
			Fptr[6][i] = (E[i]+p[i])*u[i];
		}
		free(u);
		free(p);
	}
}

void flux_LF(const unsigned int Nn, const unsigned int Nel, double *WIn, double *WOut, double *nFluxNum,
             double *nIn, const unsigned int d, const unsigned int Neq)
{
	// Standard datatypes
	unsigned int i, iMax, jMax, NnTotal;
	double       *rhoIn, *uIn, *vIn, *wIn, *pIn, *UIn, *rhoOut, *uOut, *vOut, *wOut, *pOut, *UOut, *FIn, *FOut, *maxV,
	             *maxVptr, *WInptr, *WOutptr, *nxptr, *nyptr, *nzptr, *nFluxNumptr,
	             *FxInptr, *FyInptr, *FzInptr, *FxOutptr, *FyOutptr, *FzOutptr;

	NnTotal = Nn*Nel;

	UIn  = malloc(NnTotal*Neq   * sizeof *UIn);  // free
	UOut = malloc(NnTotal*Neq   * sizeof *UOut); // free

	FIn  = malloc(NnTotal*Neq*d * sizeof *FIn);  // free
	FOut = malloc(NnTotal*Neq*d * sizeof *FOut); // free

	maxV = malloc(NnTotal       * sizeof *maxV); // free

	convert_variables(WIn,UIn,d,d,Nn,Nel,'c','p');
	convert_variables(WOut,UOut,d,d,Nn,Nel,'c','p');

	flux_inviscid(Nn,Nel,WIn,FIn,d,Neq);
	flux_inviscid(Nn,Nel,WOut,FOut,d,Neq);


	rhoIn = &UIn[NnTotal*0];
	uIn   = &UIn[NnTotal*1];
	pIn   = &UIn[NnTotal*(d+1)];

	rhoOut = &UOut[NnTotal*0];
	uOut   = &UOut[NnTotal*1];
	pOut   = &UOut[NnTotal*(d+1)];

	if (d == 3) {
		vIn = &UIn[NnTotal*2];
		wIn = &UIn[NnTotal*3];

		vOut = &UOut[NnTotal*2];
		wOut = &UOut[NnTotal*3];

		// Compute wave speed
		maxVptr = maxV;
		for (iMax = NnTotal; iMax--; ) {
			*maxVptr = max(sqrt((*uIn)*(*uIn)+(*vIn)*(*vIn)+(*wIn)*(*wIn))       + sqrt(GAMMA*(*pIn)/(*rhoIn)),
			               sqrt((*uOut)*(*uOut)+(*vOut)*(*vOut)+(*wOut)*(*wOut)) + sqrt(GAMMA*(*pOut)/(*rhoOut)));

			maxVptr++;
			rhoIn++; uIn++; vIn++; wIn++; pIn++;
			rhoOut++; uOut++; vOut++; wOut++; pOut++;
		}

		// Compute n (dot) FluxNum
		nFluxNumptr = nFluxNum;
		WInptr      = WIn;
		WOutptr     = WOut;
		for (i = 0; i < Neq; i++) {
			maxVptr  = maxV;

			nxptr = &nIn[0];
			nyptr = &nIn[1];
			nzptr = &nIn[2];

			FxInptr  = &FIn[NnTotal*(i*d+0)];
			FyInptr  = &FIn[NnTotal*(i*d+1)];
			FzInptr  = &FIn[NnTotal*(i*d+2)];
			FxOutptr = &FOut[NnTotal*(i*d+0)];
			FyOutptr = &FOut[NnTotal*(i*d+1)];
			FzOutptr = &FOut[NnTotal*(i*d+2)];

			for (jMax = NnTotal; jMax--; ) {
if (i == 0) {
	printf("% .4e % .4e % .4e % .4e\n",*nxptr,*FxInptr,*FxOutptr,(*nxptr)*((*FxInptr)+(*FxOutptr)));
}
				*nFluxNumptr = 0.5*(  (*nxptr)*((*FxInptr)+(*FxOutptr))
				                    + (*nyptr)*((*FyInptr)+(*FyOutptr))
				                    + (*nzptr)*((*FzInptr)+(*FzOutptr)) + (*maxVptr)*((*WInptr)-(*WOutptr)));

				nFluxNumptr++;
				nxptr += d; nyptr += d; nzptr += d;
				FxInptr++; FxOutptr++;
				FyInptr++; FyOutptr++;
				FzInptr++; FzOutptr++;

				maxVptr++;
			}
		}
	}

	free(UIn);
	free(UOut);
	free(FIn);
	free(FOut);
	free(maxV);
}
