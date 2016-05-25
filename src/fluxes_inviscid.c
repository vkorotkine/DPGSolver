// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Compute inviscid fluxes from input W in conservative form.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void flux_inviscid (const unsigned int Nn, const unsigned int Nel, double *W, double *F, const unsigned int d,
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
