// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "variable_functions.h"

#include <stdlib.h>
#include <stdio.h>

#include "Parameters.h"
#include "Macros.h"

/*
 *	Purpose:
 *		Provide functions relating to supported variables.
 *
 *	Comments:
 *		It is assumed that inputs: W, n and outputs: F, nFluxNum are vectorized (i.e. the memory ordering is by equation
 *		and not by element).
 *
 *	Notation:
 *
 *	References:
 */

void convert_variables(const double *const VarIn, double *const VarOut, const unsigned int dIn, const unsigned int dOut,
                       const unsigned int Nn, const unsigned int Nel, const char TypeIn, const char TypeOut)
{
	/*
	 *	Purpose:
	 *		Convert between different variable types.
	 *
	 *	Comments:
	 *		Currently, only primitive and conservative variables are supported.
	 *		Potentially add in entropy variables later. (ToBeDeleted)
	 *
	 *		Variables:
	 *			Primitive:    U = [rho u v w P]
	 *			Conservative: U = [rho*[1 u v w] E]; E = P/(GAMMA-1) + 0.5*rho*V^2
	 */

	// Standard datatypes
	unsigned int n, NnTotal, varInMax = dIn + 1, varOutMax = dOut+1;
	double       *U, *W;

	NnTotal = Nn*Nel;

	double zeros[NnTotal], rhoV2[NnTotal];

	for (n = 0; n < NnTotal; n++)
		zeros[n] = 0.0;

	switch(TypeIn) {
	case 'p': {
		const double *rho, *u, *v, *w, *p;

		// silence
		v = w = NULL;

		rho = &VarIn[NnTotal*0];
		u   = &VarIn[NnTotal*1];
		p   = &VarIn[NnTotal*varInMax];
		if (dIn == 3) {
			v = &VarIn[NnTotal*2];
			w = &VarIn[NnTotal*3];
		} else if (dIn == 2) {
			v = &VarIn[NnTotal*2];
			w = zeros;
		} else if (dIn == 1) {
			v = zeros;
			w = zeros;
		}
		switch(TypeOut) {
		case 'c':
			W = VarOut;

			for (n = 0; n < NnTotal; n++)
				rhoV2[n] = rho[n]*(u[n]*u[n]+v[n]*v[n]+w[n]*w[n]);

			for (n = 0; n < NnTotal; n++) {
				W[0*NnTotal+n]         = rho[n];
				W[1*NnTotal+n]         = rho[n]*u[n];
				W[varOutMax*NnTotal+n] = p[n]/GM1 + 0.5*rhoV2[n];
			}

			if (dOut == 3) {
				for (n = 0; n < NnTotal; n++) {
					W[2*NnTotal+n] = rho[n]*v[n];
					W[3*NnTotal+n] = rho[n]*w[n];
				}
			} else if (dOut == 2) {
				for (n = 0; n < NnTotal; n++) {
					W[2*NnTotal+n] = rho[n]*v[n];
				}
			}
			break;
		default:
			printf("Error: Unsupported TypeIn/Out combination in convert_variables.\n"), exit(1);
			break;
		}
		break;
	} case 'c': {
		const double *rho, *rhou, *rhov, *rhow, *E;

		// silence
		rhov = rhow = NULL;

		rho  = &VarIn[NnTotal*0];
		rhou = &VarIn[NnTotal*1];
		E    = &VarIn[NnTotal*varInMax];
		if (dIn == 3) {
			rhov = &VarIn[NnTotal*2];
			rhow = &VarIn[NnTotal*3];
		} else if (dIn == 2) {
			rhov = &VarIn[NnTotal*2];
			rhow = zeros;
		} else if (dIn == 1) {
			rhov = zeros;
			rhow = zeros;
		}
		switch(TypeOut) {
		case 'p':
			U = VarOut;

			for (n = 0; n < NnTotal; n++)
				rhoV2[n] = (rhou[n]*rhou[n]+rhov[n]*rhov[n]+rhow[n]*rhow[n])/rho[n];

			for (n = 0; n < NnTotal; n++) {
				U[0*NnTotal+n]         = rho[n];
				U[1*NnTotal+n]         = rhou[n]/rho[n];
				U[varOutMax*NnTotal+n] = GM1*(E[n]-0.5*rhoV2[n]);
			}

			if (dOut == 3) {
				for (n = 0; n < NnTotal; n++) {
					U[2*NnTotal+n] = rhov[n]/rho[n];
					U[3*NnTotal+n] = rhow[n]/rho[n];
				}
			} else if (dOut == 2) {
				for (n = 0; n < NnTotal; n++) {
					U[2*NnTotal+n] = rhov[n]/rho[n];
				}
			}
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
		}
		break;
	} default: {
		EXIT_UNSUPPORTED;
		break;
	}}
}

void compute_pressure(double *VarIn, double *p, const unsigned int d, const unsigned int Nn, const unsigned int Nel,
                      const char TypeIn)
{
	if (d < 2)
		printf("Error: Unsupported.\n"), EXIT_MSG;

	unsigned int NnTotal = Nn*Nel;

	switch(TypeIn) {
	case 'p': {
		printf("Error: No need to compute pressure for the pressure primitive variables.\n"), EXIT_MSG;
		break;
	} case 'c': {
		double *rho, *rhou, *rhov, *E;
		rho  = &VarIn[NnTotal*0];
		rhou = &VarIn[NnTotal*1];
		rhov = &VarIn[NnTotal*2];
		E    = &VarIn[NnTotal*(d+1)];

		if (d == 2) {
			for (size_t n = 0; n < NnTotal; n++) {
				double rhoV2 = rhou[n]*rhou[n] + rhov[n]*rhov[n];
				p[n] = GM1*(E[n]-0.5*rhoV2/rho[n]);
			}
		} else if (d == 3) {
			double *rhow = &VarIn[NnTotal*3];
			for (size_t n = 0; n < NnTotal; n++) {
				double rhoV2 = rhou[n]*rhou[n] + rhov[n]*rhov[n] + rhow[n]*rhow[n];
				p[n] = GM1*(E[n]-0.5*rhoV2/rho[n]);
			}
		}
		break;
	} default: {
		printf("Error: Unsupported.\n"), EXIT_MSG;
		break;
	}}
}
