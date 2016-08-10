// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "variable_functions.h"

#include <stdlib.h>
#include <stdio.h>

#include "Parameters.h"

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

void convert_variables(double *VarIn, double *VarOut, const unsigned int dIn, const unsigned int dOut,
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
	double       *rho, *u, *v, *w, *p, *rhou, *rhov, *rhow, *E, zeros[Nn], rhoV2[Nn], *U, *W;

	// silence
	u = NULL; v = NULL; w = NULL;
	rhov = NULL; rhow = NULL;

	NnTotal = Nn*Nel;

	for (n = 0; n < NnTotal; n++)
		zeros[n] = 0.0;

	switch(TypeIn) {
	case 'p':
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
	case 'c':
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
			printf("Error: Unsupported TypeIn/Out combination in convert_variables.\n"), exit(1);
			break;
		}
		break;
	default:
		printf("Error: Unsupported TypeIn in convert_variables.\n"), exit(1);
		break;
	}
}
