// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "fluxes_inviscid_c.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

#include "Parameters.h"

/*
 *	Purpose:
 *		Identical functions to fluxes_inviscid using complex variables used for complex step verification.
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
