// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "jacobian_fluxes_inviscid.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Parameters.h"
#include "Macros.h"

#include "fluxes_structs.h"
#include "matrix_structs.h"
#include "fluxes_inviscid.h"
#include "solver_Advection_functions.h"

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

static void jacobian_flux_Euler     (struct S_FLUX *const FLUXDATA);
static void jacobian_flux_Advection (struct S_FLUX *const FLUXDATA);

void jacobian_flux_inviscid(struct S_FLUX *const FLUXDATA)
{
	switch(FLUXDATA->PDE_index) {
		case PDE_ADVECTION:    jacobian_flux_Advection(FLUXDATA); break;
		case PDE_EULER:        // fallthrough
		case PDE_NAVIERSTOKES: jacobian_flux_Euler(FLUXDATA);     break;
		default:               EXIT_UNSUPPORTED;                  break;
	}
}

static void jacobian_flux_LF     (struct S_NUMERICALFLUX *const NUMFLUXDATA);
static void jacobian_flux_Roe    (struct S_NUMERICALFLUX *const NUMFLUXDATA);
static void jacobian_flux_upwind (struct S_NUMERICALFLUX *const NUMFLUXDATA);

void jacobian_flux_num_inviscid(struct S_NUMERICALFLUX *const NUMFLUXDATA)
{
	switch(NUMFLUXDATA->NumFluxInviscid_index) {
		case FLUX_LF:     jacobian_flux_LF(NUMFLUXDATA);     break;
		case FLUX_ROE:    jacobian_flux_Roe(NUMFLUXDATA);    break;
		case FLUX_UPWIND: jacobian_flux_upwind(NUMFLUXDATA); break;
		default:          EXIT_UNSUPPORTED;                  break;
	}
}

void jacobian_flux_inviscid_M (struct S_FLUX_M *const FLUXDATA_M)
{
	struct S_FLUX FLUXDATA;

	FLUXDATA.d = FLUXDATA_M->d;
	FLUXDATA.Nn  = FLUXDATA_M->W->extents[0];
	FLUXDATA.Nel = 1;

	FLUXDATA.W    = FLUXDATA_M->W->data;
	FLUXDATA.F    = FLUXDATA_M->F->data;
	FLUXDATA.dFdW = FLUXDATA_M->dFdW->data;

	FLUXDATA.XYZ = FLUXDATA_M->XYZ->data;

	jacobian_flux_inviscid(&FLUXDATA);
}

static void jacobian_flux_Euler(struct S_FLUX *const FLUXDATA)
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
		for (eq = 0; eq < Neq; eq++) {
		for (dim = 0; dim < d; dim++) {
			F_ptr[eq*DMAX+dim] = &F[(eq*d+dim)*NnTotal];
		}}
	}

	for (eq = 0; eq < Neq; eq++) {
	for (var = 0; var < Nvar; var++) {
	for (dim = 0; dim < d; dim++) {
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

static void jacobian_flux_Advection(struct S_FLUX *const FLUXDATA)
{
	/*
	 *	Comments:
	 *		Implicitly assumed that Neq = Nvar = 1.
	 */

	unsigned int const d       = FLUXDATA->d,
	                   Nn      = FLUXDATA->Nn,
	                   Nel     = FLUXDATA->Nel,
	                   NnTotal = Nn*Nel;

	double const *const W    = FLUXDATA->W;
	double       *const F    = FLUXDATA->F,
	             *const dFdW = FLUXDATA->dFdW;

	double *F_ptr[d];
	if (F != NULL) {
		for (size_t dim = 0; dim < d; dim++)
			F_ptr[dim] = &F[dim*NnTotal];
	}

	double *dFdW_ptr[d];
	for (size_t dim = 0; dim < d; dim++)
		dFdW_ptr[dim] = &dFdW[dim*NnTotal];

	double const *const b = compute_b_Advection(NnTotal,FLUXDATA->XYZ); // free

	for (size_t n = 0; n < NnTotal; n++) {
		if (F != NULL) {
			for (size_t dim = 0; dim < d; dim++) {
				*F_ptr[dim] = b[dim*NnTotal+n]*W[n];
				F_ptr[dim]++;
			}
		}

		for (size_t dim = 0; dim < d; dim++) {
			*dFdW_ptr[dim] = b[dim*NnTotal+n];
			dFdW_ptr[dim]++;
		}
	}
	free((double *) b);
}


static void jacobian_flux_LF(struct S_NUMERICALFLUX *const NUMFLUXDATA)
{
	unsigned int const d       = NUMFLUXDATA->d,
	                   Neq     = d+2,
	                   Nvar    = d+2,
	                   Nn      = NUMFLUXDATA->Nn,
	                   Nel     = NUMFLUXDATA->Nel,
	                   NnTotal = Nn*Nel;

	double const *const nL = NUMFLUXDATA->nL;

	double const *const WL = NUMFLUXDATA->WL,
	             *const WR = NUMFLUXDATA->WR;

	double *const nFluxNum = NUMFLUXDATA->nFluxNum,
	       *const dnFdWL   = NUMFLUXDATA->dnFluxNumdWL,
	       *const dnFdWR   = NUMFLUXDATA->dnFluxNumdWR;

	double const *rhoL_ptr  = &WL[NnTotal*0],
	             *rhouL_ptr = &WL[NnTotal*1],
	             *EL_ptr    = &WL[NnTotal*(d+1)],

	             *rhoR_ptr  = &WR[NnTotal*0],
	             *rhouR_ptr = &WR[NnTotal*1],
	             *ER_ptr    = &WR[NnTotal*(d+1)];

	double const *n_ptr = nL;

	double *nF_ptr[Neq];
	if (nFluxNum != NULL) {
		for (size_t eq = 0; eq < Neq; eq++)
			nF_ptr[eq] = &nFluxNum[eq*NnTotal];
	}

	double *dnFdWL_ptr[Neq*Neq];
	for (size_t eq = 0; eq < Neq; eq++) {
	for (size_t var = 0; var < Nvar; var++) {
		dnFdWL_ptr[eq*Nvar+var] = &dnFdWL[(eq*Nvar+var)*NnTotal];
	}}

	double *dnFdWR_ptr[Neq*Neq];
	if (dnFdWR != NULL) {
		for (size_t eq = 0; eq < Neq; eq++) {
		for (size_t var = 0; var < Nvar; var++) {
			dnFdWR_ptr[eq*Nvar+var] = &dnFdWR[(eq*Nvar+var)*NnTotal];
		}}
	}

	struct S_FLUX *const FLUXDATA = malloc(sizeof *FLUXDATA); // free
	FLUXDATA->d   = d;
	FLUXDATA->Nn  = 1;
	FLUXDATA->Nel = 1;
	FLUXDATA->F   = NULL;

	if (d == 3) {
		double const *rhovL_ptr = &WL[NnTotal*2],
		             *rhowL_ptr = &WL[NnTotal*3],

		             *rhovR_ptr = &WR[NnTotal*2],
		             *rhowR_ptr = &WR[NnTotal*3];

		for (size_t n = 0; n < NnTotal; n++) {
			// Left VOLUME
			double const rhoL  = *rhoL_ptr++,
			             rhouL = *rhouL_ptr++,
			             rhovL = *rhovL_ptr++,
			             rhowL = *rhowL_ptr++,
			             EL    = *EL_ptr++,

			             rhoL_inv = 1.0/rhoL,
			             uL = rhouL*rhoL_inv,
			             vL = rhovL*rhoL_inv,
			             wL = rhowL*rhoL_inv,

			             V2L = uL*uL+vL*vL+wL*wL,
			             VL  = sqrt(V2L),

			             pL  = GM1*(EL-0.5*rhoL*V2L),
			             cL  = sqrt(GAMMA*pL/rhoL);

			// Right VOLUME
			double const rhoR  = *rhoR_ptr++,
			             rhouR = *rhouR_ptr++,
			             rhovR = *rhovR_ptr++,
			             rhowR = *rhowR_ptr++,
			             ER    = *ER_ptr++,

			             rhoR_inv = 1.0/rhoR,
			             uR = rhouR*rhoR_inv,
			             vR = rhovR*rhoR_inv,
			             wR = rhowR*rhoR_inv,

			             V2R = uR*uR+vR*vR+wR*wR,
			             VR  = sqrt(V2R),

			             pR  = GM1*(ER-0.5*rhoR*V2R),
			             cR  = sqrt(GAMMA*pR/rhoR);


			double const maxlL = VL+cL,
			             maxlR = VR+cR;

			char sideMaxV;
			double maxV;
			if (maxlL > maxlR) {
				sideMaxV = 'L';
				maxV = maxlL;
			} else {
				sideMaxV = 'R';
				maxV = maxlR;
			}

			double const n1 = *n_ptr++,
			             n2 = *n_ptr++,
			             n3 = *n_ptr++;

			if (nFluxNum != NULL) {
				double WLn[] = {rhoL, rhouL, rhovL, rhowL, EL},
				       FLn[Neq*d];
				FLUXDATA->W = WLn;
				FLUXDATA->F = FLn;
				flux_Euler(FLUXDATA);

				double const *FL1_ptr = FLn,
				             *FL2_ptr = FL1_ptr+1,
				             *FL3_ptr = FL2_ptr+1;

				double WRn[] = {rhoR, rhouR, rhovR, rhowR, ER},
				       FRn[Neq*d];
				FLUXDATA->W = WRn;
				FLUXDATA->F = FRn;
				flux_Euler(FLUXDATA);

				double const *FR1_ptr = FRn,
				             *FR2_ptr = FR1_ptr+1,
				             *FR3_ptr = FR2_ptr+1;

				size_t IndnF = 0;
				for (size_t eq = 0; eq < Neq; eq++) {
					*nF_ptr[IndnF++]++ = 0.5*( n1*((*FL1_ptr)+(*FR1_ptr))+n2*((*FL2_ptr)+(*FR2_ptr))
					                          +n3*((*FL3_ptr)+(*FR3_ptr)) + maxV*(WLn[eq]-WRn[eq]));

					FL1_ptr += d;
					FL2_ptr += d;
					FL3_ptr += d;
					FR1_ptr += d;
					FR2_ptr += d;
					FR3_ptr += d;
				}
			}

			double dFdWn[Nvar*Nvar*d];

			// Flux term
			double const WLn[] = {rhoL, rhouL, rhovL, rhowL, EL};

			FLUXDATA->W    = WLn;
			FLUXDATA->dFdW = dFdWn;
			jacobian_flux_Euler(FLUXDATA);

			double const *dF1dW_ptr = dFdWn,
			             *dF2dW_ptr = dF1dW_ptr+1,
			             *dF3dW_ptr = dF2dW_ptr+1;

			size_t InddnFdW = 0;
			for (size_t eq = 0; eq < Neq; eq++) {
			for (size_t var = 0; var < Nvar; var++) {
				*dnFdWL_ptr[InddnFdW++] = 0.5*(n1*(*dF1dW_ptr)+n2*(*dF2dW_ptr)+n3*(*dF3dW_ptr));

				dF1dW_ptr += d;
				dF2dW_ptr += d;
				dF3dW_ptr += d;
			}}

			// Dissipation term
			InddnFdW = 0;
			for (size_t eq = 0; eq < Neq; eq++) {
			for (size_t var = 0; var < Nvar; var++) {
				if (var == eq)
					*dnFdWL_ptr[InddnFdW] += 0.5*maxV;
				InddnFdW++;
			}}
			if (sideMaxV == 'L') {
				double const WRn[]    = { rhoR,         rhouR,     rhovR,     rhowR,    ER},
				             dudW[]   = {-uL*rhoL_inv,  rhoL_inv,  0.0,       0.0,      0.0},
				             dvdW[]   = {-vL*rhoL_inv,  0.0,       rhoL_inv,  0.0,      0.0},
				             dwdW[]   = {-wL*rhoL_inv,  0.0,       0.0,       rhoL_inv, 0.0},
				             drhodW[] = { 1.0,          0.0,       0.0,       0.0,      0.0},
				             dpdW[]   = { GM1*0.5*V2L, -GM1*uL,   -GM1*vL,   -GM1*wL,   GM1};

				double dmaxVdW[Nvar];
				for (size_t var = 0; var < Nvar; var++) {
					double const dVdW = 1.0/VL*(uL*dudW[var]+vL*dvdW[var]+wL*dwdW[var]),
					             dcdW = GAMMA/(2.0*cL*rhoL*rhoL)*(dpdW[var]*rhoL-pL*drhodW[var]);
					dmaxVdW[var] = dVdW+dcdW;
				}

				InddnFdW = 0;
				for (size_t eq = 0; eq < Neq; eq++) {
				for (size_t var = 0; var < Nvar; var++) {
					*dnFdWL_ptr[InddnFdW++] += 0.5*dmaxVdW[var]*(WLn[eq]-WRn[eq]);
				}}
			}

			for (size_t i = 0, iMax = Neq*Nvar; i < iMax; i++)
				dnFdWL_ptr[i]++;

			if (dnFdWR != NULL) {
				// Flux term
				double const WRn[] = {rhoR, rhouR, rhovR, rhowR, ER};

				FLUXDATA->W    = WRn;
				FLUXDATA->dFdW = dFdWn;
				jacobian_flux_Euler(FLUXDATA);

				double const *dF1dW_ptr = dFdWn,
				             *dF2dW_ptr = dF1dW_ptr+1,
				             *dF3dW_ptr = dF2dW_ptr+1;

				size_t InddnFdW = 0;
				for (size_t eq = 0; eq < Neq; eq++) {
				for (size_t var = 0; var < Nvar; var++) {
					*dnFdWR_ptr[InddnFdW++] = 0.5*(n1*(*dF1dW_ptr)+n2*(*dF2dW_ptr)+n3*(*dF3dW_ptr));

					dF1dW_ptr += d;
					dF2dW_ptr += d;
					dF3dW_ptr += d;
				}}

				// Dissipation term
				InddnFdW = 0;
				for (size_t eq = 0; eq < Neq; eq++) {
				for (size_t var = 0; var < Nvar; var++) {
					if (var == eq)
						*dnFdWR_ptr[InddnFdW] -= 0.5*maxV;
					InddnFdW++;
				}}
				if (sideMaxV == 'R') {
					double const WLn[]    = { rhoL,         rhouL,     rhovL,     rhowL,    EL},
					             dudW[]   = {-uR*rhoR_inv,  rhoR_inv,  0.0,       0.0,      0.0},
					             dvdW[]   = {-vR*rhoR_inv,  0.0,       rhoR_inv,  0.0,      0.0},
					             dwdW[]   = {-wR*rhoR_inv,  0.0,       0.0,       rhoR_inv, 0.0},
					             drhodW[] = { 1.0,          0.0,       0.0,       0.0,      0.0},
					             dpdW[]   = { GM1*0.5*V2R, -GM1*uR,   -GM1*vR,   -GM1*wR,   GM1};

					double dmaxVdW[Nvar];
					for (size_t var = 0; var < Nvar; var++) {
						double const dVdW = 1.0/VR*(uR*dudW[var]+vR*dvdW[var]+wR*dwdW[var]),
						             dcdW = GAMMA/(2.0*cR*rhoR*rhoR)*(dpdW[var]*rhoR-pR*drhodW[var]);
						dmaxVdW[var] = dVdW+dcdW;
					}

					InddnFdW = 0;
					for (size_t eq = 0; eq < Neq; eq++) {
					for (size_t var = 0; var < Nvar; var++) {
						*dnFdWR_ptr[InddnFdW++] += 0.5*dmaxVdW[var]*(WLn[eq]-WRn[eq]);
					}}
				}

				for (size_t i = 0, iMax = Neq*Nvar; i < iMax; i++)
					dnFdWR_ptr[i]++;
			}
		}
	} else if (d == 2) {
		double const *rhovL_ptr = &WL[NnTotal*2],

		             *rhovR_ptr = &WR[NnTotal*2];

		for (size_t n = 0; n < NnTotal; n++) {
			// Left VOLUME
			double const rhoL  = *rhoL_ptr++,
			             rhouL = *rhouL_ptr++,
			             rhovL = *rhovL_ptr++,
			             EL    = *EL_ptr++,

			             rhoL_inv = 1.0/rhoL,
			             uL = rhouL*rhoL_inv,
			             vL = rhovL*rhoL_inv,

			             V2L = uL*uL+vL*vL,
			             VL  = sqrt(V2L),

			             pL  = GM1*(EL-0.5*rhoL*V2L),
			             cL  = sqrt(GAMMA*pL/rhoL);

			// Right VOLUME
			double const rhoR  = *rhoR_ptr++,
			             rhouR = *rhouR_ptr++,
			             rhovR = *rhovR_ptr++,
			             ER    = *ER_ptr++,

			             rhoR_inv = 1.0/rhoR,
			             uR = rhouR*rhoR_inv,
			             vR = rhovR*rhoR_inv,

			             V2R = uR*uR+vR*vR,
			             VR  = sqrt(V2R),

			             pR  = GM1*(ER-0.5*rhoR*V2R),
			             cR  = sqrt(GAMMA*pR/rhoR);

			double const maxlL = VL+cL,
			             maxlR = VR+cR;

			char sideMaxV;
			double maxV;
			if (maxlL > maxlR) {
				sideMaxV = 'L';
				maxV = maxlL;
			} else {
				sideMaxV = 'R';
				maxV = maxlR;
			}

			double const n1 = *n_ptr++,
			             n2 = *n_ptr++;

			if (nFluxNum != NULL) {
				double WLn[] = {rhoL, rhouL, rhovL, EL},
				       FLn[Neq*d];
				FLUXDATA->W = WLn;
				FLUXDATA->F = FLn;
				flux_Euler(FLUXDATA);

				double const *FL1_ptr = FLn,
				             *FL2_ptr = FL1_ptr+1;

				double WRn[] = {rhoR, rhouR, rhovR, ER},
				       FRn[Neq*d];
				FLUXDATA->W = WRn;
				FLUXDATA->F = FRn;
				flux_Euler(FLUXDATA);

				double const *FR1_ptr = FRn,
				             *FR2_ptr = FR1_ptr+1;

				size_t IndnF = 0;
				for (size_t eq = 0; eq < Neq; eq++) {
					*nF_ptr[IndnF++]++ = 0.5*( n1*((*FL1_ptr)+(*FR1_ptr))+n2*((*FL2_ptr)+(*FR2_ptr))
					                          + maxV*(WLn[eq]-WRn[eq]));

					FL1_ptr += d;
					FL2_ptr += d;
					FR1_ptr += d;
					FR2_ptr += d;
				}
			}

			double dFdWn[Nvar*Nvar*d];

			// Flux term
			double WLn[] = {rhoL, rhouL, rhovL, EL};

			FLUXDATA->W    = WLn;
			FLUXDATA->dFdW = dFdWn;
			jacobian_flux_Euler(FLUXDATA);

			double const *dF1dW_ptr = dFdWn,
			             *dF2dW_ptr = dF1dW_ptr+1;

			size_t InddnFdW = 0;
			for (size_t eq = 0; eq < Neq; eq++) {
			for (size_t var = 0; var < Nvar; var++) {
				*dnFdWL_ptr[InddnFdW++] = 0.5*(n1*(*dF1dW_ptr)+n2*(*dF2dW_ptr));

				dF1dW_ptr += d;
				dF2dW_ptr += d;
			}}

			// Dissipation term
			InddnFdW = 0;
			for (size_t eq = 0; eq < Neq; eq++) {
			for (size_t var = 0; var < Nvar; var++) {
				if (var == eq)
					*dnFdWL_ptr[InddnFdW] += 0.5*maxV;
				InddnFdW++;
			}}
			if (sideMaxV == 'L') {
				double const WRn[]    = { rhoR,         rhouR,     rhovR,    ER},
				             dudW[]   = {-uL*rhoL_inv,  rhoL_inv,  0.0,      0.0},
				             dvdW[]   = {-vL*rhoL_inv,  0.0,       rhoL_inv, 0.0},
				             drhodW[] = { 1.0,          0.0,       0.0,      0.0},
				             dpdW[]   = { GM1*0.5*V2L, -GM1*uL,   -GM1*vL,   GM1};

				double dmaxVdW[Nvar];
				for (size_t var = 0; var < Nvar; var++) {
					double const dVdW = 1.0/VL*(uL*dudW[var]+vL*dvdW[var]),
					             dcdW = GAMMA/(2.0*cL*rhoL*rhoL)*(dpdW[var]*rhoL-pL*drhodW[var]);
					dmaxVdW[var] = dVdW+dcdW;
				}

				InddnFdW = 0;
				for (size_t eq = 0; eq < Neq; eq++) {
				for (size_t var = 0; var < Nvar; var++) {
					*dnFdWL_ptr[InddnFdW++] += 0.5*dmaxVdW[var]*(WLn[eq]-WRn[eq]);
				}}
			}

			for (size_t i = 0, iMax = Neq*Nvar; i < iMax; i++)
				dnFdWL_ptr[i]++;

			if (dnFdWR != NULL) {
				// Flux term
				double const WRn[] = {rhoR, rhouR, rhovR, ER};

				FLUXDATA->W    = WRn;
				FLUXDATA->dFdW = dFdWn;
				jacobian_flux_Euler(FLUXDATA);

				double const *dF1dW_ptr = dFdWn,
				             *dF2dW_ptr = dF1dW_ptr+1;

				size_t InddnFdW = 0;
				for (size_t eq = 0; eq < Neq; eq++) {
				for (size_t var = 0; var < Nvar; var++) {
					*dnFdWR_ptr[InddnFdW++] = 0.5*(n1*(*dF1dW_ptr)+n2*(*dF2dW_ptr));

					dF1dW_ptr += d;
					dF2dW_ptr += d;
				}}

				// Dissipation term
				InddnFdW = 0;
				for (size_t eq = 0; eq < Neq; eq++) {
				for (size_t var = 0; var < Nvar; var++) {
					if (var == eq)
						*dnFdWR_ptr[InddnFdW] -= 0.5*maxV;
					InddnFdW++;
				}}
				if (sideMaxV == 'R') {
					double const WLn[]    = { rhoL,         rhouL,     rhovL,    EL},
					             dudW[]   = {-uR*rhoR_inv,  rhoR_inv,  0.0,      0.0},
					             dvdW[]   = {-vR*rhoR_inv,  0.0,       rhoR_inv, 0.0},
					             drhodW[] = { 1.0,          0.0,       0.0,      0.0},
					             dpdW[]   = { GM1*0.5*V2R, -GM1*uR,   -GM1*vR,   GM1};

					double dmaxVdW[Nvar];
					for (size_t var = 0; var < Nvar; var++) {
						double const dVdW = 1.0/VR*(uR*dudW[var]+vR*dvdW[var]),
						             dcdW = GAMMA/(2.0*cR*rhoR*rhoR)*(dpdW[var]*rhoR-pR*drhodW[var]);
						dmaxVdW[var] = dVdW+dcdW;
					}

					InddnFdW = 0;
					for (size_t eq = 0; eq < Neq; eq++) {
					for (size_t var = 0; var < Nvar; var++) {
						*dnFdWR_ptr[InddnFdW++] += 0.5*dmaxVdW[var]*(WLn[eq]-WRn[eq]);
					}}
				}

				for (size_t i = 0, iMax = Neq*Nvar; i < iMax; i++)
					dnFdWR_ptr[i]++;
			}
		}
	} else if (d == 1) {
		for (size_t n = 0; n < NnTotal; n++) {
			// Left VOLUME
			double const rhoL  = *rhoL_ptr++,
			             rhouL = *rhouL_ptr++,
			             EL    = *EL_ptr++,

			             rhoL_inv = 1.0/rhoL,
			             uL = rhouL*rhoL_inv,

			             V2L = uL*uL,
			             VL  = sqrt(V2L),

			             pL  = GM1*(EL-0.5*rhoL*V2L),
			             cL  = sqrt(GAMMA*pL/rhoL);

			// Right VOLUME
			double const rhoR  = *rhoR_ptr++,
			             rhouR = *rhouR_ptr++,
			             ER    = *ER_ptr++,

			             rhoR_inv = 1.0/rhoR,
			             uR = rhouR*rhoR_inv,

			             V2R = uR*uR,
			             VR  = sqrt(V2R),

			             pR  = GM1*(ER-0.5*rhoR*V2R),
			             cR  = sqrt(GAMMA*pR/rhoR);


			double const maxlL = VL+cL,
			             maxlR = VR+cR;

			char sideMaxV;
			double maxV;
			if (maxlL > maxlR) {
				sideMaxV = 'L';
				maxV = maxlL;
			} else {
				sideMaxV = 'R';
				maxV = maxlR;
			}

			double const n1 = *n_ptr++;

			if (nFluxNum != NULL) {
				double WLn[] = {rhoL, rhouL, EL},
				       FLn[Neq*d];
				FLUXDATA->W = WLn;
				FLUXDATA->F = FLn;
				flux_Euler(FLUXDATA);

				double const *FL1_ptr = FLn;

				double WRn[] = {rhoR, rhouR, ER},
				       FRn[Neq*d];
				FLUXDATA->W = WRn;
				FLUXDATA->F = FRn;
				flux_Euler(FLUXDATA);

				double const *FR1_ptr = FRn;

				size_t IndnF = 0;
				for (size_t eq = 0; eq < Neq; eq++) {
					*nF_ptr[IndnF++]++ = 0.5*( n1*((*FL1_ptr)+(*FR1_ptr)) + maxV*(WLn[eq]-WRn[eq]));

					FL1_ptr += d;
					FR1_ptr += d;
				}
			}

			double dFdWn[Nvar*Nvar*d];

			// Flux term
			double WLn[] = {rhoL, rhouL, EL};

			FLUXDATA->W    = WLn;
			FLUXDATA->dFdW = dFdWn;
			jacobian_flux_Euler(FLUXDATA);

			double const *dF1dW_ptr = dFdWn;

			size_t InddnFdW = 0;
			for (size_t eq = 0; eq < Neq; eq++) {
			for (size_t var = 0; var < Nvar; var++) {
				*dnFdWL_ptr[InddnFdW++] = 0.5*(n1*(*dF1dW_ptr));

				dF1dW_ptr += d;
			}}

			// Dissipation term
			InddnFdW = 0;
			for (size_t eq = 0; eq < Neq; eq++) {
			for (size_t var = 0; var < Nvar; var++) {
				if (var == eq)
					*dnFdWL_ptr[InddnFdW] += 0.5*maxV;
				InddnFdW++;
			}}
			if (sideMaxV == 'L') {
				double const WRn[]    = { rhoR,         rhouR,    ER},
				             dudW[]   = {-uL*rhoL_inv,  rhoL_inv, 0.0},
				             drhodW[] = { 1.0,          0.0,      0.0},
				             dpdW[]   = { GM1*0.5*V2L, -GM1*uL,   GM1};

				double dmaxVdW[Nvar];
				for (size_t var = 0; var < Nvar; var++) {
					double const dVdW = 1.0/VL*(uL*dudW[var]),
					             dcdW = GAMMA/(2.0*cL*rhoL*rhoL)*(dpdW[var]*rhoL-pL*drhodW[var]);
					dmaxVdW[var] = dVdW+dcdW;
				}

				InddnFdW = 0;
				for (size_t eq = 0; eq < Neq; eq++) {
				for (size_t var = 0; var < Nvar; var++) {
					*dnFdWL_ptr[InddnFdW++] += 0.5*dmaxVdW[var]*(WLn[eq]-WRn[eq]);
				}}
			}

			for (size_t i = 0, iMax = Neq*Nvar; i < iMax; i++)
				dnFdWL_ptr[i]++;

			if (dnFdWR != NULL) {
				// Flux term
				double const WRn[] = {rhoR, rhouR, ER};

				FLUXDATA->W    = WRn;
				FLUXDATA->dFdW = dFdWn;
				jacobian_flux_Euler(FLUXDATA);

				double const *dF1dW_ptr = dFdWn;

				size_t InddnFdW = 0;
				for (size_t eq = 0; eq < Neq; eq++) {
				for (size_t var = 0; var < Nvar; var++) {
					*dnFdWR_ptr[InddnFdW++] = 0.5*(n1*(*dF1dW_ptr));

					dF1dW_ptr += d;
				}}

				// Dissipation term
				InddnFdW = 0;
				for (size_t eq = 0; eq < Neq; eq++) {
				for (size_t var = 0; var < Nvar; var++) {
					if (var == eq)
						*dnFdWR_ptr[InddnFdW] -= 0.5*maxV;
					InddnFdW++;
				}}
				if (sideMaxV == 'R') {
					double const WLn[]    = { rhoL,         rhouL,    EL},
					             dudW[]   = {-uR*rhoR_inv,  rhoR_inv, 0.0},
					             drhodW[] = { 1.0,          0.0,      0.0},
					             dpdW[]   = { GM1*0.5*V2R, -GM1*uR,   GM1};

					double dmaxVdW[Nvar];
					for (size_t var = 0; var < Nvar; var++) {
						double const dVdW = 1.0/VR*(uR*dudW[var]),
						             dcdW = GAMMA/(2.0*cR*rhoR*rhoR)*(dpdW[var]*rhoR-pR*drhodW[var]);
						dmaxVdW[var] = dVdW+dcdW;
					}

					InddnFdW = 0;
					for (size_t eq = 0; eq < Neq; eq++) {
					for (size_t var = 0; var < Nvar; var++) {
						*dnFdWR_ptr[InddnFdW++] += 0.5*dmaxVdW[var]*(WLn[eq]-WRn[eq]);
					}}
				}

				for (size_t i = 0, iMax = Neq*Nvar; i < iMax; i++)
					dnFdWR_ptr[i]++;
			}
		}
	}
	free(FLUXDATA);
}

static void jacobian_flux_Roe(struct S_NUMERICALFLUX *const NUMFLUXDATA)
{
	unsigned int const d       = NUMFLUXDATA->d,
	                   Neq     = d+2,
	                   Nvar    = d+2,
	                   Nn      = NUMFLUXDATA->Nn,
	                   Nel     = NUMFLUXDATA->Nel,
	                   NnTotal = Nn*Nel;

	double const *const nL = NUMFLUXDATA->nL;

	double const *const WL = NUMFLUXDATA->WL,
	             *const WR = NUMFLUXDATA->WR;

	double *const nFluxNum = NUMFLUXDATA->nFluxNum,
	       *const dnFdWL   = NUMFLUXDATA->dnFluxNumdWL,
	       *const dnFdWR   = NUMFLUXDATA->dnFluxNumdWR;

	double const *rhoL_ptr  = &WL[NnTotal*0],
	             *rhouL_ptr = &WL[NnTotal*1],
	             *EL_ptr    = &WL[NnTotal*(d+1)],

	             *rhoR_ptr  = &WR[NnTotal*0],
	             *rhouR_ptr = &WR[NnTotal*1],
	             *ER_ptr    = &WR[NnTotal*(d+1)];

	double const *n_ptr = nL;

	double *nF_ptr[Neq];
	if (nFluxNum != NULL) {
		for (size_t eq = 0; eq < Neq; eq++)
			nF_ptr[eq] = &nFluxNum[eq*NnTotal];
	}

	double *dnFdWL_ptr[Neq*Neq];
	for (size_t eq = 0; eq < Neq; eq++) {
	for (size_t var = 0; var < Nvar; var++) {
		dnFdWL_ptr[eq*Nvar+var] = &dnFdWL[(eq*Nvar+var)*NnTotal];
	}}

	double *dnFdWR_ptr[Neq*Neq];
	if (dnFdWR != NULL) {
		for (size_t eq = 0; eq < Neq; eq++) {
		for (size_t var = 0; var < Nvar; var++) {
			dnFdWR_ptr[eq*Nvar+var] = &dnFdWR[(eq*Nvar+var)*NnTotal];
		}}
	}

	if (d == 3) {
		double const *rhovL_ptr = &WL[NnTotal*2],
		             *rhowL_ptr = &WL[NnTotal*3],

		             *rhovR_ptr = &WR[NnTotal*2],
		             *rhowR_ptr = &WR[NnTotal*3];

		for (size_t n = 0; n < NnTotal; n++) {
			double const n1 = *n_ptr++,
			             n2 = *n_ptr++,
			             n3 = *n_ptr++;

			// Left VOLUME
			double const rhoL  = *rhoL_ptr++,
			             rhouL = *rhouL_ptr++,
			             rhovL = *rhovL_ptr++,
			             rhowL = *rhowL_ptr++,
			             EL    = *EL_ptr++,

			             rhoL_inv  = 1.0/rhoL,
			             rhoLr_inv = sqrt(rhoL_inv),
			             uL = rhouL*rhoL_inv,
			             vL = rhovL*rhoL_inv,
			             wL = rhowL*rhoL_inv,

			             V2L = uL*uL+vL*vL+wL*wL,
			             VnL = n1*uL+n2*vL+n3*wL,

			             pL  = GM1*(EL-0.5*rhoL*V2L),
			             HL  = (EL+pL)*rhoL_inv;

			// Right VOLUME
			double const rhoR  = *rhoR_ptr++,
			             rhouR = *rhouR_ptr++,
			             rhovR = *rhovR_ptr++,
			             rhowR = *rhowR_ptr++,
			             ER    = *ER_ptr++,

			             rhoR_inv  = 1.0/rhoR,
			             rhoRr_inv = sqrt(rhoR_inv),
			             uR = rhouR*rhoR_inv,
			             vR = rhovR*rhoR_inv,
			             wR = rhowR*rhoR_inv,

			             V2R = uR*uR+vR*vR+wR*wR,
			             VnR = n1*uR+n2*vR+n3*wR,

			             pR  = GM1*(ER-0.5*rhoR*V2R),
			             HR  = (ER+pR)*rhoR_inv;

			// Roe-averaged states
			double const r   = sqrt(rhoR/rhoL),
			             rP1 = r+1.0,

			             rho = r*rhoL,
			             u   = (r*uR+uL)/rP1,
			             v   = (r*vR+vL)/rP1,
			             w   = (r*wR+wL)/rP1,
			             H   = (r*HR+HL)/rP1,
			             Vn  = n1*u+n2*v+n3*w,
			             V2  = u*u+v*v+w*w,
			             c2  = GM1*(H-0.5*V2),
			             c   = sqrt(c2),

			             Den  = sqrt(rhoL) + sqrt(rhoR);

			// Compute eigenvalues (with entropy fix)
			unsigned int case_l1;
			double const l1L     = VnL-c;
			double       l1      = Vn-c,
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

			unsigned int case_l5;
			double const l5R     = VnR+c;
			double       l5      = Vn+c,
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

			double sign_l234 = 1.0;
			if (Vn < 0.0)
				sign_l234 = -1.0;
			double const l234 = sign_l234*(Vn);


			double const drho  = rhoR-rhoL,
			             drhou = rhoR*uR-rhoL*uL,
			             drhov = rhoR*vR-rhoL*vL,
			             drhow = rhoR*wR-rhoL*wL,
			             dE    = ER-EL,
			             dp    = pR-pL,
			             dVn   = VnR-VnL,

			             lc1 = 0.5*(l5+l1) - l234,
			             lc2 = 0.5*(l5-l1),

			             disInter1 = lc1*dp/(c*c) + lc2*rho*dVn/c,
			             disInter2 = lc1*rho*dVn  + lc2*dp/c;

			if (nFluxNum != NULL) {
				double const rhoVnL = rhoL*VnL,
				             rhoVnR = rhoR*VnR,
				             pLR    = pL + pR;

				double const nF1 = rhoVnL      + rhoVnR,
				             nF2 = rhoVnL*uL   + rhoVnR*uR  + n1*pLR,
				             nF3 = rhoVnL*vL   + rhoVnR*vR  + n2*pLR,
				             nF4 = rhoVnL*wL   + rhoVnR*wR  + n3*pLR,
				             nF5 = VnL*(EL+pL) + VnR*(ER+pR),

				             dis1 = l234*drho  + disInter1,
				             dis2 = l234*drhou + disInter1*u + disInter2*n1,
				             dis3 = l234*drhov + disInter1*v + disInter2*n2,
				             dis4 = l234*drhow + disInter1*w + disInter2*n3,
				             dis5 = l234*dE    + disInter1*H + disInter2*(Vn);

				// Assemble components
				size_t IndnF = 0;
				*nF_ptr[IndnF++]++ = 0.5*(nF1 - dis1);
				*nF_ptr[IndnF++]++ = 0.5*(nF2 - dis2);
				*nF_ptr[IndnF++]++ = 0.5*(nF3 - dis3);
				*nF_ptr[IndnF++]++ = 0.5*(nF4 - dis4);
				*nF_ptr[IndnF++]++ = 0.5*(nF5 - dis5);
			}

			double dnF1dW[Neq],  dnF2dW[Neq],  dnF3dW[Neq],  dnF4dW[Neq],  dnF5dW[Neq],
			       ddis1dW[Neq], ddis2dW[Neq], ddis3dW[Neq], ddis4dW[Neq], ddis5dW[Neq];

			// Flux term
			double const duLdW[]   = {-uL*rhoL_inv,  rhoL_inv,  0.0,       0.0,      0.0},
			             dvLdW[]   = {-vL*rhoL_inv,  0.0,       rhoL_inv,  0.0,      0.0},
			             dwLdW[]   = {-wL*rhoL_inv,  0.0,       0.0,       rhoL_inv, 0.0},
			             drhoLdW[] = { 1.0,          0.0,       0.0,       0.0,      0.0},
			             dELdW[]   = { 0.0,          0.0,       0.0,       0.0,      1.0},
			             dpLdW[]   = { GM1*0.5*V2L, -GM1*uL,   -GM1*vL,   -GM1*wL,   GM1};

			double const rhoVn = rhoL*VnL;
			double dVnLdW[Nvar];
			for (size_t var = 0; var < Nvar; var++) {
				dVnLdW[var] = n1*duLdW[var]+n2*dvLdW[var]+n3*dwLdW[var];
				double const drhoVndW = drhoLdW[var]*VnL + rhoL*dVnLdW[var];

				dnF1dW[var] = drhoVndW;
				dnF2dW[var] = drhoVndW*uL + rhoVn*duLdW[var] + n1*dpLdW[var];
				dnF3dW[var] = drhoVndW*vL + rhoVn*dvLdW[var] + n2*dpLdW[var];
				dnF4dW[var] = drhoVndW*wL + rhoVn*dwLdW[var] + n3*dpLdW[var];
				dnF5dW[var] = dVnLdW[var]*(EL+pL) + VnL*(dELdW[var]+dpLdW[var]);
			}

			// Dissipation term
			double const mult = rhoLr_inv/Den;
			double const drhodW[] = { 0.5*rhoR/rho,             0.0,          0.0,          0.0,         0.0},
			             dudW[]   = {-0.5*(uL+u)*mult,          mult,         0.0,          0.0,         0.0},
			             dvdW[]   = {-0.5*(vL+v)*mult,          0.0,          mult,         0.0,         0.0},
			             dwdW[]   = {-0.5*(wL+w)*mult,          0.0,          0.0,          mult,        0.0},
			             dHdW[]   = {-0.5*(HL+H-GM1*V2L)*mult, -GM1*uL*mult, -GM1*vL*mult, -GM1*wL*mult, GAMMA*mult};

			for (size_t var = 0; var < Nvar; var++) {
				double const dcdW  = 0.5*GM1/c*(dHdW[var]-(u*dudW[var]+v*dvdW[var]+w*dwdW[var])),
				             dVndW = n1*dudW[var]+n2*dvdW[var]+n3*dwdW[var];

				double dl1dW;
				if (case_l1)
					dl1dW = sign_l1*(dVndW-dcdW);
				else
					dl1dW = sign_l1*(dVnLdW[var]-dcdW);

				double dl5dW;
				if (case_l5)
					dl5dW = sign_l5*(dVndW+dcdW);
				else
					dl5dW = sign_l5*(dcdW);

				double const dl234dW = sign_l234*dVndW,
				             dlc1dW = 0.5*(dl5dW+dl1dW) - dl234dW,
				             dlc2dW = 0.5*(dl5dW-dl1dW),

				             ddisInter1dW = (dlc1dW*dp-lc1*dpLdW[var])/c2-(2.0*lc1*dp*dcdW)/(c*c2) +
				                            (dlc2dW*rho*dVn+lc2*drhodW[var]*dVn-lc2*rho*dVnLdW[var])/c -
				                            (lc2*rho*dVn*dcdW)/c2,
				             ddisInter2dW = dlc1dW*rho*dVn+lc1*drhodW[var]*dVn-lc1*rho*dVnLdW[var] +
				                            (dlc2dW*dp-lc2*dpLdW[var])/c-(lc2*dp*dcdW)/c2,

				             drhoudW = drhoLdW[var]*uL+rhoL*duLdW[var],
				             drhovdW = drhoLdW[var]*vL+rhoL*dvLdW[var],
				             drhowdW = drhoLdW[var]*wL+rhoL*dwLdW[var];

				ddis1dW[var] = dl234dW*drho - l234*drhoLdW[var] + ddisInter1dW;
				ddis2dW[var] = dl234dW*drhou - l234*drhoudW + ddisInter1dW*u + disInter1*dudW[var] + ddisInter2dW*n1;
				ddis3dW[var] = dl234dW*drhov - l234*drhovdW + ddisInter1dW*v + disInter1*dvdW[var] + ddisInter2dW*n2;
				ddis4dW[var] = dl234dW*drhow - l234*drhowdW + ddisInter1dW*w + disInter1*dwdW[var] + ddisInter2dW*n3;
				ddis5dW[var] = dl234dW*dE - l234*dELdW[var] + ddisInter1dW*H  + disInter1*dHdW[var]
				                                            + ddisInter2dW*Vn + disInter2*dVndW;
			}

			size_t InddnFdW = 0;
			for (size_t var = 0; var < Nvar; var++) *dnFdWL_ptr[InddnFdW++]++ = 0.5*(dnF1dW[var]-ddis1dW[var]);
			for (size_t var = 0; var < Nvar; var++) *dnFdWL_ptr[InddnFdW++]++ = 0.5*(dnF2dW[var]-ddis2dW[var]);
			for (size_t var = 0; var < Nvar; var++) *dnFdWL_ptr[InddnFdW++]++ = 0.5*(dnF3dW[var]-ddis3dW[var]);
			for (size_t var = 0; var < Nvar; var++) *dnFdWL_ptr[InddnFdW++]++ = 0.5*(dnF4dW[var]-ddis4dW[var]);
			for (size_t var = 0; var < Nvar; var++) *dnFdWL_ptr[InddnFdW++]++ = 0.5*(dnF5dW[var]-ddis5dW[var]);

			if (dnFdWR != NULL) {
				// Flux term
				double const duRdW[]   = {-uR*rhoR_inv,  rhoR_inv,  0.0,       0.0,      0.0},
				             dvRdW[]   = {-vR*rhoR_inv,  0.0,       rhoR_inv,  0.0,      0.0},
				             dwRdW[]   = {-wR*rhoR_inv,  0.0,       0.0,       rhoR_inv, 0.0},
				             drhoRdW[] = { 1.0,          0.0,       0.0,       0.0,      0.0},
				             dERdW[]   = { 0.0,          0.0,       0.0,       0.0,      1.0},
				             dpRdW[]   = { GM1*0.5*V2R, -GM1*uR,   -GM1*vR,   -GM1*wR,   GM1};

				double const rhoVn = rhoR*VnR;
				double dVnRdW[Nvar];
				for (size_t var = 0; var < Nvar; var++) {
					dVnRdW[var] = n1*duRdW[var]+n2*dvRdW[var]+n3*dwRdW[var];
					double const drhoVndW = drhoRdW[var]*VnR + rhoR*dVnRdW[var];

					dnF1dW[var] = drhoVndW;
					dnF2dW[var] = drhoVndW*uR + rhoVn*duRdW[var] + n1*dpRdW[var];
					dnF3dW[var] = drhoVndW*vR + rhoVn*dvRdW[var] + n2*dpRdW[var];
					dnF4dW[var] = drhoVndW*wR + rhoVn*dwRdW[var] + n3*dpRdW[var];
					dnF5dW[var] = dVnRdW[var]*(ER+pR) + VnR*(dERdW[var]+dpRdW[var]);
				}

				// Dissipation term
				double const mult = rhoRr_inv/Den;
				double const drhodW[] = { 0.5*rhoL/rho,             0.0,          0.0,          0.0,         0.0},
				             dudW[]   = {-0.5*(uR+u)*mult,          mult,         0.0,          0.0,         0.0},
				             dvdW[]   = {-0.5*(vR+v)*mult,          0.0,          mult,         0.0,         0.0},
				             dwdW[]   = {-0.5*(wR+w)*mult,          0.0,          0.0,          mult,        0.0},
				             dHdW[]   = {-0.5*(HR+H-GM1*V2R)*mult, -GM1*uR*mult, -GM1*vR*mult, -GM1*wR*mult, GAMMA*mult};

				for (size_t var = 0; var < Nvar; var++) {
					double const dcdW  = 0.5*GM1/c*(dHdW[var]-(u*dudW[var]+v*dvdW[var]+w*dwdW[var])),
					             dVndW = n1*dudW[var]+n2*dvdW[var]+n3*dwdW[var];

					double dl1dW;
					if (case_l1)
						dl1dW = sign_l1*(dVndW-dcdW);
					else
						dl1dW = sign_l1*(-dcdW);

					double dl5dW;
					if (case_l5)
						dl5dW = sign_l5*(dVndW+dcdW);
					else
						dl5dW = sign_l5*(dVnRdW[var]+dcdW);

					double const dl234dW = sign_l234*dVndW,
					             dlc1dW = 0.5*(dl5dW+dl1dW) - dl234dW,
					             dlc2dW = 0.5*(dl5dW-dl1dW),

					             ddisInter1dW = (dlc1dW*dp+lc1*dpRdW[var])/c2-(2.0*lc1*dp*dcdW)/(c*c2) +
					                            (dlc2dW*rho*dVn+lc2*drhodW[var]*dVn+lc2*rho*dVnRdW[var])/c -
					                            (lc2*rho*dVn*dcdW)/c2,
					             ddisInter2dW = dlc1dW*rho*dVn+lc1*drhodW[var]*dVn+lc1*rho*dVnRdW[var] +
					                            (dlc2dW*dp+lc2*dpRdW[var])/c-(lc2*dp*dcdW)/c2,

					             drhoudW = drhoRdW[var]*uR+rhoR*duRdW[var],
					             drhovdW = drhoRdW[var]*vR+rhoR*dvRdW[var],
					             drhowdW = drhoRdW[var]*wR+rhoR*dwRdW[var];

					ddis1dW[var] = dl234dW*drho + l234*drhoRdW[var] + ddisInter1dW;
					ddis2dW[var] = dl234dW*drhou + l234*drhoudW + ddisInter1dW*u + disInter1*dudW[var] + ddisInter2dW*n1;
					ddis3dW[var] = dl234dW*drhov + l234*drhovdW + ddisInter1dW*v + disInter1*dvdW[var] + ddisInter2dW*n2;
					ddis4dW[var] = dl234dW*drhow + l234*drhowdW + ddisInter1dW*w + disInter1*dwdW[var] + ddisInter2dW*n3;
					ddis5dW[var] = dl234dW*dE + l234*dERdW[var] + ddisInter1dW*H  + disInter1*dHdW[var]
					                                            + ddisInter2dW*Vn + disInter2*dVndW;
				}

				size_t InddnFdW = 0;
				for (size_t var = 0; var < Nvar; var++) *dnFdWR_ptr[InddnFdW++]++ = 0.5*(dnF1dW[var]-ddis1dW[var]);
				for (size_t var = 0; var < Nvar; var++) *dnFdWR_ptr[InddnFdW++]++ = 0.5*(dnF2dW[var]-ddis2dW[var]);
				for (size_t var = 0; var < Nvar; var++) *dnFdWR_ptr[InddnFdW++]++ = 0.5*(dnF3dW[var]-ddis3dW[var]);
				for (size_t var = 0; var < Nvar; var++) *dnFdWR_ptr[InddnFdW++]++ = 0.5*(dnF4dW[var]-ddis4dW[var]);
				for (size_t var = 0; var < Nvar; var++) *dnFdWR_ptr[InddnFdW++]++ = 0.5*(dnF5dW[var]-ddis5dW[var]);
			}
		}
	} else if (d == 2) {
		double const *rhovL_ptr = &WL[NnTotal*2],

		             *rhovR_ptr = &WR[NnTotal*2];

		for (size_t n = 0; n < NnTotal; n++) {
			double const n1 = *n_ptr++,
			             n2 = *n_ptr++;

			// Left VOLUME
			double const rhoL  = *rhoL_ptr++,
			             rhouL = *rhouL_ptr++,
			             rhovL = *rhovL_ptr++,
			             EL    = *EL_ptr++,

			             rhoL_inv  = 1.0/rhoL,
			             rhoLr_inv = sqrt(rhoL_inv),
			             uL = rhouL*rhoL_inv,
			             vL = rhovL*rhoL_inv,

			             V2L = uL*uL+vL*vL,
			             VnL = n1*uL+n2*vL,

			             pL  = GM1*(EL-0.5*rhoL*V2L),
			             HL  = (EL+pL)*rhoL_inv;

			// Right VOLUME
			double const rhoR  = *rhoR_ptr++,
			             rhouR = *rhouR_ptr++,
			             rhovR = *rhovR_ptr++,
			             ER    = *ER_ptr++,

			             rhoR_inv  = 1.0/rhoR,
			             rhoRr_inv = sqrt(rhoR_inv),
			             uR = rhouR*rhoR_inv,
			             vR = rhovR*rhoR_inv,

			             V2R = uR*uR+vR*vR,
			             VnR = n1*uR+n2*vR,

			             pR  = GM1*(ER-0.5*rhoR*V2R),
			             HR  = (ER+pR)*rhoR_inv;

			// Roe-averaged states
			double const r   = sqrt(rhoR/rhoL),
			             rP1 = r+1.0,

			             rho = r*rhoL,
			             u   = (r*uR+uL)/rP1,
			             v   = (r*vR+vL)/rP1,
			             H   = (r*HR+HL)/rP1,
			             Vn  = n1*u+n2*v,
			             V2  = u*u+v*v,
			             c2  = GM1*(H-0.5*V2),
			             c   = sqrt(c2),

			             Den  = sqrt(rhoL) + sqrt(rhoR);

			// Compute eigenvalues (with entropy fix)
			unsigned int case_l1;
			double const l1L     = VnL-c;
			double       l1      = Vn-c,
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

			unsigned int case_l5;
			double const l5R     = VnR+c;
			double       l5      = Vn+c,
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

			double sign_l234 = 1.0;
			if (Vn < 0.0)
				sign_l234 = -1.0;
			double const l234 = sign_l234*(Vn);


			double const drho  = rhoR-rhoL,
			             drhou = rhoR*uR-rhoL*uL,
			             drhov = rhoR*vR-rhoL*vL,
			             dE    = ER-EL,
			             dp    = pR-pL,
			             dVn   = VnR-VnL,

			             lc1 = 0.5*(l5+l1) - l234,
			             lc2 = 0.5*(l5-l1),

			             disInter1 = lc1*dp/(c*c) + lc2*rho*dVn/c,
			             disInter2 = lc1*rho*dVn  + lc2*dp/c;

			if (nFluxNum != NULL) {
				double const rhoVnL = rhoL*VnL,
				             rhoVnR = rhoR*VnR,
				             pLR    = pL + pR;

				double const nF1 = rhoVnL      + rhoVnR,
				             nF2 = rhoVnL*uL   + rhoVnR*uR  + n1*pLR,
				             nF3 = rhoVnL*vL   + rhoVnR*vR  + n2*pLR,
				             nF5 = VnL*(EL+pL) + VnR*(ER+pR),

				             dis1 = l234*drho  + disInter1,
				             dis2 = l234*drhou + disInter1*u + disInter2*n1,
				             dis3 = l234*drhov + disInter1*v + disInter2*n2,
				             dis5 = l234*dE    + disInter1*H + disInter2*(Vn);

				// Assemble components
				size_t IndnF = 0;
				*nF_ptr[IndnF++]++ = 0.5*(nF1 - dis1);
				*nF_ptr[IndnF++]++ = 0.5*(nF2 - dis2);
				*nF_ptr[IndnF++]++ = 0.5*(nF3 - dis3);
				*nF_ptr[IndnF++]++ = 0.5*(nF5 - dis5);
			}

			double dnF1dW[Neq],  dnF2dW[Neq],  dnF3dW[Neq],  dnF5dW[Neq],
			       ddis1dW[Neq], ddis2dW[Neq], ddis3dW[Neq], ddis5dW[Neq];

			// Flux term
			double const duLdW[]   = {-uL*rhoL_inv,  rhoL_inv,  0.0,      0.0},
			             dvLdW[]   = {-vL*rhoL_inv,  0.0,       rhoL_inv, 0.0},
			             drhoLdW[] = { 1.0,          0.0,       0.0,      0.0},
			             dELdW[]   = { 0.0,          0.0,       0.0,      1.0},
			             dpLdW[]   = { GM1*0.5*V2L, -GM1*uL,   -GM1*vL,   GM1};

			double const rhoVn = rhoL*VnL;
			double dVnLdW[Nvar];
			for (size_t var = 0; var < Nvar; var++) {
				dVnLdW[var] = n1*duLdW[var]+n2*dvLdW[var];
				double const drhoVndW = drhoLdW[var]*VnL + rhoL*dVnLdW[var];

				dnF1dW[var] = drhoVndW;
				dnF2dW[var] = drhoVndW*uL + rhoVn*duLdW[var] + n1*dpLdW[var];
				dnF3dW[var] = drhoVndW*vL + rhoVn*dvLdW[var] + n2*dpLdW[var];
				dnF5dW[var] = dVnLdW[var]*(EL+pL) + VnL*(dELdW[var]+dpLdW[var]);
			}

			// Dissipation term
			double const mult = rhoLr_inv/Den;
			double const drhodW[] = { 0.5*rhoR/rho,             0.0,          0.0,         0.0},
			             dudW[]   = {-0.5*(uL+u)*mult,          mult,         0.0,         0.0},
			             dvdW[]   = {-0.5*(vL+v)*mult,          0.0,          mult,        0.0},
			             dHdW[]   = {-0.5*(HL+H-GM1*V2L)*mult, -GM1*uL*mult, -GM1*vL*mult, GAMMA*mult};

			for (size_t var = 0; var < Nvar; var++) {
				double const dcdW  = 0.5*GM1/c*(dHdW[var]-(u*dudW[var]+v*dvdW[var])),
				             dVndW = n1*dudW[var]+n2*dvdW[var];

				double dl1dW;
				if (case_l1)
					dl1dW = sign_l1*(dVndW-dcdW);
				else
					dl1dW = sign_l1*(dVnLdW[var]-dcdW);

				double dl5dW;
				if (case_l5)
					dl5dW = sign_l5*(dVndW+dcdW);
				else
					dl5dW = sign_l5*(dcdW);

				double const dl234dW = sign_l234*dVndW,
				             dlc1dW = 0.5*(dl5dW+dl1dW) - dl234dW,
				             dlc2dW = 0.5*(dl5dW-dl1dW),

				             ddisInter1dW = (dlc1dW*dp-lc1*dpLdW[var])/c2-(2.0*lc1*dp*dcdW)/(c*c2) +
				                            (dlc2dW*rho*dVn+lc2*drhodW[var]*dVn-lc2*rho*dVnLdW[var])/c -
				                            (lc2*rho*dVn*dcdW)/c2,
				             ddisInter2dW = dlc1dW*rho*dVn+lc1*drhodW[var]*dVn-lc1*rho*dVnLdW[var] +
				                            (dlc2dW*dp-lc2*dpLdW[var])/c-(lc2*dp*dcdW)/c2,

				             drhoudW = drhoLdW[var]*uL+rhoL*duLdW[var],
				             drhovdW = drhoLdW[var]*vL+rhoL*dvLdW[var];

				ddis1dW[var] = dl234dW*drho - l234*drhoLdW[var] + ddisInter1dW;
				ddis2dW[var] = dl234dW*drhou - l234*drhoudW + ddisInter1dW*u + disInter1*dudW[var] + ddisInter2dW*n1;
				ddis3dW[var] = dl234dW*drhov - l234*drhovdW + ddisInter1dW*v + disInter1*dvdW[var] + ddisInter2dW*n2;
				ddis5dW[var] = dl234dW*dE - l234*dELdW[var] + ddisInter1dW*H  + disInter1*dHdW[var]
				                                            + ddisInter2dW*Vn + disInter2*dVndW;
			}

			size_t InddnFdW = 0;
			for (size_t var = 0; var < Nvar; var++) *dnFdWL_ptr[InddnFdW++]++ = 0.5*(dnF1dW[var]-ddis1dW[var]);
			for (size_t var = 0; var < Nvar; var++) *dnFdWL_ptr[InddnFdW++]++ = 0.5*(dnF2dW[var]-ddis2dW[var]);
			for (size_t var = 0; var < Nvar; var++) *dnFdWL_ptr[InddnFdW++]++ = 0.5*(dnF3dW[var]-ddis3dW[var]);
			for (size_t var = 0; var < Nvar; var++) *dnFdWL_ptr[InddnFdW++]++ = 0.5*(dnF5dW[var]-ddis5dW[var]);

			if (dnFdWR != NULL) {
				// Flux term
				double const duRdW[]   = {-uR*rhoR_inv,  rhoR_inv,  0.0,      0.0},
				             dvRdW[]   = {-vR*rhoR_inv,  0.0,       rhoR_inv, 0.0},
				             drhoRdW[] = { 1.0,          0.0,       0.0,      0.0},
				             dERdW[]   = { 0.0,          0.0,       0.0,      1.0},
				             dpRdW[]   = { GM1*0.5*V2R, -GM1*uR,   -GM1*vR,   GM1};

				double const rhoVn = rhoR*VnR;
				double dVnRdW[Nvar];
				for (size_t var = 0; var < Nvar; var++) {
					dVnRdW[var] = n1*duRdW[var]+n2*dvRdW[var];
					double const drhoVndW = drhoRdW[var]*VnR + rhoR*dVnRdW[var];

					dnF1dW[var] = drhoVndW;
					dnF2dW[var] = drhoVndW*uR + rhoVn*duRdW[var] + n1*dpRdW[var];
					dnF3dW[var] = drhoVndW*vR + rhoVn*dvRdW[var] + n2*dpRdW[var];
					dnF5dW[var] = dVnRdW[var]*(ER+pR) + VnR*(dERdW[var]+dpRdW[var]);
				}

				// Dissipation term
				double const mult = rhoRr_inv/Den;
				double const drhodW[] = { 0.5*rhoL/rho,             0.0,          0.0,         0.0},
				             dudW[]   = {-0.5*(uR+u)*mult,          mult,         0.0,         0.0},
				             dvdW[]   = {-0.5*(vR+v)*mult,          0.0,          mult,        0.0},
				             dHdW[]   = {-0.5*(HR+H-GM1*V2R)*mult, -GM1*uR*mult, -GM1*vR*mult, GAMMA*mult};

				for (size_t var = 0; var < Nvar; var++) {
					double const dcdW  = 0.5*GM1/c*(dHdW[var]-(u*dudW[var]+v*dvdW[var])),
					             dVndW = n1*dudW[var]+n2*dvdW[var];

					double dl1dW;
					if (case_l1)
						dl1dW = sign_l1*(dVndW-dcdW);
					else
						dl1dW = sign_l1*(-dcdW);

					double dl5dW;
					if (case_l5)
						dl5dW = sign_l5*(dVndW+dcdW);
					else
						dl5dW = sign_l5*(dVnRdW[var]+dcdW);

					double const dl234dW = sign_l234*dVndW,
					             dlc1dW = 0.5*(dl5dW+dl1dW) - dl234dW,
					             dlc2dW = 0.5*(dl5dW-dl1dW),

					             ddisInter1dW = (dlc1dW*dp+lc1*dpRdW[var])/c2-(2.0*lc1*dp*dcdW)/(c*c2) +
					                            (dlc2dW*rho*dVn+lc2*drhodW[var]*dVn+lc2*rho*dVnRdW[var])/c -
					                            (lc2*rho*dVn*dcdW)/c2,
					             ddisInter2dW = dlc1dW*rho*dVn+lc1*drhodW[var]*dVn+lc1*rho*dVnRdW[var] +
					                            (dlc2dW*dp+lc2*dpRdW[var])/c-(lc2*dp*dcdW)/c2,

					             drhoudW = drhoRdW[var]*uR+rhoR*duRdW[var],
					             drhovdW = drhoRdW[var]*vR+rhoR*dvRdW[var];

					ddis1dW[var] = dl234dW*drho + l234*drhoRdW[var] + ddisInter1dW;
					ddis2dW[var] = dl234dW*drhou + l234*drhoudW + ddisInter1dW*u + disInter1*dudW[var] + ddisInter2dW*n1;
					ddis3dW[var] = dl234dW*drhov + l234*drhovdW + ddisInter1dW*v + disInter1*dvdW[var] + ddisInter2dW*n2;
					ddis5dW[var] = dl234dW*dE + l234*dERdW[var] + ddisInter1dW*H  + disInter1*dHdW[var]
					                                            + ddisInter2dW*Vn + disInter2*dVndW;
				}

				size_t InddnFdW = 0;
				for (size_t var = 0; var < Nvar; var++) *dnFdWR_ptr[InddnFdW++]++ = 0.5*(dnF1dW[var]-ddis1dW[var]);
				for (size_t var = 0; var < Nvar; var++) *dnFdWR_ptr[InddnFdW++]++ = 0.5*(dnF2dW[var]-ddis2dW[var]);
				for (size_t var = 0; var < Nvar; var++) *dnFdWR_ptr[InddnFdW++]++ = 0.5*(dnF3dW[var]-ddis3dW[var]);
				for (size_t var = 0; var < Nvar; var++) *dnFdWR_ptr[InddnFdW++]++ = 0.5*(dnF5dW[var]-ddis5dW[var]);
			}
		}
	} else if (d == 1) {
		for (size_t n = 0; n < NnTotal; n++) {
			double const n1 = *n_ptr++;

			// Left VOLUME
			double const rhoL  = *rhoL_ptr++,
			             rhouL = *rhouL_ptr++,
			             EL    = *EL_ptr++,

			             rhoL_inv  = 1.0/rhoL,
			             rhoLr_inv = sqrt(rhoL_inv),
			             uL = rhouL*rhoL_inv,

			             V2L = uL*uL,
			             VnL = n1*uL,

			             pL  = GM1*(EL-0.5*rhoL*V2L),
			             HL  = (EL+pL)*rhoL_inv;

			// Right VOLUME
			double const rhoR  = *rhoR_ptr++,
			             rhouR = *rhouR_ptr++,
			             ER    = *ER_ptr++,

			             rhoR_inv  = 1.0/rhoR,
			             rhoRr_inv = sqrt(rhoR_inv),
			             uR = rhouR*rhoR_inv,

			             V2R = uR*uR,
			             VnR = n1*uR,

			             pR  = GM1*(ER-0.5*rhoR*V2R),
			             HR  = (ER+pR)*rhoR_inv;

			// Roe-averaged states
			double const r   = sqrt(rhoR/rhoL),
			             rP1 = r+1.0,

			             rho = r*rhoL,
			             u   = (r*uR+uL)/rP1,
			             H   = (r*HR+HL)/rP1,
			             Vn  = n1*u,
			             V2  = u*u,
			             c2  = GM1*(H-0.5*V2),
			             c   = sqrt(c2),

			             Den  = sqrt(rhoL) + sqrt(rhoR);

			// Compute eigenvalues (with entropy fix)
			unsigned int case_l1;
			double const l1L     = VnL-c;
			double       l1      = Vn-c,
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

			unsigned int case_l5;
			double const l5R     = VnR+c;
			double       l5      = Vn+c,
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

			double sign_l234 = 1.0;
			if (Vn < 0.0)
				sign_l234 = -1.0;
			double const l234 = sign_l234*(Vn);


			double const drho  = rhoR-rhoL,
			             drhou = rhoR*uR-rhoL*uL,
			             dE    = ER-EL,
			             dp    = pR-pL,
			             dVn   = VnR-VnL,

			             lc1 = 0.5*(l5+l1) - l234,
			             lc2 = 0.5*(l5-l1),

			             disInter1 = lc1*dp/(c*c) + lc2*rho*dVn/c,
			             disInter2 = lc1*rho*dVn  + lc2*dp/c;

			if (nFluxNum != NULL) {
				double const rhoVnL = rhoL*VnL,
				             rhoVnR = rhoR*VnR,
				             pLR    = pL + pR;

				double const nF1 = rhoVnL      + rhoVnR,
				             nF2 = rhoVnL*uL   + rhoVnR*uR  + n1*pLR,
				             nF5 = VnL*(EL+pL) + VnR*(ER+pR),

				             dis1 = l234*drho  + disInter1,
				             dis2 = l234*drhou + disInter1*u + disInter2*n1,
				             dis5 = l234*dE    + disInter1*H + disInter2*(Vn);

				// Assemble components
				size_t IndnF = 0;
				*nF_ptr[IndnF++]++ = 0.5*(nF1 - dis1);
				*nF_ptr[IndnF++]++ = 0.5*(nF2 - dis2);
				*nF_ptr[IndnF++]++ = 0.5*(nF5 - dis5);
			}

			double dnF1dW[Neq],  dnF2dW[Neq],  dnF5dW[Neq],
			       ddis1dW[Neq], ddis2dW[Neq], ddis5dW[Neq];

			// Flux term
			double const duLdW[]   = {-uL*rhoL_inv,  rhoL_inv, 0.0},
			             drhoLdW[] = { 1.0,          0.0,      0.0},
			             dELdW[]   = { 0.0,          0.0,      1.0},
			             dpLdW[]   = { GM1*0.5*V2L, -GM1*uL,   GM1};

			double const rhoVn = rhoL*VnL;
			double dVnLdW[Nvar];
			for (size_t var = 0; var < Nvar; var++) {
				dVnLdW[var] = n1*duLdW[var];
				double const drhoVndW = drhoLdW[var]*VnL + rhoL*dVnLdW[var];

				dnF1dW[var] = drhoVndW;
				dnF2dW[var] = drhoVndW*uL + rhoVn*duLdW[var] + n1*dpLdW[var];
				dnF5dW[var] = dVnLdW[var]*(EL+pL) + VnL*(dELdW[var]+dpLdW[var]);
			}

			// Dissipation term
			double const mult = rhoLr_inv/Den;
			double const drhodW[] = { 0.5*rhoR/rho,             0.0,         0.0},
			             dudW[]   = {-0.5*(uL+u)*mult,          mult,        0.0},
			             dHdW[]   = {-0.5*(HL+H-GM1*V2L)*mult, -GM1*uL*mult, GAMMA*mult};

			for (size_t var = 0; var < Nvar; var++) {
				double const dcdW  = 0.5*GM1/c*(dHdW[var]-(u*dudW[var])),
				             dVndW = n1*dudW[var];

				double dl1dW;
				if (case_l1)
					dl1dW = sign_l1*(dVndW-dcdW);
				else
					dl1dW = sign_l1*(dVnLdW[var]-dcdW);

				double dl5dW;
				if (case_l5)
					dl5dW = sign_l5*(dVndW+dcdW);
				else
					dl5dW = sign_l5*(dcdW);

				double const dl234dW = sign_l234*dVndW,
				             dlc1dW = 0.5*(dl5dW+dl1dW) - dl234dW,
				             dlc2dW = 0.5*(dl5dW-dl1dW),

				             ddisInter1dW = (dlc1dW*dp-lc1*dpLdW[var])/c2-(2.0*lc1*dp*dcdW)/(c*c2) +
				                            (dlc2dW*rho*dVn+lc2*drhodW[var]*dVn-lc2*rho*dVnLdW[var])/c -
				                            (lc2*rho*dVn*dcdW)/c2,
				             ddisInter2dW = dlc1dW*rho*dVn+lc1*drhodW[var]*dVn-lc1*rho*dVnLdW[var] +
				                            (dlc2dW*dp-lc2*dpLdW[var])/c-(lc2*dp*dcdW)/c2,

				             drhoudW = drhoLdW[var]*uL+rhoL*duLdW[var];

				ddis1dW[var] = dl234dW*drho - l234*drhoLdW[var] + ddisInter1dW;
				ddis2dW[var] = dl234dW*drhou - l234*drhoudW + ddisInter1dW*u + disInter1*dudW[var] + ddisInter2dW*n1;
				ddis5dW[var] = dl234dW*dE - l234*dELdW[var] + ddisInter1dW*H  + disInter1*dHdW[var]
				                                            + ddisInter2dW*Vn + disInter2*dVndW;
			}

			size_t InddnFdW = 0;
			for (size_t var = 0; var < Nvar; var++) *dnFdWL_ptr[InddnFdW++]++ = 0.5*(dnF1dW[var]-ddis1dW[var]);
			for (size_t var = 0; var < Nvar; var++) *dnFdWL_ptr[InddnFdW++]++ = 0.5*(dnF2dW[var]-ddis2dW[var]);
			for (size_t var = 0; var < Nvar; var++) *dnFdWL_ptr[InddnFdW++]++ = 0.5*(dnF5dW[var]-ddis5dW[var]);

			if (dnFdWR != NULL) {
				// Flux term
				double const duRdW[]   = {-uR*rhoR_inv,  rhoR_inv, 0.0},
				             drhoRdW[] = { 1.0,          0.0,      0.0},
				             dERdW[]   = { 0.0,          0.0,      1.0},
				             dpRdW[]   = { GM1*0.5*V2R, -GM1*uR,   GM1};

				double const rhoVn = rhoR*VnR;
				double dVnRdW[Nvar];
				for (size_t var = 0; var < Nvar; var++) {
					dVnRdW[var] = n1*duRdW[var];
					double const drhoVndW = drhoRdW[var]*VnR + rhoR*dVnRdW[var];

					dnF1dW[var] = drhoVndW;
					dnF2dW[var] = drhoVndW*uR + rhoVn*duRdW[var] + n1*dpRdW[var];
					dnF5dW[var] = dVnRdW[var]*(ER+pR) + VnR*(dERdW[var]+dpRdW[var]);
				}

				// Dissipation term
				double const mult = rhoRr_inv/Den;
				double const drhodW[] = { 0.5*rhoL/rho,             0.0,         0.0},
				             dudW[]   = {-0.5*(uR+u)*mult,          mult,        0.0},
				             dHdW[]   = {-0.5*(HR+H-GM1*V2R)*mult, -GM1*uR*mult, GAMMA*mult};

				for (size_t var = 0; var < Nvar; var++) {
					double const dcdW  = 0.5*GM1/c*(dHdW[var]-(u*dudW[var])),
					             dVndW = n1*dudW[var];

					double dl1dW;
					if (case_l1)
						dl1dW = sign_l1*(dVndW-dcdW);
					else
						dl1dW = sign_l1*(-dcdW);

					double dl5dW;
					if (case_l5)
						dl5dW = sign_l5*(dVndW+dcdW);
					else
						dl5dW = sign_l5*(dVnRdW[var]+dcdW);

					double const dl234dW = sign_l234*dVndW,
					             dlc1dW = 0.5*(dl5dW+dl1dW) - dl234dW,
					             dlc2dW = 0.5*(dl5dW-dl1dW),

					             ddisInter1dW = (dlc1dW*dp+lc1*dpRdW[var])/c2-(2.0*lc1*dp*dcdW)/(c*c2) +
					                            (dlc2dW*rho*dVn+lc2*drhodW[var]*dVn+lc2*rho*dVnRdW[var])/c -
					                            (lc2*rho*dVn*dcdW)/c2,
					             ddisInter2dW = dlc1dW*rho*dVn+lc1*drhodW[var]*dVn+lc1*rho*dVnRdW[var] +
					                            (dlc2dW*dp+lc2*dpRdW[var])/c-(lc2*dp*dcdW)/c2,

					             drhoudW = drhoRdW[var]*uR+rhoR*duRdW[var];

					ddis1dW[var] = dl234dW*drho + l234*drhoRdW[var] + ddisInter1dW;
					ddis2dW[var] = dl234dW*drhou + l234*drhoudW + ddisInter1dW*u + disInter1*dudW[var] + ddisInter2dW*n1;
					ddis5dW[var] = dl234dW*dE + l234*dERdW[var] + ddisInter1dW*H  + disInter1*dHdW[var]
					                                            + ddisInter2dW*Vn + disInter2*dVndW;
				}

				size_t InddnFdW = 0;
				for (size_t var = 0; var < Nvar; var++) *dnFdWR_ptr[InddnFdW++]++ = 0.5*(dnF1dW[var]-ddis1dW[var]);
				for (size_t var = 0; var < Nvar; var++) *dnFdWR_ptr[InddnFdW++]++ = 0.5*(dnF2dW[var]-ddis2dW[var]);
				for (size_t var = 0; var < Nvar; var++) *dnFdWR_ptr[InddnFdW++]++ = 0.5*(dnF5dW[var]-ddis5dW[var]);
			}
		}
	}

}

static void jacobian_flux_upwind(struct S_NUMERICALFLUX *const NUMFLUXDATA)
{
	unsigned int const d       = NUMFLUXDATA->d,
	                   Nn      = NUMFLUXDATA->Nn,
	                   Nel     = NUMFLUXDATA->Nel,
	                   NnTotal = Nn*Nel;

	double const *const nL = NUMFLUXDATA->nL;

	double const *const WL = NUMFLUXDATA->WL,
	             *const WR = NUMFLUXDATA->WR;

	double       *const nFluxNum     = NUMFLUXDATA->nFluxNum,
	             *const dnFluxNumdWL = NUMFLUXDATA->dnFluxNumdWL,
	             *const dnFluxNumdWR = NUMFLUXDATA->dnFluxNumdWR;

	double const *const b = compute_b_Advection(NnTotal,NUMFLUXDATA->XYZ); // free

	for (size_t n = 0; n < NnTotal; n++) {
		double b_dot_n = 0.0;
		for (size_t dim = 0; dim < d; dim++)
			b_dot_n += b[dim*NnTotal+n]*nL[n*d+dim];

		if (nFluxNum != NULL) {
			if (b_dot_n >= 0.0)
				nFluxNum[n] = b_dot_n*WL[n];
			else
				nFluxNum[n] = b_dot_n*WR[n];
		}

		if (b_dot_n >= 0.0)
			dnFluxNumdWL[n] = b_dot_n;
		else
			dnFluxNumdWL[n] = 0.0;

		if (dnFluxNumdWR != NULL) {
			if (b_dot_n >= 0.0)
				dnFluxNumdWR[n] = 0.0;
			else
				dnFluxNumdWR[n] = b_dot_n;
		}
	}
	free((double *) b);
}
