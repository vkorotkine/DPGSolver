// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "jacobian_fluxes_inviscid.h"

/*
 *	Purpose:
 *		Compute inviscid flux jacobians from input W in conservative form.
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

void jacobian_flux_inviscid(const unsigned int Nn, const unsigned int Nel, double *W, double *dFdW,
                            const unsigned int d, const unsigned int Neq)
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

	// Standard datatypes
	unsigned int i, n, eq, var, dim, iMax, Nvar, NnTotal, InddFdW;
	double *rho_ptr, *rhou_ptr, *rhov_ptr, *rhow_ptr, *E_ptr,
	       rho, u, v, w, u2, uv, uw, v2, vw, w2, V2, E, p, H, alpha, beta, *dFdW_ptr[DMAX*Neq*Neq];

	Nvar    = Neq;
	NnTotal = Nn*Nel;

	rho_ptr  = &W[NnTotal*0];
	rhou_ptr = &W[NnTotal*1];
	E_ptr    = &W[NnTotal*(d+1)];

	for (eq  = 0; eq  < Neq;  eq++)  {
	for (var = 0; var < Nvar; var++) {
	for (dim = 0; dim < d;    dim++) {
		dFdW_ptr[(eq*Nvar+var)*DMAX+dim] = &dFdW[((eq*Nvar+var)*d+dim)*NnTotal];
	}}}

	if (d == 3) {
		rhov_ptr = &W[NnTotal*2];
		rhow_ptr = &W[NnTotal*3];

		for (n = 0; n < NnTotal; n++) {
			rho = *rho_ptr;
			u   = (*rhou_ptr)/rho;
			v   = (*rhov_ptr)/rho;
			w   = (*rhow_ptr)/rho;
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
			u   = (*rhou_ptr)/rho;
			v   = (*rhov_ptr)/rho;
			E   = *E_ptr;

			u2 = u*u;
			uv = u*v;
			v2 = v*v;

			V2 = u2+v2;
			p  = GM1*(E-0.5*rho*V2);
			H  = (E+p)/rho;

			alpha = 0.5*GM1*V2;
			beta  = alpha-H;

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
			u   = (*rhou_ptr)/rho;
			E   = *E_ptr;

			u2 = u*u;

			V2 = u2;
			p  = GM1*(E-0.5*rho*V2);
			H  = (E+p)/rho;

			alpha = 0.5*GM1*V2;
			beta  = alpha-H;

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
