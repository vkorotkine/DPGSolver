// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "fluxes_viscous.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"

#include "fluxes_structs.h"

/*
 *	Purpose:
 *		Compute viscous fluxes from inputs W and Q in conservative form.
 *
 *	Comments:
 *		It is assumed that inputs: W, Q, n and outputs: F, nFluxNum are vectorized (i.e. the memory ordering is by
 *		equation and not by element).
 *
 *	Notation:
 *
 *	References:
 */

void flux_viscous(struct S_FLUX *const FLUXDATA)
{
	/*
	 *	Comments:
	 *		The negated viscous flux is returned such that the same operators can be used as for the inviscid
	 *		contribution.
	 *
	 *		The storage ordering of the fluxes (Node, then dimension, then equation) is chosen such that memory stride
	 *		is minimized when converting from physical to reference space.
	 *
	 *		Assumptions:
	 *      	Coefficient of bulk viscosity = 0.0
	 *
	 *			kappa = mu*Cp/Pr  (1)
	 *			Cp/Rg = GAMMA/GM1 (2)
	 *
	 *		(1), (2) -> q = -kappa*Grad(T) == -mu/Pr*GAMMA*Grad(E/rho-0.5*V^2)
	 */

	unsigned int const d   = FLUXDATA->d,
	                   Neq = d+2,
	                   Nn  = FLUXDATA->Nn,
	                   Nel = FLUXDATA->Nel;

	double const *const W        = FLUXDATA->W,
	             *const *const Q = FLUXDATA->Q;
	double       *const F        = FLUXDATA->F;

	const double Pr = DB.Pr;

	if (!(d == 2 || d == 3))
		EXIT_UNSUPPORTED;

	if (DB.Pr == 0.0 || DB.mu == 0.0)
		EXIT_UNSUPPORTED;

	const unsigned int NnTotal = Nn*Nel;

	const double *rho_ptr  = &W[NnTotal*0],
	             *rhou_ptr = &W[NnTotal*1],
	             *rhov_ptr = &W[NnTotal*2],
	             *E_ptr    = &W[NnTotal*(d+1)],
	             *drho_ptr[DMAX], *drhou_ptr[DMAX], *drhov_ptr[DMAX], *dE_ptr[DMAX];

	for (size_t dim = 0; dim < d; dim++) {
		drho_ptr[dim]  = &Q[dim][NnTotal*0];
		drhou_ptr[dim] = &Q[dim][NnTotal*1];
		drhov_ptr[dim] = &Q[dim][NnTotal*2];
		dE_ptr[dim]    = &Q[dim][NnTotal*(d+1)];
	}

	double *F_ptr[DMAX*Neq];
	for (size_t eq = 0; eq < Neq; eq++) {
	for (size_t dim = 0; dim < d; dim++) {
		F_ptr[eq*DMAX+dim] = &F[(eq*d+dim)*NnTotal];
	}}

	if (d == 3) {
		const double *rhow_ptr = &W[NnTotal*d],
		             *drhow_ptr[DMAX];
		for (size_t dim = 0; dim < d; dim++)
			drhow_ptr[dim] = &Q[dim][NnTotal*3];

		for (size_t n = 0; n < NnTotal; n++) {
			const double rho      = *rho_ptr++,
			             rho_inv  = 1.0/rho,
			             rho_inv2 = rho_inv*rho_inv;

			const double u = (*rhou_ptr++)*rho_inv,
			             v = (*rhov_ptr++)*rho_inv,
			             w = (*rhow_ptr++)*rho_inv,
			             E = *E_ptr++;

			const double drho[DMAX]  = { *drho_ptr[0]++,  *drho_ptr[1]++,  *drho_ptr[2]++,  },
			             drhou[DMAX] = { *drhou_ptr[0]++, *drhou_ptr[1]++, *drhou_ptr[2]++, },
			             drhov[DMAX] = { *drhov_ptr[0]++, *drhov_ptr[1]++, *drhov_ptr[2]++, },
			             drhow[DMAX] = { *drhow_ptr[0]++, *drhow_ptr[1]++, *drhow_ptr[2]++, },
			             dE[DMAX]    = { *dE_ptr[0]++,    *dE_ptr[1]++,    *dE_ptr[2]++,    };

			const double du[DMAX] = { rho_inv*(drhou[0]-drho[0]*u),
			                          rho_inv*(drhou[1]-drho[1]*u),
			                          rho_inv*(drhou[2]-drho[2]*u), },
			             dv[DMAX] = { rho_inv*(drhov[0]-drho[0]*v),
			                          rho_inv*(drhov[1]-drho[1]*v),
			                          rho_inv*(drhov[2]-drho[2]*v), },
			             dw[DMAX] = { rho_inv*(drhow[0]-drho[0]*w),
			                          rho_inv*(drhow[1]-drho[1]*w),
			                          rho_inv*(drhow[2]-drho[2]*w), };

			const double divV = du[0]+dv[1]+dw[2];

			double mu = 0.0;
			if (DB.Const_mu) {
				mu = DB.mu;
			} else {
				// Sutherland's formula
				// Be sure to test both mu configurations.
				EXIT_UNSUPPORTED;
			}

			double tau[d][d];
			tau[0][0] = mu*2.0*(du[0]-divV/3.0);
			tau[0][1] = mu*(dv[0]+du[1]);
			tau[0][2] = mu*(dw[0]+du[2]);
			tau[1][0] = tau[0][1];
			tau[1][1] = mu*2.0*(dv[1]-divV/3.0);
			tau[1][2] = mu*(dw[1]+dv[2]);
			tau[2][0] = tau[0][2];
			tau[2][1] = tau[1][2];
			tau[2][2] = mu*2.0*(dw[2]-divV/3.0);

			const double dEoRho[DMAX] = { rho_inv2*(dE[0]*rho-E*drho[0]),
			                              rho_inv2*(dE[1]*rho-E*drho[1]),
			                              rho_inv2*(dE[2]*rho-E*drho[2]), },
			             dV2[DMAX]    = { 2.0*(u*du[0]+v*dv[0]+w*dw[0]),
			                              2.0*(u*du[1]+v*dv[1]+w*dw[1]),
			                              2.0*(u*du[2]+v*dv[2]+w*dw[2]), };

			 // (T)emperature (s)caled
			const double dTs[DMAX] = { dEoRho[0]-0.5*dV2[0], dEoRho[1]-0.5*dV2[1], dEoRho[2]-0.5*dV2[2], };

			size_t IndF = 0;
			// eq 1
			for (size_t dim = 0; dim < d; dim++)
				*F_ptr[IndF++] = 0.0;

			// eq 2
			for (size_t dim = 0; dim < d; dim++)
				*F_ptr[IndF++] = -(tau[0][dim]);

			// eq 3
			for (size_t dim = 0; dim < d; dim++)
				*F_ptr[IndF++] = -(tau[1][dim]);

			// eq 4
			for (size_t dim = 0; dim < d; dim++)
				*F_ptr[IndF++] = -(tau[2][dim]);

			// eq 5
			for (size_t dim = 0; dim < d; dim++)
				*F_ptr[IndF++] = -(u*tau[dim][0]+v*tau[dim][1]+w*tau[dim][2] + mu*GAMMA/Pr*dTs[dim]);

			for (size_t i = 0, iMax = Neq*DMAX; i < iMax; i++)
				F_ptr[i]++;
		}
	} else if (d == 2) {
		for (size_t n = 0; n < NnTotal; n++) {
			double const rho      = *rho_ptr++,
			             rho_inv  = 1.0/rho,
			             rho_inv2 = rho_inv*rho_inv,
			             u = (*rhou_ptr++)*rho_inv,
			             v = (*rhov_ptr++)*rho_inv,
			             E = *E_ptr++;

			double const drho[DMAX]  = { *drho_ptr[0]++,  *drho_ptr[1]++,  0.0, },
			             drhou[DMAX] = { *drhou_ptr[0]++, *drhou_ptr[1]++, 0.0, },
			             drhov[DMAX] = { *drhov_ptr[0]++, *drhov_ptr[1]++, 0.0, },
			             dE[DMAX]    = { *dE_ptr[0]++,    *dE_ptr[1]++,    0.0, };

			const double du[DMAX] = { rho_inv*(drhou[0]-drho[0]*u), rho_inv*(drhou[1]-drho[1]*u), 0.0, },
			             dv[DMAX] = { rho_inv*(drhov[0]-drho[0]*v), rho_inv*(drhov[1]-drho[1]*v), 0.0, };

			const double divV = du[0]+dv[1];

			double mu = 0.0;
			if (DB.Const_mu) {
				mu = DB.mu;
			} else {
				// Sutherland's formula
				// Be sure to test both mu configurations.
				EXIT_UNSUPPORTED;
			}

			double tau[d][d];
			tau[0][0] = mu*2.0*(du[0]-divV/3.0);
			tau[0][1] = mu*(dv[0]+du[1]);
			tau[1][0] = tau[0][1];
			tau[1][1] = mu*2.0*(dv[1]-divV/3.0);

			const double dEoRho[DMAX] = { rho_inv2*(dE[0]*rho-E*drho[0]), rho_inv2*(dE[1]*rho-E*drho[1]), },
			             dV2[DMAX]    = { 2.0*(u*du[0]+v*dv[0]),          2.0*(u*du[1]+v*dv[1]), };

			const double dTs[DMAX] = { dEoRho[0]-0.5*dV2[0], dEoRho[1]-0.5*dV2[1], };

			size_t IndF = 0;
			// eq 1
			for (size_t dim = 0; dim < d; dim++)
				*F_ptr[IndF++] = 0.0;
			IndF += 1;

			// eq 2
			for (size_t dim = 0; dim < d; dim++)
				*F_ptr[IndF++] = -(tau[0][dim]);
			IndF += 1;

			// eq 3
			for (size_t dim = 0; dim < d; dim++)
				*F_ptr[IndF++] = -(tau[1][dim]);
			IndF += 1;

			// eq 4
			for (size_t dim = 0; dim < d; dim++)
				*F_ptr[IndF++] = -(u*tau[dim][0]+v*tau[dim][1] + mu*GAMMA/Pr*dTs[dim]);

			for (size_t i = 0, iMax = Neq*DMAX; i < iMax; i++)
				F_ptr[i]++;
		}
	} else {
		EXIT_UNSUPPORTED;
	}
}
