// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "fluxes_viscous.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"

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

void flux_viscous(const unsigned int Nn, const unsigned int Nel, const double *const W, const double *const *const Q,
                  double *const F)
{
	/*
	 *	Comments:
	 *		The storage ordering of the fluxes (Node, then dimension, then equation) is chosen such that memory stride
	 *		is minimized when converting from physical to reference space.
	 *
	 *		Assumptions:
	 *			Stokes hypothesis: Coefficient of bulk viscosity = -2/3*mu
	 *
	 *			kappa = mu*Cp/Pr  (1)
	 *			Cp/Rg = GAMMA/GM1 (2)
	 *
	 *		(1), (2) -> q = -kappa*Grad(T) == -mu/Pr*GAMMA*Grad(E/rho-0.5*V^2)
	 */

	const unsigned int d   = DB.d,
	                   Neq = DB.Neq;
	const double       Pr  = DB.Pr;

	if (!(d == 2 || d == 3))
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
			const double rho      = *rho_ptr++;
			const double rho_inv  = 1.0/rho;
			const double rho_inv2 = rho_inv*rho_inv;

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
				EXIT_UNSUPPORTED;
			}

			const double tau[DMAX][DMAX] = { { mu*2.0*(du[0]-divV/3.0), mu*(dv[0]+du[1]),        mu*(dw[0]+du[2]), },
			                                 { mu*(du[1]+dv[0]),        mu*2.0*(dv[1]-divV/3.0), mu*(dw[1]+dv[2]), },
			                                 { mu*(du[2]+dw[0]),        mu*(dv[2]+dw[1]),        mu*2.0*(dw[2]-divV/3.0), }, };

			const double dEoRho[DMAX] = { rho_inv2*(dE[0]*rho-E*drho[0]),
			                              rho_inv2*(dE[1]*rho-E*drho[1]),
			                              rho_inv2*(dE[2]*rho-E*drho[2]), },
			             dV2[DMAX]    = { 2.0*(u*du[0]+v*dv[0]+w*dw[0]),
			                              2.0*(u*du[1]+v*dv[1]+w*dw[1]),
			                              2.0*(u*du[2]+v*dv[2]+w*dw[2]), };

			const double dTs[DMAX] = { dEoRho[0]-0.5*dV2[0], // (T)emperature (s)caled
			                           dEoRho[1]-0.5*dV2[1],
			                           dEoRho[2]-0.5*dV2[2], };

			size_t IndF = 0;
			// eq 1
			*F_ptr[IndF++] = 0.0;
			*F_ptr[IndF++] = 0.0;
			*F_ptr[IndF++] = 0.0;

			// eq 2
			*F_ptr[IndF++] = tau[0][0];
			*F_ptr[IndF++] = tau[0][1];
			*F_ptr[IndF++] = tau[0][2];

			// eq 3
			*F_ptr[IndF++] = tau[1][0];
			*F_ptr[IndF++] = tau[1][1];
			*F_ptr[IndF++] = tau[1][2];

			// eq 4
			*F_ptr[IndF++] = tau[2][0];
			*F_ptr[IndF++] = tau[2][1];
			*F_ptr[IndF++] = tau[2][2];

			// eq 5
			*F_ptr[IndF++] = u*tau[0][0]+v*tau[0][1]+w*tau[0][2] + mu*GAMMA/Pr*dTs[0];
			*F_ptr[IndF++] = u*tau[1][0]+v*tau[1][1]+w*tau[1][2] + mu*GAMMA/Pr*dTs[1];
			*F_ptr[IndF++] = u*tau[2][0]+v*tau[2][1]+w*tau[2][2] + mu*GAMMA/Pr*dTs[2];

			for (size_t i = 0, iMax = Neq*DMAX; i < iMax; i++)
				F_ptr[i]++;
		}
	} else if (d == 2) {
		for (size_t n = 0; n < NnTotal; n++) {
			const double rho      = *rho_ptr++;
			const double rho_inv  = 1.0/rho;
			const double rho_inv2 = rho_inv*rho_inv;

			const double u = (*rhou_ptr++)*rho_inv,
			             v = (*rhov_ptr++)*rho_inv,
			             E = *E_ptr++;

			const double drho[DMAX]  = { *drho_ptr[0]++,  *drho_ptr[1]++,  },
			             drhou[DMAX] = { *drhou_ptr[0]++, *drhou_ptr[1]++, },
			             drhov[DMAX] = { *drhov_ptr[0]++, *drhov_ptr[1]++, },
			             dE[DMAX]    = { *dE_ptr[0]++,    *dE_ptr[1]++,    };

			const double du[DMAX] = { rho_inv*(drhou[0]-drho[0]*u), rho_inv*(drhou[1]-drho[1]*u), },
			             dv[DMAX] = { rho_inv*(drhov[0]-drho[0]*v), rho_inv*(drhov[1]-drho[1]*v), };

			const double divV = du[0]+dv[1];

			double mu = 0.0;
			if (DB.Const_mu) {
				mu = DB.mu;
			} else {
				// Sutherland's formula
				EXIT_UNSUPPORTED;
			}

			const double tau[DMAX][DMAX] = { { mu*2.0*(du[0]-divV/3.0), mu*(dv[0]+du[1]),        },
			                                 { mu*(du[1]+dv[0]),        mu*2.0*(dv[1]-divV/3.0), }, };

			const double dEoRho[DMAX] = { rho_inv2*(dE[0]*rho-E*drho[0]), rho_inv2*(dE[1]*rho-E*drho[1]), },
			             dV2[DMAX]    = { 2.0*(u*du[0]+v*dv[0]),          2.0*(u*du[1]+v*dv[1]), };

			const double dTs[DMAX] = { dEoRho[0]-0.5*dV2[0], dEoRho[1]-0.5*dV2[1], };

			size_t IndF = 0;
			// eq 1
			*F_ptr[IndF++] = 0.0;
			*F_ptr[IndF++] = 0.0;
			IndF += 1;

			// eq 2
			*F_ptr[IndF++] = tau[0][0];
			*F_ptr[IndF++] = tau[0][1];
			IndF += 1;

			// eq 3
			*F_ptr[IndF++] = tau[1][0];
			*F_ptr[IndF++] = tau[1][1];
			IndF += 1;

			// eq 4
			*F_ptr[IndF++] = u*tau[0][0]+v*tau[0][1]+ mu*GAMMA/Pr*dTs[0];
			*F_ptr[IndF++] = u*tau[1][0]+v*tau[1][1]+ mu*GAMMA/Pr*dTs[1];

			for (size_t i = 0, iMax = Neq*DMAX; i < iMax; i++)
				F_ptr[i]++;
		}
	} else {
		EXIT_UNSUPPORTED;
	}
}
