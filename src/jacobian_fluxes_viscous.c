// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "jacobian_fluxes_viscous.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"

#include "fluxes_structs.h"

/*
 *	Purpose:
 *		Compute viscous flux jacobians from inputs W, Q in conservative form.
 *
 *	Comments:
 *
 *		*** IMPORTANT ***   *** IMPORTANT ***   *** IMPORTANT ***
 *
 *		Returned viscous fluxes are negated. This is done so that the same functions can be used for the inviscid and
 *		viscous contributions.
 *
 *		*** IMPORTANT ***   *** IMPORTANT ***   *** IMPORTANT ***
 *
 *		See comments of fluxes_viscous.
 *		See comments of jacobian_fluxes_inviscid.
 *
 *		Jacobians with respect to W and Q are both returned.
 *
 *	Notation:
 *
 *	References:
 */

static void jacobian_flux_Poisson      (struct S_FLUX *const FLUXDATA);
static void jacobian_flux_NavierStokes (struct S_FLUX *const FLUXDATA);

void jacobian_flux_viscous(struct S_FLUX *const FLUXDATA)
{
	switch(FLUXDATA->PDE_index) {
		case PDE_POISSON:      jacobian_flux_Poisson(FLUXDATA);      break;
		case PDE_NAVIERSTOKES: jacobian_flux_NavierStokes(FLUXDATA); break;
		default:
			printf("%d\n",FLUXDATA->PDE_index);
			EXIT_UNSUPPORTED;
		break;
	}
}

static void jacobian_flux_Poisson (struct S_FLUX *const FLUXDATA)
{
	/*
	 *	Comments:
	 *		Implicitly assumed that Neq = Nvar = 1.
	 *		F(W,Q) == Q => dFdW = 0, dFdQ = 1.
	 */

	unsigned int const d       = FLUXDATA->d,
	                   Nn      = FLUXDATA->Nn,
	                   Nel     = FLUXDATA->Nel,
	                   NnTotal = Nn*Nel;

	double const *const *const Q    = FLUXDATA->Q;
	double       *const F           = FLUXDATA->F,
	             *const dFdW        = FLUXDATA->dFdW,
	             *const *const dFdQ = FLUXDATA->dFdQ;

	double *F_ptr[DMAX];
	if (F != NULL) {
		for (size_t dim = 0; dim < d; dim++)
			F_ptr[dim] = &F[dim*NnTotal];
	}

	// dFdW == 0
	if (dFdW != NULL)
		EXIT_UNSUPPORTED;

	double *dFdQ_ptr[d][DMAX];
	if (dFdQ != NULL) {
		for (size_t dim1 = 0; dim1 < d; dim1++) {
		for (size_t dim = 0; dim < d; dim++) {
			dFdQ_ptr[dim1][dim] = &dFdQ[dim1][dim*NnTotal];
		}}
	}

	for (size_t n = 0; n < NnTotal; n++) {
		// ***************************************** F ***************************************** //
		if (F != NULL) {
			size_t IndF = 0;
			for (size_t dim = 0; dim < d; dim++)
				*F_ptr[IndF++] = -Q[dim][n];

			for (size_t i = 0, iMax = DMAX; i < iMax; i++)
				F_ptr[i]++;
		}

		// ***************************************** dFdQ ***************************************** //
		if (dFdQ != NULL) {
			for (size_t dim1 = 0; dim1 < d; dim1++) {
				size_t InddFdQ = 0;
				for (size_t dim = 0; dim < d; dim++) {
					if (dim == dim1)
						*dFdQ_ptr[dim1][InddFdQ++] = -1.0;
					else
						*dFdQ_ptr[dim1][InddFdQ++] = 0.0;
				}

				for (size_t i = 0, iMax = DMAX; i < iMax; i++)
					dFdQ_ptr[dim1][i]++;
			}
		}
	}
}

static void jacobian_flux_NavierStokes (struct S_FLUX *const FLUXDATA)
{
	unsigned int const d       = FLUXDATA->d,
	                   Neq     = d+2,
	                   Nvar    = d+2,
	                   Nn      = FLUXDATA->Nn,
	                   Nel     = FLUXDATA->Nel,
	                   NnTotal = Nn*Nel;

	double const *const W           = FLUXDATA->W,
	             *const *const Q    = FLUXDATA->Q;
	double       *const F           = FLUXDATA->F,
	             *const dFdW        = FLUXDATA->dFdW,
	             *const *const dFdQ = FLUXDATA->dFdQ;

	const double Pr = DB.Pr;

	if (!(d == 2 || d == 3))
		EXIT_UNSUPPORTED;

	if (DB.Pr == 0.0 || DB.mu == 0.0)
		EXIT_UNSUPPORTED;

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
	if (F != NULL) {
		for (size_t eq = 0; eq < Neq; eq++) {
		for (size_t dim = 0; dim < d; dim++) {
			F_ptr[eq*DMAX+dim] = &F[(eq*d+dim)*NnTotal];
		}}
	}

	double *dFdW_ptr[DMAX*Neq*Nvar];
	if (dFdW != NULL) {
		for (size_t eq  = 0; eq  < Neq;  eq++)  {
		for (size_t var = 0; var < Nvar; var++) {
		for (size_t dim = 0; dim < d;    dim++) {
			dFdW_ptr[(eq*Nvar+var)*DMAX+dim] = &dFdW[((eq*Nvar+var)*d+dim)*NnTotal];
		}}}
	}

	double *dFdQ_ptr[d][DMAX*Neq*Nvar];
	if (dFdQ != NULL) {
		for (size_t dim1 = 0; dim1 < d; dim1++) {
		for (size_t eq  = 0; eq  < Neq;  eq++)  {
		for (size_t var = 0; var < Nvar; var++) {
		for (size_t dim = 0; dim < d;    dim++) {
			dFdQ_ptr[dim1][(eq*Nvar+var)*DMAX+dim] = &dFdQ[dim1][((eq*Nvar+var)*d+dim)*NnTotal];
		}}}}
	}

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

			const double dTs[DMAX] = { dEoRho[0]-0.5*dV2[0], dEoRho[1]-0.5*dV2[1], dEoRho[2]-0.5*dV2[2], };


			// ***************************************** F ***************************************** //
			if (F != NULL) {
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

			// ***************************************** dFdW ***************************************** //
			if (dFdW != NULL) {
				double drhodW[]      = { 1.0,                  0.0, 0.0, 0.0, 0.0 },
				       drho_invdW[]  = {-rho_inv2,             0.0, 0.0, 0.0, 0.0 },
				       drho_inv2dW[] = {-2.0*rho_inv2*rho_inv, 0.0, 0.0, 0.0, 0.0 },

				       dudW[] = { -u*rho_inv, rho_inv, 0.0,     0.0,     0.0 },
				       dvdW[] = { -v*rho_inv, 0.0,     rho_inv, 0.0,     0.0 },
				       dwdW[] = { -w*rho_inv, 0.0,     0.0,     rho_inv, 0.0 },

				       dEdW[] = { 0.0, 0.0, 0.0, 0.0, 1.0 };

				double ddudW[d][Nvar], ddvdW[d][Nvar], ddwdW[d][Nvar];
				for (size_t var = 0; var < Nvar; var++) {
				for (size_t dim = 0; dim < d; dim++) {
					ddudW[dim][var] = drho_invdW[var]*(drhou[dim]-drho[dim]*u)+rho_inv*(0.0-drho[dim]*dudW[var]);
					ddvdW[dim][var] = drho_invdW[var]*(drhov[dim]-drho[dim]*v)+rho_inv*(0.0-drho[dim]*dvdW[var]);
					ddwdW[dim][var] = drho_invdW[var]*(drhow[dim]-drho[dim]*w)+rho_inv*(0.0-drho[dim]*dwdW[var]);
				}}

				double ddivVdW[Nvar];
				for (size_t var = 0; var < Nvar; var++)
					ddivVdW[var] = ddudW[0][var]+ddvdW[1][var]+ddwdW[2][var];

				double dtaudW[d][d][Nvar];
				for (size_t var = 0; var < Nvar; var++) {
					dtaudW[0][0][var] = mu*2.0*(ddudW[0][var]-ddivVdW[var]/3.0);
					dtaudW[0][1][var] = mu*(ddvdW[0][var]+ddudW[1][var]);
					dtaudW[0][2][var] = mu*(ddwdW[0][var]+ddudW[2][var]);
					dtaudW[1][0][var] = dtaudW[0][1][var];
					dtaudW[1][1][var] = mu*2.0*(ddvdW[1][var]-ddivVdW[var]/3.0);
					dtaudW[1][2][var] = mu*(ddwdW[1][var]+ddvdW[2][var]);
					dtaudW[2][0][var] = dtaudW[0][2][var];
					dtaudW[2][1][var] = dtaudW[1][2][var];
					dtaudW[2][2][var] = mu*2.0*(ddwdW[2][var]-ddivVdW[var]/3.0);
				}

				double dmudW[NVAR3D] = {0.0};
				if (!DB.Const_mu) {
					EXIT_UNSUPPORTED;

					for (size_t var = 0; var < Nvar; var++) {
						dtaudW[0][0][var] += dmudW[var]*2.0*(du[0]-divV/3.0);
						dtaudW[0][1][var] += dmudW[var]*(dv[0]+du[1]);
						dtaudW[0][2][var] += dmudW[var]*(dw[0]+du[2]);
						dtaudW[1][0][var]  = dtaudW[0][1][var];
						dtaudW[1][1][var] += dmudW[var]*2.0*(dv[1]-divV/3.0);
						dtaudW[1][2][var] += dmudW[var]*(dw[1]+dv[2]);
						dtaudW[2][0][var]  = dtaudW[0][2][var];
						dtaudW[2][1][var]  = dtaudW[1][2][var];
						dtaudW[2][2][var] += dmudW[var]*2.0*(dw[2]-divV/3.0);
					}
				}

				double ddTsdW[d][Nvar];
				for (size_t var = 0; var < Nvar; var++) {
				for (size_t dim = 0; dim < d; dim++) {
					double const ddEoRhodW = drho_inv2dW[var]*(dE[dim]*rho-E*drho[dim])
					                        +rho_inv2*((dE[dim]*drhodW[var])-(dEdW[var]*drho[dim])),
					             ddV2dW    = 2.0*( dudW[var]*du[dim]+dvdW[var]*dv[dim]+dwdW[var]*dw[dim]
					                              +u*ddudW[dim][var]+v*ddvdW[dim][var]+w*ddwdW[dim][var]);
					ddTsdW[dim][var] = ddEoRhodW-0.5*ddV2dW;
				}}


				size_t InddFdW = 0;
				// *** eq 1 ***
				for (size_t var = 0; var < Nvar; var++) {
				for (size_t dim = 0; dim < d; dim++) {
					*dFdW_ptr[InddFdW++] = 0.0;
				}}

				// *** eq 2 ***
				for (size_t var = 0; var < Nvar; var++) {
				for (size_t dim = 0; dim < d; dim++) {
					*dFdW_ptr[InddFdW++] = -(dtaudW[0][dim][var]);
				}}

				// *** eq 3 ***
				for (size_t var = 0; var < Nvar; var++) {
				for (size_t dim = 0; dim < d; dim++) {
					*dFdW_ptr[InddFdW++] = -(dtaudW[1][dim][var]);
				}}

				// *** eq 4 ***
				for (size_t var = 0; var < Nvar; var++) {
				for (size_t dim = 0; dim < d; dim++) {
					*dFdW_ptr[InddFdW++] = -(dtaudW[2][dim][var]);
				}}

				// *** eq 5 ***
				for (size_t var = 0; var < Nvar; var++) {
				for (size_t dim = 0; dim < d; dim++) {
					*dFdW_ptr[InddFdW++] = -( ( (dudW[var]*tau[dim][0]+dvdW[var]*tau[dim][1]+dwdW[var]*tau[dim][2])
					                           +(u*dtaudW[dim][0][var]+v*dtaudW[dim][1][var]+w*dtaudW[dim][2][var]) )
					                         +(mu*GAMMA/Pr*ddTsdW[dim][var]                                     ) );
				}}

				if (!DB.Const_mu) {
					InddFdW -= DMAX*Nvar;
					for (size_t var = 0; var < Nvar; var++) {
					for (size_t dim = 0; dim < d; dim++) {
						*dFdW_ptr[InddFdW++] -= dmudW[var]*GAMMA/Pr*dTs[dim];
					}}
				}

				for (size_t i = 0, iMax = Neq*Nvar*DMAX; i < iMax; i++)
					dFdW_ptr[i]++;
			}

			// ***************************************** dFdQ ***************************************** //
			if (dFdQ != NULL) {
				for (size_t dim1 = 0; dim1 < d; dim1++) {
					const double ddrhodQ[]  = { 1.0, 0.0, 0.0, 0.0, 0.0 },
					             ddrhoudQ[] = { 0.0, 1.0, 0.0, 0.0, 0.0 },
					             ddrhovdQ[] = { 0.0, 0.0, 1.0, 0.0, 0.0 },
					             ddrhowdQ[] = { 0.0, 0.0, 0.0, 1.0, 0.0 },
					             ddEdQ[]    = { 0.0, 0.0, 0.0, 0.0, 1.0 };

					double ddudQ[3][NVAR3D] = {{0.0}},
					       ddvdQ[3][NVAR3D] = {{0.0}},
					       ddwdQ[3][NVAR3D] = {{0.0}};
					for (size_t var = 0; var < Nvar; var++) {
						ddudQ[dim1][var] = rho_inv*(ddrhoudQ[var]-ddrhodQ[var]*u);
						ddvdQ[dim1][var] = rho_inv*(ddrhovdQ[var]-ddrhodQ[var]*v);
						ddwdQ[dim1][var] = rho_inv*(ddrhowdQ[var]-ddrhodQ[var]*w);
					}

					double ddivVdQ[Nvar];
					for (size_t var = 0; var < Nvar; var++)
						ddivVdQ[var] = ddudQ[0][var]+ddvdQ[1][var]+ddwdQ[2][var];

					double dtaudQ[d][d][Nvar];
					for (size_t var = 0; var < Nvar; var++) {
						dtaudQ[0][0][var] = mu*2.0*(ddudQ[0][var]-ddivVdQ[var]/3.0);
						dtaudQ[0][1][var] = mu*(ddvdQ[0][var]+ddudQ[1][var]);
						dtaudQ[0][2][var] = mu*(ddwdQ[0][var]+ddudQ[2][var]);
						dtaudQ[1][0][var] = dtaudQ[0][1][var];
						dtaudQ[1][1][var] = mu*2.0*(ddvdQ[1][var]-ddivVdQ[var]/3.0);
						dtaudQ[1][2][var] = mu*(ddwdQ[1][var]+ddvdQ[2][var]);
						dtaudQ[2][0][var] = dtaudQ[0][2][var];
						dtaudQ[2][1][var] = dtaudQ[1][2][var];
						dtaudQ[2][2][var] = mu*2.0*(ddwdQ[2][var]-ddivVdQ[var]/3.0);
					}

					if (!DB.Const_mu) {
						EXIT_UNSUPPORTED; // Ensure that mu does not depend on solution gradients for this to be OK.
					}

					double ddTsdQ[3][NVAR3D] = {{0.0}};
					for (size_t var = 0; var < Nvar; var++) {
						double const ddEoRhodQ = rho_inv2*(ddEdQ[var]*rho-E*ddrhodQ[var]),
						             ddV2dQ    = 2.0*(u*ddudQ[dim1][var]+v*ddvdQ[dim1][var]+w*ddwdQ[dim1][var]);
						ddTsdQ[dim1][var] = ddEoRhodQ - 0.5*ddV2dQ;
					}

					size_t InddFdQ = 0;
					// *** eq 1 ***
					for (size_t var = 0; var < Nvar; var++) {
					for (size_t dim = 0; dim < d; dim++) {
						*dFdQ_ptr[dim1][InddFdQ++] = 0.0;
					}}

					// *** eq 2 ***
					for (size_t var = 0; var < Nvar; var++) {
					for (size_t dim = 0; dim < d; dim++) {
						*dFdQ_ptr[dim1][InddFdQ++] = -(dtaudQ[0][dim][var]);
					}}

					// *** eq 3 ***
					for (size_t var = 0; var < Nvar; var++) {
					for (size_t dim = 0; dim < d; dim++) {
						*dFdQ_ptr[dim1][InddFdQ++] = -(dtaudQ[1][dim][var]);
					}}

					// *** eq 4 ***
					for (size_t var = 0; var < Nvar; var++) {
					for (size_t dim = 0; dim < d; dim++) {
						*dFdQ_ptr[dim1][InddFdQ++] = -(dtaudQ[2][dim][var]);
					}}

					// *** eq 5 ***
					for (size_t var = 0; var < Nvar; var++) {
					for (size_t dim = 0; dim < d; dim++) {
						*dFdQ_ptr[dim1][InddFdQ++] = -( u*dtaudQ[dim][0][var]+v*dtaudQ[dim][1][var]+w*dtaudQ[dim][2][var]
						                               +mu*GAMMA/Pr*ddTsdQ[dim][var] );
					}}

					if (!DB.Const_mu) {
						EXIT_UNSUPPORTED; // Ensure that mu does not depend on solution gradients for this to be OK.
					}

					for (size_t i = 0, iMax = Neq*Nvar*DMAX; i < iMax; i++)
						dFdQ_ptr[dim1][i]++;
				}
			}
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

			// ***************************************** F ***************************************** //
			if (F != NULL) {
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

			// ***************************************** dFdW ***************************************** //
			if (dFdW != NULL) {
				double drhodW[]      = { 1.0,                  0.0, 0.0, 0.0 },
				       drho_invdW[]  = {-rho_inv2,             0.0, 0.0, 0.0 },
				       drho_inv2dW[] = {-2.0*rho_inv2*rho_inv, 0.0, 0.0, 0.0 },

				       dudW[] = { -u*rho_inv, rho_inv, 0.0,     0.0 },
				       dvdW[] = { -v*rho_inv, 0.0,     rho_inv, 0.0 },

				       dEdW[] = { 0.0, 0.0, 0.0, 1.0 };

				double ddudW[d][Nvar], ddvdW[d][Nvar];
				for (size_t var = 0; var < Nvar; var++) {
				for (size_t dim = 0; dim < d; dim++) {
					ddudW[dim][var] = drho_invdW[var]*(drhou[dim]-drho[dim]*u)+rho_inv*(0.0-drho[dim]*dudW[var]);
					ddvdW[dim][var] = drho_invdW[var]*(drhov[dim]-drho[dim]*v)+rho_inv*(0.0-drho[dim]*dvdW[var]);
				}}

				double ddivVdW[Nvar];
				for (size_t var = 0; var < Nvar; var++)
					ddivVdW[var] = ddudW[0][var]+ddvdW[1][var];

				double dtaudW[d][d][Nvar];
				for (size_t var = 0; var < Nvar; var++) {
					dtaudW[0][0][var] = mu*2.0*(ddudW[0][var]-ddivVdW[var]/3.0);
					dtaudW[0][1][var] = mu*(ddvdW[0][var]+ddudW[1][var]);
					dtaudW[1][0][var] = dtaudW[0][1][var];
					dtaudW[1][1][var] = mu*2.0*(ddvdW[1][var]-ddivVdW[var]/3.0);
				}

				double dmudW[NVAR2D] = {0.0};
				if (!DB.Const_mu) {
					EXIT_UNSUPPORTED;

					for (size_t var = 0; var < Nvar; var++) {
						dtaudW[0][0][var] += dmudW[var]*2.0*(du[0]-divV/3.0);
						dtaudW[0][1][var] += dmudW[var]*(dv[0]+du[1]);
						dtaudW[1][0][var]  = dtaudW[0][1][var];
						dtaudW[1][1][var] += dmudW[var]*2.0*(dv[1]-divV/3.0);
					}
				}

				double ddTsdW[d][Nvar];
				for (size_t var = 0; var < Nvar; var++) {
				for (size_t dim = 0; dim < d; dim++) {
					double const ddEoRhodW = drho_inv2dW[var]*(dE[dim]*rho-E*drho[dim])
					                        +rho_inv2*((dE[dim]*drhodW[var])-(dEdW[var]*drho[dim])),
					             ddV2dW    = 2.0*( dudW[var]*du[dim]+dvdW[var]*dv[dim]+u*ddudW[dim][var]+v*ddvdW[dim][var]);
					ddTsdW[dim][var] = ddEoRhodW-0.5*ddV2dW;
				}}


				size_t InddFdW = 0;
				// *** eq 1 ***
				for (size_t var = 0; var < Nvar; var++) {
					for (size_t dim = 0; dim < d; dim++)
						*dFdW_ptr[InddFdW++] = 0.0;
					InddFdW += 1;
				}

				// *** eq 2 ***
				for (size_t var = 0; var < Nvar; var++) {
					for (size_t dim = 0; dim < d; dim++)
						*dFdW_ptr[InddFdW++] = -(dtaudW[0][dim][var]);
					InddFdW += 1;
				}

				// *** eq 3 ***
				for (size_t var = 0; var < Nvar; var++) {
					for (size_t dim = 0; dim < d; dim++)
						*dFdW_ptr[InddFdW++] = -(dtaudW[1][dim][var]);
					InddFdW += 1;
				}

				// *** eq 4 ***
				for (size_t var = 0; var < Nvar; var++) {
					for (size_t dim = 0; dim < d; dim++) {
						*dFdW_ptr[InddFdW++] = -( ( (dudW[var]*tau[dim][0]+dvdW[var]*tau[dim][1])
						                           +(u*dtaudW[dim][0][var]+v*dtaudW[dim][1][var]) )
						                         +(mu*GAMMA/Pr*ddTsdW[dim][var]                   ) );
					}
					InddFdW += 1;
				}

				if (!DB.Const_mu) {
					InddFdW -= DMAX*Nvar;
					for (size_t var = 0; var < Nvar; var++) {
						for (size_t dim = 0; dim < d; dim++)
							*dFdW_ptr[InddFdW++] -= dmudW[var]*GAMMA/Pr*dTs[dim];
						InddFdW += 1;
					}
				}

				for (size_t i = 0, iMax = Neq*Nvar*DMAX; i < iMax; i++)
					dFdW_ptr[i]++;
			}

			// ***************************************** dFdQ ***************************************** //
			if (dFdQ != NULL) {
				for (size_t dim1 = 0; dim1 < d; dim1++) {
					const double ddrhodQ[]  = { 1.0, 0.0, 0.0, 0.0 },
					             ddrhoudQ[] = { 0.0, 1.0, 0.0, 0.0 },
					             ddrhovdQ[] = { 0.0, 0.0, 1.0, 0.0 },
					             ddEdQ[]    = { 0.0, 0.0, 0.0, 1.0 };

					double ddudQ[2][NVAR2D] = {{0.0}},
					       ddvdQ[2][NVAR2D] = {{0.0}};
					for (size_t var = 0; var < Nvar; var++) {
						ddudQ[dim1][var] = rho_inv*(ddrhoudQ[var]-ddrhodQ[var]*u);
						ddvdQ[dim1][var] = rho_inv*(ddrhovdQ[var]-ddrhodQ[var]*v);
					}

					double ddivVdQ[Nvar];
					for (size_t var = 0; var < Nvar; var++)
						ddivVdQ[var] = ddudQ[0][var]+ddvdQ[1][var];

					double dtaudQ[d][d][Nvar];
					for (size_t var = 0; var < Nvar; var++) {
						dtaudQ[0][0][var] = mu*2.0*(ddudQ[0][var]-ddivVdQ[var]/3.0);
						dtaudQ[0][1][var] = mu*(ddvdQ[0][var]+ddudQ[1][var]);
						dtaudQ[1][0][var] = dtaudQ[0][1][var];
						dtaudQ[1][1][var] = mu*2.0*(ddvdQ[1][var]-ddivVdQ[var]/3.0);
					}

					if (!DB.Const_mu) {
						EXIT_UNSUPPORTED; // Ensure that mu does not depend on solution gradients for this to be OK.
					}

					double ddTsdQ[2][NVAR2D] = {{0.0}};
					for (size_t var = 0; var < Nvar; var++) {
						double const ddEoRhodQ = rho_inv2*(ddEdQ[var]*rho-E*ddrhodQ[var]),
						             ddV2dQ    = 2.0*(u*ddudQ[dim1][var]+v*ddvdQ[dim1][var]);
						ddTsdQ[dim1][var] = ddEoRhodQ - 0.5*ddV2dQ;
					}

					size_t InddFdQ = 0;
					// *** eq 1 ***
					for (size_t var = 0; var < Nvar; var++) {
						for (size_t dim = 0; dim < d; dim++)
							*dFdQ_ptr[dim1][InddFdQ++] = 0.0;
						InddFdQ += 1;
					}

					// *** eq 2 ***
					for (size_t var = 0; var < Nvar; var++) {
						for (size_t dim = 0; dim < d; dim++)
							*dFdQ_ptr[dim1][InddFdQ++] = -(dtaudQ[0][dim][var]);
						InddFdQ += 1;
					}

					// *** eq 3 ***
					for (size_t var = 0; var < Nvar; var++) {
						for (size_t dim = 0; dim < d; dim++)
							*dFdQ_ptr[dim1][InddFdQ++] = -(dtaudQ[1][dim][var]);
						InddFdQ += 1;
					}

					// *** eq 4 ***
					for (size_t var = 0; var < Nvar; var++) {
						for (size_t dim = 0; dim < d; dim++) {
							*dFdQ_ptr[dim1][InddFdQ++] = -( u*dtaudQ[dim][0][var]+v*dtaudQ[dim][1][var]
							                               +mu*GAMMA/Pr*ddTsdQ[dim][var] );
						}
						InddFdQ += 1;
					}

					if (!DB.Const_mu) {
						EXIT_UNSUPPORTED; // Ensure that mu does not depend on solution gradients for this to be OK.
					}

					for (size_t i = 0, iMax = Neq*Nvar*DMAX; i < iMax; i++)
						dFdQ_ptr[dim1][i]++;
				}
			}
		}
	}
}
