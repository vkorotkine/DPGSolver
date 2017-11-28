/* {{{
This file is part of DPGSolver.

DPGSolver is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or any later version.

DPGSolver is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along with DPGSolver.  If not, see
<http://www.gnu.org/licenses/>.
}}} */
/** \file
 *  \todo Clean-up.
 */

#include "numerical_flux_euler.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gsl/gsl_math.h"

#include "macros.h"
#include "definitions_core.h"
#include "definitions_test_case.h"

#include "multiarray.h"

#include "const_cast.h"
#include "flux.h"
#include "flux_euler.h"
#include "numerical_flux.h"

// Static function declarations ************************************************************************************* //

#define NEQ  NEQ_EULER  ///< Number of equations.
#define NVAR NVAR_EULER ///< Number of variables.

// Interface functions ********************************************************************************************** //

void compute_Numerical_Flux_euler_lax_friedrichs
	(const struct Numerical_Flux_Input* num_flux_i, struct mutable_Numerical_Flux* num_flux)
{
	const ptrdiff_t NnTotal = num_flux_i->bv_l.s->extents[0];

	double const *const nL = num_flux_i->bv_l.normals->data;

	double const *const WL = num_flux_i->bv_l.s->data,
	             *const WR = num_flux_i->bv_r.s->data;

	double       *const nFluxNum = num_flux->nnf->data;

	double const *rhoL_ptr  = &WL[NnTotal*0],
	             *rhouL_ptr = &WL[NnTotal*1],
	             *EL_ptr    = &WL[NnTotal*(DIM+1)],

	             *rhoR_ptr  = &WR[NnTotal*0],
	             *rhouR_ptr = &WR[NnTotal*1],
	             *ER_ptr    = &WR[NnTotal*(DIM+1)];

	double const *n_ptr = nL;

	struct Flux_Input* flux_i = malloc(sizeof *flux_i); // free
	flux_i->s = constructor_move_const_Multiarray_d_d('C',2,(ptrdiff_t[]){1,NEQ},false,NULL); // destructed

	struct mutable_Flux flux;
	flux.f = constructor_move_Multiarray_d_d('C',3,(ptrdiff_t[]){1,DIM,NEQ},false,NULL); // destructed

	double *nF_ptr[NEQ];
	for (int eq = 0; eq < NEQ; eq++)
		nF_ptr[eq] = &nFluxNum[eq*NnTotal];

	if (DIM == 3) {
		double const *rhovL_ptr = &WL[NnTotal*2],
		             *rhowL_ptr = &WL[NnTotal*3],

		             *rhovR_ptr = &WR[NnTotal*2],
		             *rhowR_ptr = &WR[NnTotal*3];

		for (ptrdiff_t n = 0; n < NnTotal; n++) {
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

			double maxV;
			if (maxlL > maxlR)
				maxV = maxlL;
			else
				maxV = maxlR;

			double const n1 = *n_ptr++,
			             n2 = *n_ptr++,
			             n3 = *n_ptr++;

			double WLn[] = {rhoL, rhouL, rhovL, rhowL, EL},
			       FLn[NEQ*DIM];
			const_cast_d1(&flux_i->s->data,WLn);
			flux.f->data   = FLn;
			compute_Flux_euler(flux_i,&flux);

			double const *FL1_ptr = FLn,
			             *FL2_ptr = FL1_ptr+1,
			             *FL3_ptr = FL2_ptr+1;

			double WRn[] = {rhoR, rhouR, rhovR, rhowR, ER},
			       FRn[NEQ*DIM];
			const_cast_d1(&flux_i->s->data,WRn);
			flux.f->data   = FRn;
			compute_Flux_euler(flux_i,&flux);

			double const *FR1_ptr = FRn,
			             *FR2_ptr = FR1_ptr+1,
			             *FR3_ptr = FR2_ptr+1;

			for (int eq = 0; eq < NEQ; eq++) {
				*nF_ptr[eq]++ = 0.5*( n1*((*FL1_ptr)+(*FR1_ptr))+n2*((*FL2_ptr)+(*FR2_ptr))+n3*((*FL3_ptr)+(*FR3_ptr))
				                     + maxV*(WLn[eq]-WRn[eq]));

				FL1_ptr += DIM;
				FL2_ptr += DIM;
				FL3_ptr += DIM;
				FR1_ptr += DIM;
				FR2_ptr += DIM;
				FR3_ptr += DIM;
			}
		}
	} else if (DIM == 2) {
		double const *rhovL_ptr = &WL[NnTotal*2],

		             *rhovR_ptr = &WR[NnTotal*2];

		for (ptrdiff_t n = 0; n < NnTotal; n++) {
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

			double maxV;
			if (maxlL > maxlR)
				maxV = maxlL;
			else
				maxV = maxlR;

			double const n1 = *n_ptr++,
			             n2 = *n_ptr++;

			double WLn[] = {rhoL, rhouL, rhovL, EL},
			       FLn[NEQ*DIM];
			const_cast_d1(&flux_i->s->data,WLn);
			flux.f->data   = FLn;
			compute_Flux_euler(flux_i,&flux);

			double const *FL1_ptr = FLn,
			             *FL2_ptr = FL1_ptr+1;

			double WRn[] = {rhoR, rhouR, rhovR, ER},
			       FRn[NEQ*DIM];
			const_cast_d1(&flux_i->s->data,WRn);
			flux.f->data   = FRn;
			compute_Flux_euler(flux_i,&flux);

			double const *FR1_ptr = FRn,
			             *FR2_ptr = FR1_ptr+1;

			for (int eq = 0; eq < NEQ; eq++) {
				*nF_ptr[eq]++ = 0.5*( n1*((*FL1_ptr)+(*FR1_ptr))+n2*((*FL2_ptr)+(*FR2_ptr)) + maxV*(WLn[eq]-WRn[eq]));

				FL1_ptr += DIM;
				FL2_ptr += DIM;
				FR1_ptr += DIM;
				FR2_ptr += DIM;
			}
		}
	} else if (DIM == 1) {
		for (ptrdiff_t n = 0; n < NnTotal; n++) {
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

			double maxV;
			if (maxlL > maxlR)
				maxV = maxlL;
			else
				maxV = maxlR;

			double const n1 = *n_ptr++;

			double WLn[] = {rhoL, rhouL, EL},
			       FLn[NEQ*DIM];
			const_cast_d1(&flux_i->s->data,WLn);
			flux.f->data   = FLn;
			compute_Flux_euler(flux_i,&flux);

			double const *FL1_ptr = FLn;

			double WRn[] = {rhoR, rhouR, ER},
			       FRn[NEQ*DIM];
			const_cast_d1(&flux_i->s->data,WRn);
			flux.f->data   = FRn;
			compute_Flux_euler(flux_i,&flux);

			double const *FR1_ptr = FRn;

			for (int eq = 0; eq < NEQ; eq++) {
				*nF_ptr[eq]++ = 0.5*(n1*((*FL1_ptr)+(*FR1_ptr)) + maxV*(WLn[eq]-WRn[eq]));

				FL1_ptr += DIM;
				FR1_ptr += DIM;
			}
		}
	}
	destructor_const_Multiarray_d(flux_i->s);
	destructor_Multiarray_d(flux.f);
	free(flux_i);
}

void compute_Numerical_Flux_euler_roe_pike
	(const struct Numerical_Flux_Input* num_flux_i, struct mutable_Numerical_Flux* num_flux)
{
	/// The simple entropy fix is taken from (eq. (35), \cite Qu2015).

	const ptrdiff_t NnTotal = num_flux_i->bv_l.s->extents[0];

	double const *const nL = num_flux_i->bv_l.normals->data;

	double const *const WL = num_flux_i->bv_l.s->data,
	             *const WR = num_flux_i->bv_r.s->data;

	double       *const nFluxNum = num_flux->nnf->data;

	double const *rhoL_ptr  = &WL[NnTotal*0],
	             *rhouL_ptr = &WL[NnTotal*1],
	             *EL_ptr    = &WL[NnTotal*(DIM+1)],

	             *rhoR_ptr  = &WR[NnTotal*0],
	             *rhouR_ptr = &WR[NnTotal*1],
	             *ER_ptr    = &WR[NnTotal*(DIM+1)];

	double const *n_ptr = nL;

	double *nF_ptr[NEQ];
	for (int eq = 0; eq < NEQ; eq++)
		nF_ptr[eq] = &nFluxNum[eq*NnTotal];

	if (DIM == 3) {
		double const *rhovL_ptr = &WL[NnTotal*2],
		             *rhowL_ptr = &WL[NnTotal*3],

		             *rhovR_ptr = &WR[NnTotal*2],
		             *rhowR_ptr = &WR[NnTotal*3];

		for (ptrdiff_t n = 0; n < NnTotal; n++) {
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
			             c   = sqrt(c2);

			// Compute eigenvalues (with entropy fix)
			double const l1   = GSL_MIN(fabs(Vn-c),fabs(VnL-c)),
			             l234 = fabs(Vn),
			             l5   = GSL_MAX(fabs(Vn+c),fabs(VnR+c));

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

			// Compute contribution of normal flux components (multiplied by 0.5 below)
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
			int IndnF = 0;
			*nF_ptr[IndnF++]++ = 0.5*(nF1 - dis1);
			*nF_ptr[IndnF++]++ = 0.5*(nF2 - dis2);
			*nF_ptr[IndnF++]++ = 0.5*(nF3 - dis3);
			*nF_ptr[IndnF++]++ = 0.5*(nF4 - dis4);
			*nF_ptr[IndnF++]++ = 0.5*(nF5 - dis5);
		}
	} else if (DIM == 2) {
		double const *rhovL_ptr = &WL[NnTotal*2],

		             *rhovR_ptr = &WR[NnTotal*2];

		for (ptrdiff_t n = 0; n < NnTotal; n++) {
			double const n1 = *n_ptr++,
			             n2 = *n_ptr++;

			// Left VOLUME
			double const rhoL  = *rhoL_ptr++,
			             rhouL = *rhouL_ptr++,
			             rhovL = *rhovL_ptr++,
			             EL    = *EL_ptr++,

			             rhoL_inv  = 1.0/rhoL,
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
			             c   = sqrt(c2);

			// Compute eigenvalues (with entropy fix)
			double const l1   = GSL_MIN(fabs(Vn-c),fabs(VnL-c)),
			             l234 = fabs(Vn),
			             l5   = GSL_MAX(fabs(Vn+c),fabs(VnR+c));

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

			// Compute contribution of normal flux components (multiplied by 0.5 below)
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
			int IndnF = 0;
			*nF_ptr[IndnF++]++ = 0.5*(nF1 - dis1);
			*nF_ptr[IndnF++]++ = 0.5*(nF2 - dis2);
			*nF_ptr[IndnF++]++ = 0.5*(nF3 - dis3);
			*nF_ptr[IndnF++]++ = 0.5*(nF5 - dis5);
		}
	} else if (DIM == 1) {
		for (ptrdiff_t n = 0; n < NnTotal; n++) {
			double const n1 = *n_ptr++;

			// Left VOLUME
			double const rhoL  = *rhoL_ptr++,
			             rhouL = *rhouL_ptr++,
			             EL    = *EL_ptr++,

			             rhoL_inv  = 1.0/rhoL,
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
			             c   = sqrt(c2);

			// Compute eigenvalues (with entropy fix)
			double const l1   = GSL_MIN(fabs(Vn-c),fabs(VnL-c)),
			             l234 = fabs(Vn),
			             l5   = GSL_MAX(fabs(Vn+c),fabs(VnR+c));

			double const drho  = rhoR-rhoL,
			             drhou = rhoR*uR-rhoL*uL,
			             dE    = ER-EL,
			             dp    = pR-pL,
			             dVn   = VnR-VnL,

			             lc1 = 0.5*(l5+l1) - l234,
			             lc2 = 0.5*(l5-l1),

			             disInter1 = lc1*dp/(c*c) + lc2*rho*dVn/c,
			             disInter2 = lc1*rho*dVn  + lc2*dp/c;

			// Compute contribution of normal flux components (multiplied by 0.5 below)
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
			int IndnF = 0;
			*nF_ptr[IndnF++]++ = 0.5*(nF1 - dis1);
			*nF_ptr[IndnF++]++ = 0.5*(nF2 - dis2);
			*nF_ptr[IndnF++]++ = 0.5*(nF5 - dis5);
		}
	}
}

void compute_Numerical_Flux_euler_roe_pike_jacobian
	(const struct Numerical_Flux_Input* num_flux_i, struct mutable_Numerical_Flux* num_flux)
{
	/// The simple entropy fix is taken from (eq. (35), \cite Qu2015).

	const ptrdiff_t NnTotal = num_flux_i->bv_l.s->extents[0];

	double const *const nL = num_flux_i->bv_l.normals->data;

	double const *const WL = num_flux_i->bv_l.s->data,
	             *const WR = num_flux_i->bv_r.s->data;

	double *const nFluxNum = (num_flux->nnf ? num_flux->nnf->data : NULL ),
	       *const dnFdWL   = num_flux->neigh_info[0].dnnf_ds->data,
	       *const dnFdWR   = num_flux->neigh_info[1].dnnf_ds->data;

	double const *rhoL_ptr  = &WL[NnTotal*0],
	             *rhouL_ptr = &WL[NnTotal*1],
	             *EL_ptr    = &WL[NnTotal*(DIM+1)],

	             *rhoR_ptr  = &WR[NnTotal*0],
	             *rhouR_ptr = &WR[NnTotal*1],
	             *ER_ptr    = &WR[NnTotal*(DIM+1)];

	double const *n_ptr = nL;

	double *nF_ptr[NEQ];
	if (nFluxNum != NULL) {
		for (int eq = 0; eq < NEQ; eq++)
			nF_ptr[eq] = &nFluxNum[eq*NnTotal];
	}

	double *dnFdWL_ptr[NEQ*NEQ];
	for (int eq = 0; eq < NEQ; eq++) {
	for (int var = 0; var < NVAR; var++) {
		dnFdWL_ptr[eq*NVAR+var] = &dnFdWL[(eq*NVAR+var)*NnTotal];
	}}

	double *dnFdWR_ptr[NEQ*NEQ];
	if (dnFdWR != NULL) {
		for (int eq = 0; eq < NEQ; eq++) {
		for (int var = 0; var < NVAR; var++) {
			dnFdWR_ptr[eq*NVAR+var] = &dnFdWR[(eq*NVAR+var)*NnTotal];
		}}
	}

	if (DIM == 3) {
		double const *rhovL_ptr = &WL[NnTotal*2],
		             *rhowL_ptr = &WL[NnTotal*3],

		             *rhovR_ptr = &WR[NnTotal*2],
		             *rhowR_ptr = &WR[NnTotal*3];

		for (ptrdiff_t n = 0; n < NnTotal; n++) {
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
				int IndnF = 0;
				*nF_ptr[IndnF++]++ = 0.5*(nF1 - dis1);
				*nF_ptr[IndnF++]++ = 0.5*(nF2 - dis2);
				*nF_ptr[IndnF++]++ = 0.5*(nF3 - dis3);
				*nF_ptr[IndnF++]++ = 0.5*(nF4 - dis4);
				*nF_ptr[IndnF++]++ = 0.5*(nF5 - dis5);
			}

			double dnF1dW[NEQ],  dnF2dW[NEQ],  dnF3dW[NEQ],  dnF4dW[NEQ],  dnF5dW[NEQ],
			       ddis1dW[NEQ], ddis2dW[NEQ], ddis3dW[NEQ], ddis4dW[NEQ], ddis5dW[NEQ];

			// Flux term
			double const duLdW[]   = {-uL*rhoL_inv,  rhoL_inv,  0.0,       0.0,      0.0},
			             dvLdW[]   = {-vL*rhoL_inv,  0.0,       rhoL_inv,  0.0,      0.0},
			             dwLdW[]   = {-wL*rhoL_inv,  0.0,       0.0,       rhoL_inv, 0.0},
			             drhoLdW[] = { 1.0,          0.0,       0.0,       0.0,      0.0},
			             dELdW[]   = { 0.0,          0.0,       0.0,       0.0,      1.0},
			             dpLdW[]   = { GM1*0.5*V2L, -GM1*uL,   -GM1*vL,   -GM1*wL,   GM1};

			double const rhoVn = rhoL*VnL;
			double dVnLdW[NVAR];
			for (int var = 0; var < NVAR; var++) {
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

			for (int var = 0; var < NVAR; var++) {
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

			int InddnFdW = 0;
			for (int var = 0; var < NVAR; var++) *dnFdWL_ptr[InddnFdW++]++ = 0.5*(dnF1dW[var]-ddis1dW[var]);
			for (int var = 0; var < NVAR; var++) *dnFdWL_ptr[InddnFdW++]++ = 0.5*(dnF2dW[var]-ddis2dW[var]);
			for (int var = 0; var < NVAR; var++) *dnFdWL_ptr[InddnFdW++]++ = 0.5*(dnF3dW[var]-ddis3dW[var]);
			for (int var = 0; var < NVAR; var++) *dnFdWL_ptr[InddnFdW++]++ = 0.5*(dnF4dW[var]-ddis4dW[var]);
			for (int var = 0; var < NVAR; var++) *dnFdWL_ptr[InddnFdW++]++ = 0.5*(dnF5dW[var]-ddis5dW[var]);

			if (dnFdWR != NULL) {
				// Flux term
				double const duRdW[]   = {-uR*rhoR_inv,  rhoR_inv,  0.0,       0.0,      0.0},
				             dvRdW[]   = {-vR*rhoR_inv,  0.0,       rhoR_inv,  0.0,      0.0},
				             dwRdW[]   = {-wR*rhoR_inv,  0.0,       0.0,       rhoR_inv, 0.0},
				             drhoRdW[] = { 1.0,          0.0,       0.0,       0.0,      0.0},
				             dERdW[]   = { 0.0,          0.0,       0.0,       0.0,      1.0},
				             dpRdW[]   = { GM1*0.5*V2R, -GM1*uR,   -GM1*vR,   -GM1*wR,   GM1};

				double const rhoVn = rhoR*VnR;
				double dVnRdW[NVAR];
				for (int var = 0; var < NVAR; var++) {
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

				for (int var = 0; var < NVAR; var++) {
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

				int InddnFdW = 0;
				for (int var = 0; var < NVAR; var++) *dnFdWR_ptr[InddnFdW++]++ = 0.5*(dnF1dW[var]-ddis1dW[var]);
				for (int var = 0; var < NVAR; var++) *dnFdWR_ptr[InddnFdW++]++ = 0.5*(dnF2dW[var]-ddis2dW[var]);
				for (int var = 0; var < NVAR; var++) *dnFdWR_ptr[InddnFdW++]++ = 0.5*(dnF3dW[var]-ddis3dW[var]);
				for (int var = 0; var < NVAR; var++) *dnFdWR_ptr[InddnFdW++]++ = 0.5*(dnF4dW[var]-ddis4dW[var]);
				for (int var = 0; var < NVAR; var++) *dnFdWR_ptr[InddnFdW++]++ = 0.5*(dnF5dW[var]-ddis5dW[var]);
			}
		}
	} else if (DIM == 2) {
		double const *rhovL_ptr = &WL[NnTotal*2],

		             *rhovR_ptr = &WR[NnTotal*2];

		for (ptrdiff_t n = 0; n < NnTotal; n++) {
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
printf("nfe: %f\n",l1L);
printf("%f\n",l1);
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
				int IndnF = 0;
				*nF_ptr[IndnF++]++ = 0.5*(nF1 - dis1);
				*nF_ptr[IndnF++]++ = 0.5*(nF2 - dis2);
				*nF_ptr[IndnF++]++ = 0.5*(nF3 - dis3);
				*nF_ptr[IndnF++]++ = 0.5*(nF5 - dis5);
			}

			double dnF1dW[NEQ],  dnF2dW[NEQ],  dnF3dW[NEQ],  dnF5dW[NEQ],
			       ddis1dW[NEQ], ddis2dW[NEQ], ddis3dW[NEQ], ddis5dW[NEQ];

			// Flux term
			double const duLdW[]   = {-uL*rhoL_inv,  rhoL_inv,  0.0,      0.0},
			             dvLdW[]   = {-vL*rhoL_inv,  0.0,       rhoL_inv, 0.0},
			             drhoLdW[] = { 1.0,          0.0,       0.0,      0.0},
			             dELdW[]   = { 0.0,          0.0,       0.0,      1.0},
			             dpLdW[]   = { GM1*0.5*V2L, -GM1*uL,   -GM1*vL,   GM1};

			double const rhoVn = rhoL*VnL;
			double dVnLdW[NVAR];
			for (int var = 0; var < NVAR; var++) {
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

			for (int var = 0; var < NVAR; var++) {
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

			int InddnFdW = 0;
			for (int var = 0; var < NVAR; var++) *dnFdWL_ptr[InddnFdW++]++ = 0.5*(dnF1dW[var]-ddis1dW[var]);
			for (int var = 0; var < NVAR; var++) *dnFdWL_ptr[InddnFdW++]++ = 0.5*(dnF2dW[var]-ddis2dW[var]);
			for (int var = 0; var < NVAR; var++) *dnFdWL_ptr[InddnFdW++]++ = 0.5*(dnF3dW[var]-ddis3dW[var]);
			for (int var = 0; var < NVAR; var++) *dnFdWL_ptr[InddnFdW++]++ = 0.5*(dnF5dW[var]-ddis5dW[var]);

			if (dnFdWR != NULL) {
				// Flux term
				double const duRdW[]   = {-uR*rhoR_inv,  rhoR_inv,  0.0,      0.0},
				             dvRdW[]   = {-vR*rhoR_inv,  0.0,       rhoR_inv, 0.0},
				             drhoRdW[] = { 1.0,          0.0,       0.0,      0.0},
				             dERdW[]   = { 0.0,          0.0,       0.0,      1.0},
				             dpRdW[]   = { GM1*0.5*V2R, -GM1*uR,   -GM1*vR,   GM1};

				double const rhoVn = rhoR*VnR;
				double dVnRdW[NVAR];
				for (int var = 0; var < NVAR; var++) {
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

				for (int var = 0; var < NVAR; var++) {
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

				int InddnFdW = 0;
				for (int var = 0; var < NVAR; var++) *dnFdWR_ptr[InddnFdW++]++ = 0.5*(dnF1dW[var]-ddis1dW[var]);
				for (int var = 0; var < NVAR; var++) *dnFdWR_ptr[InddnFdW++]++ = 0.5*(dnF2dW[var]-ddis2dW[var]);
				for (int var = 0; var < NVAR; var++) *dnFdWR_ptr[InddnFdW++]++ = 0.5*(dnF3dW[var]-ddis3dW[var]);
				for (int var = 0; var < NVAR; var++) *dnFdWR_ptr[InddnFdW++]++ = 0.5*(dnF5dW[var]-ddis5dW[var]);
			}
		}
	} else if (DIM == 1) {
		for (ptrdiff_t n = 0; n < NnTotal; n++) {
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
				int IndnF = 0;
				*nF_ptr[IndnF++]++ = 0.5*(nF1 - dis1);
				*nF_ptr[IndnF++]++ = 0.5*(nF2 - dis2);
				*nF_ptr[IndnF++]++ = 0.5*(nF5 - dis5);
			}

			double dnF1dW[NEQ],  dnF2dW[NEQ],  dnF5dW[NEQ],
			       ddis1dW[NEQ], ddis2dW[NEQ], ddis5dW[NEQ];

			// Flux term
			double const duLdW[]   = {-uL*rhoL_inv,  rhoL_inv, 0.0},
			             drhoLdW[] = { 1.0,          0.0,      0.0},
			             dELdW[]   = { 0.0,          0.0,      1.0},
			             dpLdW[]   = { GM1*0.5*V2L, -GM1*uL,   GM1};

			double const rhoVn = rhoL*VnL;
			double dVnLdW[NVAR];
			for (int var = 0; var < NVAR; var++) {
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

			for (int var = 0; var < NVAR; var++) {
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

			int InddnFdW = 0;
			for (int var = 0; var < NVAR; var++) *dnFdWL_ptr[InddnFdW++]++ = 0.5*(dnF1dW[var]-ddis1dW[var]);
			for (int var = 0; var < NVAR; var++) *dnFdWL_ptr[InddnFdW++]++ = 0.5*(dnF2dW[var]-ddis2dW[var]);
			for (int var = 0; var < NVAR; var++) *dnFdWL_ptr[InddnFdW++]++ = 0.5*(dnF5dW[var]-ddis5dW[var]);

			if (dnFdWR != NULL) {
				// Flux term
				double const duRdW[]   = {-uR*rhoR_inv,  rhoR_inv, 0.0},
				             drhoRdW[] = { 1.0,          0.0,      0.0},
				             dERdW[]   = { 0.0,          0.0,      1.0},
				             dpRdW[]   = { GM1*0.5*V2R, -GM1*uR,   GM1};

				double const rhoVn = rhoR*VnR;
				double dVnRdW[NVAR];
				for (int var = 0; var < NVAR; var++) {
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

				for (int var = 0; var < NVAR; var++) {
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

				int InddnFdW = 0;
				for (int var = 0; var < NVAR; var++) *dnFdWR_ptr[InddnFdW++]++ = 0.5*(dnF1dW[var]-ddis1dW[var]);
				for (int var = 0; var < NVAR; var++) *dnFdWR_ptr[InddnFdW++]++ = 0.5*(dnF2dW[var]-ddis2dW[var]);
				for (int var = 0; var < NVAR; var++) *dnFdWR_ptr[InddnFdW++]++ = 0.5*(dnF5dW[var]-ddis5dW[var]);
			}
		}
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
