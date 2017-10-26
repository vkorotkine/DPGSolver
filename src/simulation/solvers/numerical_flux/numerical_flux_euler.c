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
#include "definitions_test_case.h"

#include "multiarray.h"

#include "const_cast.h"
#include "flux.h"
#include "flux_euler.h"
#include "numerical_flux.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void compute_Numerical_Flux_euler_lax_friedrichs
	(const struct Numerical_Flux_Input* num_flux_i, struct mutable_Numerical_Flux* num_flux)
{
	int const d   = num_flux_i->d,
	          Neq = d+2;
	const ptrdiff_t NnTotal = num_flux_i->neigh_info[0].s->extents[0];

	double const *const nL = num_flux_i->normals_l->data;

	double const *const WL = num_flux_i->neigh_info[0].s->data,
	             *const WR = num_flux_i->neigh_info[1].s->data;

	double       *const nFluxNum = num_flux->nnf->data;

	double const *rhoL_ptr  = &WL[NnTotal*0],
	             *rhouL_ptr = &WL[NnTotal*1],
	             *EL_ptr    = &WL[NnTotal*(d+1)],

	             *rhoR_ptr  = &WR[NnTotal*0],
	             *rhouR_ptr = &WR[NnTotal*1],
	             *ER_ptr    = &WR[NnTotal*(d+1)];

	double const *n_ptr = nL;

	struct Flux_Input* flux_i = malloc(sizeof *flux_i); // free
	const_cast_i(&flux_i->d,d);
	flux_i->s = constructor_move_const_Multiarray_d_d('C',2,(ptrdiff_t[]){1,Neq},false,NULL); // destructed

	struct mutable_Flux flux;
	flux.f = constructor_move_Multiarray_d_d('C',3,(ptrdiff_t[]){1,d,Neq},false,NULL); // destructed

	double *nF_ptr[Neq];
	for (int eq = 0; eq < Neq; eq++)
		nF_ptr[eq] = &nFluxNum[eq*NnTotal];

	if (d == 3) {
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
			       FLn[Neq*d];
			const_cast_d1(&flux_i->s->data,WLn);
			flux.f->data   = FLn;
			compute_Flux_euler(flux_i,&flux);

			double const *FL1_ptr = FLn,
			             *FL2_ptr = FL1_ptr+1,
			             *FL3_ptr = FL2_ptr+1;

			double WRn[] = {rhoR, rhouR, rhovR, rhowR, ER},
			       FRn[Neq*d];
			const_cast_d1(&flux_i->s->data,WRn);
			flux.f->data   = FRn;
			compute_Flux_euler(flux_i,&flux);

			double const *FR1_ptr = FRn,
			             *FR2_ptr = FR1_ptr+1,
			             *FR3_ptr = FR2_ptr+1;

			for (int eq = 0; eq < Neq; eq++) {
				*nF_ptr[eq]++ = 0.5*( n1*((*FL1_ptr)+(*FR1_ptr))+n2*((*FL2_ptr)+(*FR2_ptr))+n3*((*FL3_ptr)+(*FR3_ptr))
				                     + maxV*(WLn[eq]-WRn[eq]));

				FL1_ptr += d;
				FL2_ptr += d;
				FL3_ptr += d;
				FR1_ptr += d;
				FR2_ptr += d;
				FR3_ptr += d;
			}
		}
	} else if (d == 2) {
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
			       FLn[Neq*d];
			const_cast_d1(&flux_i->s->data,WLn);
			flux.f->data   = FLn;
			compute_Flux_euler(flux_i,&flux);

			double const *FL1_ptr = FLn,
			             *FL2_ptr = FL1_ptr+1;

			double WRn[] = {rhoR, rhouR, rhovR, ER},
			       FRn[Neq*d];
			const_cast_d1(&flux_i->s->data,WRn);
			flux.f->data   = FRn;
			compute_Flux_euler(flux_i,&flux);

			double const *FR1_ptr = FRn,
			             *FR2_ptr = FR1_ptr+1;

			for (int eq = 0; eq < Neq; eq++) {
				*nF_ptr[eq]++ = 0.5*( n1*((*FL1_ptr)+(*FR1_ptr))+n2*((*FL2_ptr)+(*FR2_ptr)) + maxV*(WLn[eq]-WRn[eq]));

				FL1_ptr += d;
				FL2_ptr += d;
				FR1_ptr += d;
				FR2_ptr += d;
			}
		}
	} else if (d == 1) {
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
			       FLn[Neq*d];
			const_cast_d1(&flux_i->s->data,WLn);
			flux.f->data   = FLn;
			compute_Flux_euler(flux_i,&flux);

			double const *FL1_ptr = FLn;

			double WRn[] = {rhoR, rhouR, ER},
			       FRn[Neq*d];
			const_cast_d1(&flux_i->s->data,WRn);
			flux.f->data   = FRn;
			compute_Flux_euler(flux_i,&flux);

			double const *FR1_ptr = FRn;

			for (int eq = 0; eq < Neq; eq++) {
				*nF_ptr[eq]++ = 0.5*(n1*((*FL1_ptr)+(*FR1_ptr)) + maxV*(WLn[eq]-WRn[eq]));

				FL1_ptr += d;
				FR1_ptr += d;
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

	int const d   = num_flux_i->d,
	          Neq = d+2;
	const ptrdiff_t NnTotal = num_flux_i->neigh_info[0].s->extents[0];

	double const *const nL = num_flux_i->normals_l->data;

	double const *const WL = num_flux_i->neigh_info[0].s->data,
	             *const WR = num_flux_i->neigh_info[1].s->data;

	double       *const nFluxNum = num_flux->nnf->data;

	double const *rhoL_ptr  = &WL[NnTotal*0],
	             *rhouL_ptr = &WL[NnTotal*1],
	             *EL_ptr    = &WL[NnTotal*(d+1)],

	             *rhoR_ptr  = &WR[NnTotal*0],
	             *rhouR_ptr = &WR[NnTotal*1],
	             *ER_ptr    = &WR[NnTotal*(d+1)];

	double const *n_ptr = nL;

	double *nF_ptr[Neq];
	for (int eq = 0; eq < Neq; eq++)
		nF_ptr[eq] = &nFluxNum[eq*NnTotal];

	if (d == 3) {
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
	} else if (d == 2) {
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
	} else if (d == 1) {
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

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
