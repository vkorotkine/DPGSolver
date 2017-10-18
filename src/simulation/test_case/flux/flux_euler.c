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

#include "flux_euler.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"
#include "definitions_core.h"
#include "definitions_test_case.h"

#include "multiarray.h"

#include "flux.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void compute_Flux_euler (const struct Flux_Input* flux_i, struct mutable_Flux* flux)
{
	int const d   = flux_i->d,
	          Neq = flux_i->n_eq;
	const ptrdiff_t NnTotal = flux_i->s->extents[0];

	double const *const W = flux_i->s->data;
	double       *const F = flux->f->data;

	double const *rho_ptr  = &W[NnTotal*0],
	             *rhou_ptr = &W[NnTotal*1],
	             *E_ptr    = &W[NnTotal*(d+1)];

	double *F_ptr[d*Neq];
	for (int eq = 0; eq < Neq; eq++)  {
	for (int dim = 0; dim < d; dim++) {
		F_ptr[eq*d+dim] = &F[(eq*d+dim)*NnTotal];
	}}

	if (d == 3) {
		double const *rhov_ptr = &W[NnTotal*2],
		             *rhow_ptr = &W[NnTotal*3];

		for (ptrdiff_t n = 0; n < NnTotal; n++) {
			double const rho  = *rho_ptr++,
			             rhou = *rhou_ptr++,
			             rhov = *rhov_ptr++,
			             rhow = *rhow_ptr++,
			             E    = *E_ptr++,

			             u   = rhou/rho,
			             v   = rhov/rho,
			             w   = rhow/rho,

			             p = GM1*(E-0.5*rho*(u*u+v*v+w*w));

			ptrdiff_t IndF = 0;
			// eq 1
			*F_ptr[IndF++]++ += rhou;
			*F_ptr[IndF++]++ += rhov;
			*F_ptr[IndF++]++ += rhow;

			// eq 2
			*F_ptr[IndF++]++ += rhou*u + p;
			*F_ptr[IndF++]++ += rhou*v;
			*F_ptr[IndF++]++ += rhou*w;

			// eq 3
			*F_ptr[IndF++]++ += rhov*u;
			*F_ptr[IndF++]++ += rhov*v + p;
			*F_ptr[IndF++]++ += rhov*w;

			// eq 4
			*F_ptr[IndF++]++ += rhow*u;
			*F_ptr[IndF++]++ += rhow*v;
			*F_ptr[IndF++]++ += rhow*w + p;

			// eq 5
			*F_ptr[IndF++]++ += (E+p)*u;
			*F_ptr[IndF++]++ += (E+p)*v;
			*F_ptr[IndF++]++ += (E+p)*w;
		}
	} else if (d == 2) {
		double const *rhov_ptr = &W[NnTotal*2];

		for (ptrdiff_t n = 0; n < NnTotal; n++) {
			double const rho  = *rho_ptr++,
			             rhou = *rhou_ptr++,
			             rhov = *rhov_ptr++,
			             E    = *E_ptr++,

			             u   = rhou/rho,
			             v   = rhov/rho,

			             p = GM1*(E-0.5*rho*(u*u+v*v));

			ptrdiff_t IndF = 0;
			// eq 1
			*F_ptr[IndF++]++ += rhou;
			*F_ptr[IndF++]++ += rhov;

			// eq 2
			*F_ptr[IndF++]++ += rhou*u + p;
			*F_ptr[IndF++]++ += rhou*v;

			// eq 3
			*F_ptr[IndF++]++ += rhov*u;
			*F_ptr[IndF++]++ += rhov*v + p;

			// eq 4
			*F_ptr[IndF++]++ += (E+p)*u;
			*F_ptr[IndF++]++ += (E+p)*v;
		}
	} else if (d == 1) {
		for (ptrdiff_t n = 0; n < NnTotal; n++) {
			double const rho  = *rho_ptr++,
			             rhou = *rhou_ptr++,
			             E    = *E_ptr++,

			             u   = rhou/rho,

			             p = GM1*(E-0.5*rho*(u*u));

			ptrdiff_t IndF = 0;
			// eq 1
			*F_ptr[IndF++]++ += rhou;

			// eq 2
			*F_ptr[IndF++]++ += rhou*u + p;

			// eq 3
			*F_ptr[IndF++]++ += (E+p)*u;
		}
	}
}

void compute_Flux_euler_jacobian (const struct Flux_Input* flux_i, struct mutable_Flux* flux)
{
	int const d   = flux_i->d,
	          Neq = flux_i->n_eq;
	const ptrdiff_t NnTotal = flux_i->s->extents[0];

	double const *const W    = flux_i->s->data;
	double       *const F    = flux->f->data;
	double       *const dFdW = flux->df_ds->data;

	// Standard datatypes
	int i, n, eq, var, dim, iMax, Nvar, InddFdW;
	double       rho, u, v, w, u2, uv, uw, v2, vw, w2, V2, E, p, H, alpha, beta, *dFdW_ptr[DMAX*Neq*Neq];
	const double *rho_ptr, *rhou_ptr, *rhov_ptr, *rhow_ptr, *E_ptr;

	Nvar    = Neq;

	rho_ptr  = &W[NnTotal*0];
	rhou_ptr = &W[NnTotal*1];
	E_ptr    = &W[NnTotal*(d+1)];

	// Store pointers to the arrays that the data will be written into.
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
				ptrdiff_t IndF = 0;
				// eq 1
				*F_ptr[IndF++] += rhou;
				*F_ptr[IndF++] += rhov;
				*F_ptr[IndF++] += rhow;

				// eq 2
				*F_ptr[IndF++] += rhou*u + p;
				*F_ptr[IndF++] += rhou*v;
				*F_ptr[IndF++] += rhou*w;

				// eq 3
				*F_ptr[IndF++] += rhov*u;
				*F_ptr[IndF++] += rhov*v + p;
				*F_ptr[IndF++] += rhov*w;

				// eq 4
				*F_ptr[IndF++] += rhow*u;
				*F_ptr[IndF++] += rhow*v;
				*F_ptr[IndF++] += rhow*w + p;

				// eq 5
				*F_ptr[IndF++] += (E+p)*u;
				*F_ptr[IndF++] += (E+p)*v;
				*F_ptr[IndF++] += (E+p)*w;

				for (i = 0, iMax = Neq*DMAX; i < iMax; i++)
					F_ptr[i]++;
			}

			InddFdW = 0;
			// *** eq 1 ***
			// var 1
			*dFdW_ptr[InddFdW++] +=  0.0;
			*dFdW_ptr[InddFdW++] +=  0.0;
			*dFdW_ptr[InddFdW++] +=  0.0;

			// var 2
			*dFdW_ptr[InddFdW++] +=  1.0;
			*dFdW_ptr[InddFdW++] +=  0.0;
			*dFdW_ptr[InddFdW++] +=  0.0;

			// var 3
			*dFdW_ptr[InddFdW++] +=  0.0;
			*dFdW_ptr[InddFdW++] +=  1.0;
			*dFdW_ptr[InddFdW++] +=  0.0;

			// var 4
			*dFdW_ptr[InddFdW++] +=  0.0;
			*dFdW_ptr[InddFdW++] +=  0.0;
			*dFdW_ptr[InddFdW++] +=  1.0;

			// var 5
			*dFdW_ptr[InddFdW++] +=  0.0;
			*dFdW_ptr[InddFdW++] +=  0.0;
			*dFdW_ptr[InddFdW++] +=  0.0;

			// *** eq 2 ***
			// var 1
			*dFdW_ptr[InddFdW++] += -u2+alpha;
			*dFdW_ptr[InddFdW++] += -uv;
			*dFdW_ptr[InddFdW++] += -uw;

			// var 2
			*dFdW_ptr[InddFdW++] += -GM3*u;
			*dFdW_ptr[InddFdW++] +=  v;
			*dFdW_ptr[InddFdW++] +=  w;

			// var 3
			*dFdW_ptr[InddFdW++] += -GM1*v;
			*dFdW_ptr[InddFdW++] +=  u;
			*dFdW_ptr[InddFdW++] +=  0.0;

			// var 4
			*dFdW_ptr[InddFdW++] += -GM1*w;
			*dFdW_ptr[InddFdW++] +=  0.0;
			*dFdW_ptr[InddFdW++] +=  u;

			// var 5
			*dFdW_ptr[InddFdW++] +=  GM1;
			*dFdW_ptr[InddFdW++] +=  0.0;
			*dFdW_ptr[InddFdW++] +=  0.0;

			// *** eq 3 ***
			// var 1
			*dFdW_ptr[InddFdW++] += -uv;
			*dFdW_ptr[InddFdW++] += -v2+alpha;
			*dFdW_ptr[InddFdW++] += -vw;

			// var 2
			*dFdW_ptr[InddFdW++] +=  v;
			*dFdW_ptr[InddFdW++] += -GM1*u;
			*dFdW_ptr[InddFdW++] +=  0.0;

			// var 3
			*dFdW_ptr[InddFdW++] +=  u;
			*dFdW_ptr[InddFdW++] += -GM3*v;
			*dFdW_ptr[InddFdW++] +=  w;

			// var 4
			*dFdW_ptr[InddFdW++] +=  0.0;
			*dFdW_ptr[InddFdW++] += -GM1*w;
			*dFdW_ptr[InddFdW++] +=  v;

			// var 5
			*dFdW_ptr[InddFdW++] +=  0.0;
			*dFdW_ptr[InddFdW++] +=  GM1;
			*dFdW_ptr[InddFdW++] +=  0.0;

			// *** eq 4 ***
			// var 1
			*dFdW_ptr[InddFdW++] += -uw;
			*dFdW_ptr[InddFdW++] += -vw;
			*dFdW_ptr[InddFdW++] += -w2+alpha;

			// var 2
			*dFdW_ptr[InddFdW++] +=  w;
			*dFdW_ptr[InddFdW++] +=  0.0;
			*dFdW_ptr[InddFdW++] += -GM1*u;

			// var 3
			*dFdW_ptr[InddFdW++] +=  0.0;
			*dFdW_ptr[InddFdW++] +=  w;
			*dFdW_ptr[InddFdW++] += -GM1*v;

			// var 4
			*dFdW_ptr[InddFdW++] +=  u;
			*dFdW_ptr[InddFdW++] +=  v;
			*dFdW_ptr[InddFdW++] += -GM3*w;

			// var 5
			*dFdW_ptr[InddFdW++] +=  0.0;
			*dFdW_ptr[InddFdW++] +=  0.0;
			*dFdW_ptr[InddFdW++] +=  GM1;

			// *** eq 5 ***
			// var 1
			*dFdW_ptr[InddFdW++] +=  u*beta;
			*dFdW_ptr[InddFdW++] +=  v*beta;
			*dFdW_ptr[InddFdW++] +=  w*beta;

			// var 2
			*dFdW_ptr[InddFdW++] +=  H-GM1*u2;
			*dFdW_ptr[InddFdW++] += -GM1*uv;
			*dFdW_ptr[InddFdW++] += -GM1*uw;

			// var 3
			*dFdW_ptr[InddFdW++] += -GM1*uv;
			*dFdW_ptr[InddFdW++] +=  H-GM1*v2;
			*dFdW_ptr[InddFdW++] += -GM1*vw;

			// var 4
			*dFdW_ptr[InddFdW++] += -GM1*uw;
			*dFdW_ptr[InddFdW++] += -GM1*vw;
			*dFdW_ptr[InddFdW++] +=  H-GM1*w2;

			// var 5
			*dFdW_ptr[InddFdW++] +=  GAMMA*u;
			*dFdW_ptr[InddFdW++] +=  GAMMA*v;
			*dFdW_ptr[InddFdW++] +=  GAMMA*w;

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

			// F = Null in the implicit solver
			if (F != NULL) {
				ptrdiff_t IndF = 0;
				// eq 1
				*F_ptr[IndF++] += rhou;
				*F_ptr[IndF++] += rhov;
				IndF += 1;

				// eq 2
				*F_ptr[IndF++] += rhou*u + p;
				*F_ptr[IndF++] += rhou*v;
				IndF += 1;

				// eq 3
				*F_ptr[IndF++] += rhov*u;
				*F_ptr[IndF++] += rhov*v + p;
				IndF += 1;

				// eq 4
				*F_ptr[IndF++] += (E+p)*u;
				*F_ptr[IndF++] += (E+p)*v;

				for (i = 0, iMax = Neq*DMAX; i < iMax; i++)
					F_ptr[i]++;
			}

			InddFdW = 0;
			// *** eq 1 ***
			// var 1
			*dFdW_ptr[InddFdW++] +=  0.0;
			*dFdW_ptr[InddFdW++] +=  0.0;
			InddFdW += 1;

			// var 2
			*dFdW_ptr[InddFdW++] +=  1.0;
			*dFdW_ptr[InddFdW++] +=  0.0;
			InddFdW += 1;

			// var 3
			*dFdW_ptr[InddFdW++] +=  0.0;
			*dFdW_ptr[InddFdW++] +=  1.0;
			InddFdW += 1;

			// var 4
			*dFdW_ptr[InddFdW++] +=  0.0;
			*dFdW_ptr[InddFdW++] +=  0.0;
			InddFdW += 1;

			// *** eq 2 ***
			// var 1
			*dFdW_ptr[InddFdW++] += -u2+alpha;
			*dFdW_ptr[InddFdW++] += -uv;
			InddFdW += 1;

			// var 2
			*dFdW_ptr[InddFdW++] += -GM3*u;
			*dFdW_ptr[InddFdW++] +=  v;
			InddFdW += 1;

			// var 3
			*dFdW_ptr[InddFdW++] += -GM1*v;
			*dFdW_ptr[InddFdW++] +=  u;
			InddFdW += 1;

			// var 4
			*dFdW_ptr[InddFdW++] +=  GM1;
			*dFdW_ptr[InddFdW++] +=  0.0;
			InddFdW += 1;

			// *** eq 3 ***
			// var 1
			*dFdW_ptr[InddFdW++] += -uv;
			*dFdW_ptr[InddFdW++] += -v2+alpha;
			InddFdW += 1;

			// var 2
			*dFdW_ptr[InddFdW++] +=  v;
			*dFdW_ptr[InddFdW++] += -GM1*u;
			InddFdW += 1;

			// var 3
			*dFdW_ptr[InddFdW++] +=  u;
			*dFdW_ptr[InddFdW++] += -GM3*v;
			InddFdW += 1;

			// var 4
			*dFdW_ptr[InddFdW++] +=  0.0;
			*dFdW_ptr[InddFdW++] +=  GM1;
			InddFdW += 1;

			// *** eq 4 ***
			// var 1
			*dFdW_ptr[InddFdW++] +=  u*beta;
			*dFdW_ptr[InddFdW++] +=  v*beta;
			InddFdW += 1;

			// var 2
			*dFdW_ptr[InddFdW++] +=  H-GM1*u2;
			*dFdW_ptr[InddFdW++] += -GM1*uv;
			InddFdW += 1;

			// var 3
			*dFdW_ptr[InddFdW++] += -GM1*uv;
			*dFdW_ptr[InddFdW++] +=  H-GM1*v2;
			InddFdW += 1;

			// var 4
			*dFdW_ptr[InddFdW++] +=  GAMMA*u;
			*dFdW_ptr[InddFdW++] +=  GAMMA*v;

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
				ptrdiff_t IndF = 0;
				// eq 1
				*F_ptr[IndF++] += rhou;
				IndF += 2;

				// eq 2
				*F_ptr[IndF++] += rhou*u + p;
				IndF += 2;

				// eq 3
				*F_ptr[IndF++] += (E+p)*u;
				IndF += 2;

				for (i = 0, iMax = Neq*DMAX; i < iMax; i++)
					F_ptr[i]++;
			}

			InddFdW = 0;
			// *** eq 1 ***
			// var 1
			*dFdW_ptr[InddFdW++] +=  0.0;
			InddFdW += 2;

			// var 2
			*dFdW_ptr[InddFdW++] +=  1.0;
			InddFdW += 2;

			// var 3
			*dFdW_ptr[InddFdW++] +=  0.0;
			InddFdW += 2;

			// *** eq 2 ***
			// var 1
			*dFdW_ptr[InddFdW++] += -u2+alpha;
			InddFdW += 2;

			// var 2
			*dFdW_ptr[InddFdW++] += -GM3*u;
			InddFdW += 2;

			// var 3
			*dFdW_ptr[InddFdW++] +=  GM1;
			InddFdW += 2;

			// *** eq 3 ***
			// var 1
			*dFdW_ptr[InddFdW++] +=  u*beta;
			InddFdW += 2;

			// var 2
			*dFdW_ptr[InddFdW++] +=  H-GM1*u2;
			InddFdW += 2;

			// var 3
			*dFdW_ptr[InddFdW++] +=  GAMMA*u;

			rho_ptr++; rhou_ptr++; E_ptr++;
			for (i = 0, iMax = Neq*Nvar*DMAX; i < iMax; i++)
				dFdW_ptr[i]++;
		}
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
