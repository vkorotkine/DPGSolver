// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "fluxes_inviscid.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"

#include "fluxes_structs.h"
#include "matrix_structs.h"
#include "variable_functions.h"
#include "solver_Advection_functions.h"

/*
 *	Purpose:
 *		Compute inviscid fluxes from input W in conservative form.
 *
 *	Comments:
 *		It is assumed that inputs: W, n and outputs: F, nFluxNum are vectorized (i.e. the memory ordering is by equation
 *		and not by element).
 *		Flux functions were intentionally written with if statements for each dimension in order to avoid redundant
 *		calculations if d < 3.
 *		Try using BLAS calls for dot products and check if there is a speed-up. (ToBeDeleted)
 *		The Roe flux function is very difficult to parallelize due to the large amount of required memory.
 *		There is some additional discussion in Obward(2015) regarding incorrect dissipation of the Roe scheme at low
 *		mach numbers. The proposed method to correct this requires evaluation of a shock indicator on all FACEs
 *		adjacent to the inner VOLUME of the FACE on which the numerical flux is being computed (which seems quite
 *		expensive to do). In the same paper, there is a reference to Thornber(2008) which motivates a local normal
 *		velocity modification to correct the overly dissipative numerical flux. However, this is done only in the
 *		context of a finite volume scheme with MUSCL reconstruction. Perhaps look into whether this problem is important
 *		for high-order methods and add a fix similar to Thornber(2008). (ToBeModified)
 *
 *		Further discussion of Roe scheme advancements:
 *		1) Several changes to the Roe scheme have been proposed based on incorrect asymptotic scaling of certain terms.
 *		See Li(2016) x2, Qu(2015) for latest developments. Also note that the two papers by Li from 2016 recommend
 *		different modifications when using the Roe scheme in standard test cases vs. LES cases. However, the large
 *		majority of the results are tested in 1st order FV frameworks using MUSCL reconstruction and the results may
 *		hence not be as relevant for the high-order setting of this code. => Compare traditional Roe flux (only with
 *		simple entropy fix), improved Roe flux with several fixes (Li(2016)) and DPG before making conclusions.
 *		(ToBeModified).
 *		2) It was also noted that a lot of discussion regarding the superiority of Flux Vector Splitting methods over
 *		the Roe scheme was found. See Kitamura(2016), Kitamura(2013) (fig. 3) and contained references for AUSM, SLAU,
 *		LDFSS2001, and Roe scheme comparisons. If flux vector splitting is found to be acceptable, it is likely cheaper
 *		to compute than the Roe flux. => Compare with Flux Difference Splitting schemes and DPG above and make
 *		conclusions then. (ToBeModified).
 *
 *		The simple entropy fix is taken from Qu(2015), eq. 35.
 *
 *		The functions with '_M' subscripts are used to wrap calls to the original functions where the new matrix struct
 *		format is used.
 *
 *	Notation:
 *
 *	References:
 *		flux_LF : Based off of Hesthaven's NodalDG function.(https://github.com/tcew/nodal-dg/.../CFD2D/EulerLF2D)
 *		flux_Roe: Toro(2009)-Riemann_Solvers_and_Numerical_Methods_for_Fluid_Dynamics (Ch. 11.3)
 *		        : Obwald(2015)-L2Roe:_A_Low_Dissipation_Version_of_Roe's_Approximate_Riemann_Solver_for_Low_Mach_Numbers
 *				: Thornber(2008)-An_improved_reconstruction_method_for_compressible_flows_with_low_Mach_number_features
 *
 *		Additional: (ToBeModified)
 *		Kitamura(2013)-Towards_Shock-Stable_and_Accurate_Hypersonic_Heating_Computations-_A_New_Pressure_Flux_for_AUSM-family_Schemes
 *		Kitamura(2015)-Reduced_Dissipation_AUSM-family_Fluxes-_HR-SLAU2_and_HR-AUSMp-up_for_High_Resolution_Unsteady_Flow_Simulations
 *		Qu(2015)-A_New_Roe-Type_Scheme_for_All_Speeds
 *		Li(2016)-All-Speed_Roe_Scheme_for_the_Large_Eddy_Simulation_of_Homogeneous_Decaying_Turbulence
 *		Li(2016)-Cures_for_the_Expansion_Shock_and_the_Shock_Instability_of_the_Roe_Scheme
 */

       void flux_Euler     (struct S_FLUX *const FLUXDATA);
static void flux_Advection (struct S_FLUX *const FLUXDATA);

void flux_inviscid(struct S_FLUX *const FLUXDATA)
{
	switch(FLUXDATA->PDE_index) {
		case PDE_ADVECTION:    flux_Advection(FLUXDATA); break;
		case PDE_EULER:        // fallthrough
		case PDE_NAVIERSTOKES: flux_Euler(FLUXDATA);     break;
		default:               EXIT_UNSUPPORTED;         break;
	}
}

static void flux_LF     (struct S_NUMERICALFLUX *const NUMFLUXDATA);
static void flux_Roe    (struct S_NUMERICALFLUX *const NUMFLUXDATA);
static void flux_upwind (struct S_NUMERICALFLUX *const NUMFLUXDATA);

void flux_num_inviscid(struct S_NUMERICALFLUX *const NUMFLUXDATA)
{
	switch(NUMFLUXDATA->NumFluxInviscid_index) {
		case FLUX_LF:     flux_LF(NUMFLUXDATA);     break;
		case FLUX_ROE:    flux_Roe(NUMFLUXDATA);    break;
		case FLUX_UPWIND: flux_upwind(NUMFLUXDATA); break;
		default:          EXIT_UNSUPPORTED;         break;
	}
}

void flux_inviscid_MA (struct S_FLUX_MA *const FLUXDATA_MA)
{
	struct S_FLUX FLUXDATA;

	FLUXDATA.d   = FLUXDATA_MA->d;
	FLUXDATA.Nn  = FLUXDATA_MA->W->extents[0];
	FLUXDATA.Nel = 1;

	FLUXDATA.W = FLUXDATA_MA->W->data;
	FLUXDATA.F = FLUXDATA_MA->F->data;

	FLUXDATA.XYZ = FLUXDATA_MA->XYZ->data;

	flux_inviscid(&FLUXDATA);
}

void flux_Euler(struct S_FLUX *const FLUXDATA)
{
	/*
	 *	Comments:
	 *		The storage ordering of the fluxes (Node, then dimension, then equation) is chosen such that memory stride
	 *		is minimized when converting from physical to reference space.
	 */

	unsigned int const d       = FLUXDATA->d,
	                   Neq     = d+2,
	                   Nn      = FLUXDATA->Nn,
	                   Nel     = FLUXDATA->Nel,
	                   NnTotal = Nn*Nel;

	double const *const W = FLUXDATA->W;
	double       *const F = FLUXDATA->F;

	double const *rho_ptr  = &W[NnTotal*0],
	             *rhou_ptr = &W[NnTotal*1],
	             *E_ptr    = &W[NnTotal*(d+1)];

	double *F_ptr[d*Neq];
	for (size_t eq = 0; eq < Neq; eq++)  {
	for (size_t dim = 0; dim < d; dim++) {
		F_ptr[eq*d+dim] = &F[(eq*d+dim)*NnTotal];
	}}

	if (d == 3) {
		double const *rhov_ptr = &W[NnTotal*2],
		             *rhow_ptr = &W[NnTotal*3];

		for (size_t n = 0; n < NnTotal; n++) {
			double const rho  = *rho_ptr++,
			             rhou = *rhou_ptr++,
			             rhov = *rhov_ptr++,
			             rhow = *rhow_ptr++,
			             E    = *E_ptr++,

			             u   = rhou/rho,
			             v   = rhov/rho,
			             w   = rhow/rho,

			             p = GM1*(E-0.5*rho*(u*u+v*v+w*w));

			size_t IndF = 0;
			// eq 1
			*F_ptr[IndF++]++ = rhou;
			*F_ptr[IndF++]++ = rhov;
			*F_ptr[IndF++]++ = rhow;

			// eq 2
			*F_ptr[IndF++]++ = rhou*u + p;
			*F_ptr[IndF++]++ = rhou*v;
			*F_ptr[IndF++]++ = rhou*w;

			// eq 3
			*F_ptr[IndF++]++ = rhov*u;
			*F_ptr[IndF++]++ = rhov*v + p;
			*F_ptr[IndF++]++ = rhov*w;

			// eq 4
			*F_ptr[IndF++]++ = rhow*u;
			*F_ptr[IndF++]++ = rhow*v;
			*F_ptr[IndF++]++ = rhow*w + p;

			// eq 5
			*F_ptr[IndF++]++ = (E+p)*u;
			*F_ptr[IndF++]++ = (E+p)*v;
			*F_ptr[IndF++]++ = (E+p)*w;
		}
	} else if (d == 2) {
		double const *rhov_ptr = &W[NnTotal*2];

		for (size_t n = 0; n < NnTotal; n++) {
			double const rho  = *rho_ptr++,
			             rhou = *rhou_ptr++,
			             rhov = *rhov_ptr++,
			             E    = *E_ptr++,

			             u   = rhou/rho,
			             v   = rhov/rho,

			             p = GM1*(E-0.5*rho*(u*u+v*v));

			size_t IndF = 0;
			// eq 1
			*F_ptr[IndF++]++ = rhou;
			*F_ptr[IndF++]++ = rhov;

			// eq 2
			*F_ptr[IndF++]++ = rhou*u + p;
			*F_ptr[IndF++]++ = rhou*v;

			// eq 3
			*F_ptr[IndF++]++ = rhov*u;
			*F_ptr[IndF++]++ = rhov*v + p;

			// eq 4
			*F_ptr[IndF++]++ = (E+p)*u;
			*F_ptr[IndF++]++ = (E+p)*v;
		}
	} else if (d == 1) {
		for (size_t n = 0; n < NnTotal; n++) {
			double const rho  = *rho_ptr++,
			             rhou = *rhou_ptr++,
			             E    = *E_ptr++,

			             u   = rhou/rho,

			             p = GM1*(E-0.5*rho*(u*u));

			size_t IndF = 0;
			// eq 1
			*F_ptr[IndF++]++ = rhou;

			// eq 2
			*F_ptr[IndF++]++ = rhou*u + p;

			// eq 3
			*F_ptr[IndF++]++ = (E+p)*u;
		}
	}
}

static void flux_Advection(struct S_FLUX *const FLUXDATA)
{
	unsigned int const d       = FLUXDATA->d,
	                   Nn      = FLUXDATA->Nn,
	                   Nel     = FLUXDATA->Nel,
	                   NnTotal = Nn*Nel;

	double const *const W = FLUXDATA->W;
	double       *const F = FLUXDATA->F;

	double *F_ptr[d];
	for (size_t dim = 0; dim < d; dim++)
		F_ptr[dim] = &F[dim*NnTotal];

	double const *const b = compute_b_Advection(NnTotal,FLUXDATA->XYZ); // free

	for (size_t n = 0; n < NnTotal; n++) {
		for (size_t dim = 0; dim < d; dim++) {
			*F_ptr[dim] = b[dim*NnTotal+n]*W[n];
			F_ptr[dim]++;
		}
	}
	free((double *) b);
}


static void flux_LF(struct S_NUMERICALFLUX *const NUMFLUXDATA)
{
	unsigned int const d       = NUMFLUXDATA->d,
	                   Neq     = d+2,
	                   Nn      = NUMFLUXDATA->Nn,
	                   Nel     = NUMFLUXDATA->Nel,
	                   NnTotal = Nn*Nel;

	double const *const nL = NUMFLUXDATA->nL;

	double const *const WL = NUMFLUXDATA->WL,
	             *const WR = NUMFLUXDATA->WR;

	double       *const nFluxNum = NUMFLUXDATA->nFluxNum;

	double const *rhoL_ptr  = &WL[NnTotal*0],
	             *rhouL_ptr = &WL[NnTotal*1],
	             *EL_ptr    = &WL[NnTotal*(d+1)],

	             *rhoR_ptr  = &WR[NnTotal*0],
	             *rhouR_ptr = &WR[NnTotal*1],
	             *ER_ptr    = &WR[NnTotal*(d+1)];

	double const *n_ptr = nL;

	struct S_FLUX *const FLUXDATA = malloc(sizeof *FLUXDATA); // free
	FLUXDATA->d   = d;
	FLUXDATA->Nn  = 1;
	FLUXDATA->Nel = 1;

	double *nF_ptr[Neq];
	for (size_t eq = 0; eq < Neq; eq++)
		nF_ptr[eq] = &nFluxNum[eq*NnTotal];

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

			for (size_t eq = 0; eq < Neq; eq++) {
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

			for (size_t eq = 0; eq < Neq; eq++) {
				*nF_ptr[eq]++ = 0.5*( n1*((*FL1_ptr)+(*FR1_ptr))+n2*((*FL2_ptr)+(*FR2_ptr)) + maxV*(WLn[eq]-WRn[eq]));

				FL1_ptr += d;
				FL2_ptr += d;
				FR1_ptr += d;
				FR2_ptr += d;
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

			for (size_t eq = 0; eq < Neq; eq++) {
				*nF_ptr[eq]++ = 0.5*(n1*((*FL1_ptr)+(*FR1_ptr)) + maxV*(WLn[eq]-WRn[eq]));

				FL1_ptr += d;
				FR1_ptr += d;
			}
		}
	}
	free(FLUXDATA);
}

static void flux_Roe(struct S_NUMERICALFLUX *const NUMFLUXDATA)
{
	/*
	 *	Comments:
	 *		This is the Roe-Pike version of the scheme which is different from the original Roe scheme in that the wave
	 *		numbers are linearized for faster computation.
	 */

	unsigned int const d       = NUMFLUXDATA->d,
	                   Neq     = d+2,
	                   Nn      = NUMFLUXDATA->Nn,
	                   Nel     = NUMFLUXDATA->Nel,
	                   NnTotal = Nn*Nel;

	double const *const nL = NUMFLUXDATA->nL;

	double const *const WL = NUMFLUXDATA->WL,
	             *const WR = NUMFLUXDATA->WR;

	double       *const nFluxNum = NUMFLUXDATA->nFluxNum;

	double const *rhoL_ptr  = &WL[NnTotal*0],
	             *rhouL_ptr = &WL[NnTotal*1],
	             *EL_ptr    = &WL[NnTotal*(d+1)],

	             *rhoR_ptr  = &WR[NnTotal*0],
	             *rhouR_ptr = &WR[NnTotal*1],
	             *ER_ptr    = &WR[NnTotal*(d+1)];

	double const *n_ptr = nL;

	double *nF_ptr[Neq];
	for (size_t eq = 0; eq < Neq; eq++)
		nF_ptr[eq] = &nFluxNum[eq*NnTotal];

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
			double const l1   = min(fabs(Vn-c),fabs(VnL-c)),
			             l234 = fabs(Vn),
			             l5   = max(fabs(Vn+c),fabs(VnR+c));

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
			size_t IndnF = 0;
			*nF_ptr[IndnF++]++ = 0.5*(nF1 - dis1);
			*nF_ptr[IndnF++]++ = 0.5*(nF2 - dis2);
			*nF_ptr[IndnF++]++ = 0.5*(nF3 - dis3);
			*nF_ptr[IndnF++]++ = 0.5*(nF4 - dis4);
			*nF_ptr[IndnF++]++ = 0.5*(nF5 - dis5);
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
			double const l1   = min(fabs(Vn-c),fabs(VnL-c)),
			             l234 = fabs(Vn),
			             l5   = max(fabs(Vn+c),fabs(VnR+c));

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
			size_t IndnF = 0;
			*nF_ptr[IndnF++]++ = 0.5*(nF1 - dis1);
			*nF_ptr[IndnF++]++ = 0.5*(nF2 - dis2);
			*nF_ptr[IndnF++]++ = 0.5*(nF3 - dis3);
			*nF_ptr[IndnF++]++ = 0.5*(nF5 - dis5);
		}
	} else if (d == 1) {
		for (size_t n = 0; n < NnTotal; n++) {
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
			double const l1   = min(fabs(Vn-c),fabs(VnL-c)),
			             l234 = fabs(Vn),
			             l5   = max(fabs(Vn+c),fabs(VnR+c));

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
			size_t IndnF = 0;
			*nF_ptr[IndnF++]++ = 0.5*(nF1 - dis1);
			*nF_ptr[IndnF++]++ = 0.5*(nF2 - dis2);
			*nF_ptr[IndnF++]++ = 0.5*(nF5 - dis5);
		}
	}
}

static void flux_upwind(struct S_NUMERICALFLUX *const NUMFLUXDATA)
{
	unsigned int const d       = NUMFLUXDATA->d,
	                   Nn      = NUMFLUXDATA->Nn,
	                   Nel     = NUMFLUXDATA->Nel,
	                   NnTotal = Nn*Nel;

	double const *const nL = NUMFLUXDATA->nL;

	double const *const WL = NUMFLUXDATA->WL,
	             *const WR = NUMFLUXDATA->WR;

	double       *const nFluxNum = NUMFLUXDATA->nFluxNum;

	double const *const b = compute_b_Advection(NnTotal,NUMFLUXDATA->XYZ); // free

	for (size_t n = 0; n < NnTotal; n++) {
		double b_dot_n = 0.0;
		for (size_t dim = 0; dim < d; dim++)
			b_dot_n += b[dim*NnTotal+n]*nL[n*d+dim];

		if (b_dot_n >= 0.0)
			nFluxNum[n] = b_dot_n*WL[n];
		else
			nFluxNum[n] = b_dot_n*WR[n];
	}
	free((double *) b);
}
