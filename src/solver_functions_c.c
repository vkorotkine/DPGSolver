// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "solver_functions_c.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_VOLUME.h"
#include "S_FACE.h"
#include "S_OpCSR.h"

#include "solver_functions.h"
#include "matrix_functions.h"
#include "sum_factorization.h"

#include "array_swap.h"
#include "boundary_conditions_c.h"
#include "fluxes_inviscid_c.h"
#include "fluxes_viscous_c.h"

#include "array_free.h"
#include "array_print.h"

/*
 *	Purpose:
 *		Provide solver related functions for complex variables (used for complex step verification).
 *
 *	Comments:
 *		Only supports imex_type = 'E' as the complex solver functions are not used to compute linearizations.
 *
 *	Notation:
 *
 *	References:
 */

// **************************************************************************************************** //
// VOLUME functions
// **************************************************************************************************** //

void coef_to_values_vI_c(struct S_VDATA const *const VDATA, char const coef_type)
{
	if (DB.Collocated) {
		EXIT_UNSUPPORTED; // Set A_vI = VOLUME->Ahat in calling function.
	} else {
		if (!(coef_type == 'W' || coef_type == 'Q'))
			EXIT_UNSUPPORTED;

		unsigned int const d    = DB.d,
		                   Nvar = DB.Nvar;

		struct S_OPERATORS_V const *const *const OPS    = (struct S_OPERATORS_V const *const *const) VDATA->OPS;
		struct S_VOLUME      const *const        VOLUME = VDATA->VOLUME;

		if (coef_type == 'W') {
			mm_dcc(CBCM,CBT,CBNT,OPS[0]->NvnI,Nvar,OPS[0]->NvnS,1.0,0.0,OPS[0]->ChiS_vI,VOLUME->What_c,VDATA->W_vI_c);
		} else if (coef_type == 'Q') {
			for (size_t dim = 0; dim < d; dim++)
				mm_dcc(CBCM,CBT,CBNT,OPS[0]->NvnI,Nvar,OPS[0]->NvnS,1.0,0.0,OPS[0]->ChiS_vI,VOLUME->Qhat_c[dim],VDATA->Q_vI_c[dim]);
		}
	}
}

void convert_between_rp_c(unsigned int const Nn, unsigned int const Nrc, double const *const C,
                          double complex *const Ap, double complex *const Ar, char const *const conv_type)
{
	unsigned int d = DB.d;

	if (strstr(conv_type,"FluxToRef")) {
		memset(Ar,0.0,Nn*Nrc*d * sizeof *Ar);
		for (size_t col = 0; col < Nrc; col++) {
		for (size_t dim1 = 0; dim1 < d; dim1++) {
			size_t const IndAr = (Nrc*dim1+col)*Nn;
			for (size_t dim2 = 0; dim2 < d; dim2++) {
				size_t const IndAp = (col*d+dim2)*Nn,
				             IndC  = (dim1*d+dim2)*Nn;
				for (size_t n = 0; n < Nn; n++)
					Ar[IndAr+n] += Ap[IndAp+n]*C[IndC+n];
			}
		}}
	} else {
		EXIT_UNSUPPORTED;
	}
}

void finalize_VOLUME_Inviscid_Weak_c(unsigned int const Nrc, double complex const *const Ar_vI,
                                     double complex *const RLHS, char const imex_type,
                                     struct S_VDATA const *const VDATA)
{
	struct S_OPERATORS_V const *const *const OPS = (struct S_OPERATORS_V const *const *const) VDATA->OPS;

	unsigned int const d    = DB.d,
	                   NvnS = OPS[0]->NvnS,
	                   NvnI = OPS[0]->NvnI;

	if (imex_type == 'E') {
		memset(RLHS,0.0,NvnS*Nrc * sizeof *RLHS);
		for (size_t dim = 0; dim < d; dim++)
			mm_dcc(CBCM,CBT,CBNT,NvnS,Nrc,NvnI,1.0,1.0,OPS[0]->D_Weak[dim],&Ar_vI[NvnI*Nrc*dim],RLHS);
	} else {
		EXIT_UNSUPPORTED;
	}
}

void finalize_VOLUME_Viscous_Weak_c(unsigned int const Nrc, double complex *const Ar_vI, double complex *const RLHS,
                                    char const imex_type, struct S_VDATA const *const VDATA)
{
	finalize_VOLUME_Inviscid_Weak_c(Nrc,Ar_vI,RLHS,imex_type,VDATA);
}

// **************************************************************************************************** //
// FACE functions
// **************************************************************************************************** //

void coef_to_values_fI_c(struct S_FDATA const *const FDATA, char const coef_type)
{
	if (!(coef_type == 'W' || coef_type == 'Q'))
		EXIT_UNSUPPORTED;

	struct S_OPERATORS_F const *const *const OPS    = (struct S_OPERATORS_F const *const *const) FDATA->OPS;
	struct S_VOLUME      const *const        VOLUME = FDATA->VOLUME;

	unsigned int const d        = DB.d,
	                   Nvar     = DB.Nvar,
	                   Vf       = FDATA->Vf,
	                   IndFType = FDATA->IndFType,
	                   NfnI     = OPS[IndFType]->NfnI;

	if (coef_type == 'W') {
		mm_dcc(CBCM,CBT,CBNT,NfnI,Nvar,OPS[0]->NvnS,1.0,0.0,OPS[0]->ChiS_fI[Vf],VOLUME->What_c,FDATA->W_fIL_c);
	} else if (coef_type == 'Q') {
		for (size_t dim = 0; dim < d; dim++)
			mm_dcc(CBCM,CBT,CBNT,NfnI,Nvar,OPS[0]->NvnS,1.0,0.0,OPS[0]->ChiS_fI[Vf],VOLUME->QhatV_c[dim],FDATA->GradW_fIL_c[dim]);
	}
}

void compute_WR_fIL_c(struct S_FDATA const *const FDATA, double complex const *const WL_fIL,
                      double complex *const WR_fIL)
{
	struct S_OPERATORS_F const *const *const OPS  = (struct S_OPERATORS_F const *const *const) FDATA->OPS;
	struct S_FACE        const *const        FACE = FDATA->FACE;

	unsigned int const d    = DB.d,
	                   Nvar = DB.Nvar,
	                   IndFType = FDATA->IndFType,
	                   BC      = FACE->BC,
	                   NfnI    = OPS[IndFType]->NfnI,
	                   *nOrdRL = OPS[IndFType]->nOrdRL;

	double const       *const XYZ_fIL = FACE->XYZ_fI,
	                   *const n_fIL   = FACE->n_fI;

	if (BC == 0 || (BC % BC_STEP_SC > 50)) { // Internal/Periodic FACE
		coef_to_values_fI_c(FDATA,'W');
		array_rearrange_cmplx(NfnI,Nvar,nOrdRL,'C',WR_fIL);
	} else { // Boundary FACE
		struct S_BC *const BCdata = malloc(sizeof *BCdata); // free

		BCdata->d   = DB.d;
		BCdata->Nn  = NfnI;
		BCdata->Nel = 1;

		BCdata->XYZ    = XYZ_fIL;
		BCdata->nL     = n_fIL;
		BCdata->WL_c     = WL_fIL;
		BCdata->WB_c     = WR_fIL;
		BCdata->GradWL_c = NULL;
		BCdata->GradWB_c = NULL;

		if (BC % BC_STEP_SC == BC_RIEMANN) {
			boundary_Riemann_c(NfnI,1,XYZ_fIL,WL_fIL,NULL,WR_fIL,n_fIL,d);
		} else if (BC % BC_STEP_SC == BC_SLIPWALL) {
			boundary_SlipWall_c(NfnI,1,WL_fIL,WR_fIL,n_fIL,d);
		} else if (BC % BC_STEP_SC == BC_BACKPRESSURE) {
			boundary_BackPressure_c(NfnI,1,WL_fIL,WR_fIL,n_fIL,d,Nvar);
		} else if (BC % BC_STEP_SC == BC_TOTAL_TP) {
			boundary_Total_TP_c(NfnI,1,XYZ_fIL,WL_fIL,WR_fIL,n_fIL,d,Nvar);
		} else if (BC % BC_STEP_SC == BC_SUPERSONIC_IN) {
			boundary_SupersonicInflow_c(NfnI,1,XYZ_fIL,WL_fIL,WR_fIL,n_fIL,d,Nvar);
		} else if (BC % BC_STEP_SC == BC_SUPERSONIC_OUT) {
			boundary_SupersonicOutflow_c(NfnI,1,XYZ_fIL,WL_fIL,WR_fIL,n_fIL,d,Nvar);
		} else if (BC % BC_STEP_SC == BC_NOSLIP_T) {
			boundary_NoSlip_Dirichlet_c(BCdata);
		} else if (BC % BC_STEP_SC == BC_NOSLIP_ADIABATIC) {
			boundary_NoSlip_Adiabatic_c(BCdata);
		} else {
			EXIT_UNSUPPORTED;
		}
		free(BCdata);
	}
}

void compute_WR_GradWR_fIL_c(struct S_FDATA const *const FDATA, double complex const *const WL_fIL,
                             double complex *const WR_fIL, double complex const *const *const GradWL_fIL,
                             double complex *const *const GradWR_fIL)
{
	struct S_OPERATORS_F const *const *const OPS  = FDATA->OPS;
	struct S_FACE        const *const        FACE = FDATA->FACE;

	unsigned int const d        = DB.d,
	                   Nvar     = DB.Nvar,
	                   IndFType = FDATA->IndFType,
	                   BC       = FACE->BC,
	                   NfnI     = OPS[IndFType]->NfnI,
	                   *const nOrdRL = OPS[IndFType]->nOrdRL;

	if (BC == 0 || (BC % BC_STEP_SC > 50)) { // Internal/Periodic FACE
		compute_WR_fIL_c(FDATA,WL_fIL,WR_fIL);

		coef_to_values_fI_c(FDATA,'Q');
		for (size_t dim = 0; dim < d; dim++)
			array_rearrange_cmplx(NfnI,Nvar,nOrdRL,'C',GradWR_fIL[dim]);
	} else { // Boundary FACE
		struct S_BC *const BCdata = malloc(sizeof *BCdata); // free

		BCdata->d   = DB.d;
		BCdata->Nn  = NfnI;
		BCdata->Nel = 1;

		BCdata->XYZ    = FACE->XYZ_fI;
		BCdata->nL     = FACE->n_fI;
		BCdata->WL_c     = WL_fIL;
		BCdata->WB_c     = WR_fIL;
		BCdata->GradWL_c = GradWL_fIL;
		BCdata->GradWB_c = GradWR_fIL;

		if (BC % BC_STEP_SC == BC_NOSLIP_T) {
			boundary_NoSlip_Dirichlet_c(BCdata);
		} else if (BC % BC_STEP_SC == BC_NOSLIP_ADIABATIC) {
			boundary_NoSlip_Adiabatic_c(BCdata);
		}
		free(BCdata);
	}
}

void compute_numerical_flux_c(struct S_FDATA const *const FDATA, char const imex_type)
{
	if (imex_type != 'E')
		EXIT_UNSUPPORTED;

	struct S_OPERATORS_F const *const *const OPS  = (struct S_OPERATORS_F const *const *const) FDATA->OPS;
	struct S_FACE        const *const        FACE = FDATA->FACE;

	unsigned int const d        = DB.d,
	                   Neq      = DB.Neq,
	                   IndFType = FDATA->IndFType,
	                   NfnI     = OPS[IndFType]->NfnI;

	double const *const n_fIL = FACE->n_fI;

	double complex const *const WL_fIL       = FDATA->NFluxData->WL_fIL_c,
	                     *const WR_fIL       = FDATA->NFluxData->WR_fIL_c;
	double complex       *const nFluxNum_fIL = FDATA->NFluxData->nFluxNum_fI_c;

	switch (DB.InviscidFluxType) {
	case FLUX_LF:
		flux_LF_c(NfnI,1,WL_fIL,WR_fIL,nFluxNum_fIL,n_fIL,d,Neq);
		break;
	case FLUX_ROE:
		flux_Roe_c(NfnI,1,WL_fIL,WR_fIL,nFluxNum_fIL,n_fIL,d,Neq);
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
}

void compute_numerical_solution_c(struct S_FDATA const *const FDATA, char const imex_type)
{
	if (!(imex_type == 'E'))
		EXIT_UNSUPPORTED;

	struct S_OPERATORS_F const *const *const OPS  = (struct S_OPERATORS_F const *const *const) FDATA->OPS;
	struct S_FACE        const *const        FACE = FDATA->FACE;

	unsigned int const d        = DB.d,
	                   Nvar     = d+2,
	                   IndFType = FDATA->IndFType,
	                   NfnI     = OPS[IndFType]->NfnI;

	double const *const n_fIL = FACE->n_fI;

	double complex const *const WL_fIL = FDATA->NFluxData->WL_fIL_c,
	                     *const WR_fIL = FDATA->NFluxData->WR_fIL_c;

	double complex       *const *const nSolNum_fIL = FDATA->NFluxData->nSolNum_fI_c;

	for (size_t dim = 0; dim < d; dim++) {
	for (size_t var = 0; var < Nvar; var++) {
	for (size_t n = 0; n < NfnI; n++) {
		nSolNum_fIL[dim][var*NfnI+n] = n_fIL[n*d+dim]*0.5*(WL_fIL[var*NfnI+n]+WR_fIL[var*NfnI+n]);
	}}}
}

static void correct_numerical_solution_strong_c(struct S_FDATA const *const FDATA, char const imex_type, char const side,
                                                char const Form_MF1)
{
	if (Form_MF1 != 'S')
		EXIT_UNSUPPORTED;

	if (!(imex_type == 'E'))
		EXIT_UNSUPPORTED;

	struct S_OPERATORS_F const *const *const OPS  = (struct S_OPERATORS_F const *const *const) FDATA->OPS;
	struct S_FACE        const *const        FACE = FDATA->FACE;

	unsigned int const d        = DB.d,
	                   Nvar     = d+2,
	                   IndFType = FDATA->IndFType,
	                   NfnI     = OPS[IndFType]->NfnI;

	double const *const n_fIL     = FACE->n_fI,
	             *const detJF_fIL = FACE->detJF_fI;

	double complex const *const WL_fIL = FDATA->NFluxData->WL_fIL_c,
	                     *const WR_fIL = FDATA->NFluxData->WR_fIL_c;

	double complex *const *const nSolNum_fIL = FDATA->NFluxData->nSolNum_fI_c;

	if (side == 'L') {
		for (size_t dim = 0; dim < d; dim++) {
		for (size_t var = 0; var < Nvar; var++) {
		for (size_t n = 0; n < NfnI; n++) {
			nSolNum_fIL[dim][var*NfnI+n] -= n_fIL[n*d+dim]*detJF_fIL[n]*WL_fIL[var*NfnI+n];
		}}}
	} else if (side == 'R') {
		for (size_t dim = 0; dim < d; dim++) {
		for (size_t var = 0; var < Nvar; var++) {
		for (size_t n = 0; n < NfnI; n++) {
			nSolNum_fIL[dim][var*NfnI+n] += n_fIL[n*d+dim]*detJF_fIL[n]*(WL_fIL[var*NfnI+n]-WR_fIL[var*NfnI+n]);
		}}}
	} else {
		EXIT_UNSUPPORTED;
	}
}

void compute_numerical_flux_viscous_c(struct S_FDATA const *const FDATAL, struct S_FDATA const *const FDATAR,
                                      char const imex_type)
{
	unsigned int const EvalType = 1;

	if (!(imex_type == 'E'))
		EXIT_UNSUPPORTED;

	struct S_OPERATORS_F const *const *const OPSL = (struct S_OPERATORS_F const *const *const) FDATAL->OPS,
	                           *const *const OPSR = (struct S_OPERATORS_F const *const *const) FDATAR->OPS;
	struct S_FACE        const *const        FACE = FDATAL->FACE;

	unsigned int const d        = DB.d,
	                   Nvar     = d+2,
	                   Neq      = d+2,
	                   IndFType = FDATAL->IndFType,
	                   Boundary = FACE->Boundary,
	                   NfnI     = OPSL[IndFType]->NfnI,
	                   NvnSL    = OPSL[0]->NvnS,
	                   NvnSR    = OPSR[0]->NvnS,
	                   *const nOrdRL = OPSL[FDATAL->IndFType]->nOrdRL;

	double const *const n_fIL                = FACE->n_fI;
	double complex const *const WL_fIL               = FDATAL->NFluxData->WL_fIL_c,
	                     *const WR_fIL               = FDATAL->NFluxData->WR_fIL_c,
	                     *const *const GradWL_fIL    = (double complex const *const *const) FDATAL->GradW_fIL_c,
	                     *const *const GradWR_fIL    = (double complex const *const *const) FDATAR->GradW_fIL_c;

	double complex *const nFluxViscNum_fIL = FDATAL->NFluxData->nFluxViscNum_fI_c;


	double complex **Q_fIL = malloc(d * sizeof *Q_fIL); // free
	for (size_t dim = 0; dim < d; dim++)
		Q_fIL[dim] = malloc(NfnI*Nvar * sizeof *Q_fIL[dim]); // free

	double complex *const FluxViscNum_fIL = malloc(NfnI*d*Neq * sizeof *FluxViscNum_fIL); // free
	if (DB.ViscousFluxType == FLUX_CDG2) {
		double const chi = 0.5*(DB.NfMax);

		if (Boundary) {
			for (size_t dim = 0; dim < d; dim++) {
				mm_dcc(CBCM,CBT,CBNT,NfnI,Nvar,NvnSL,2.0*chi,0.0,OPSL[0]->ChiS_fI[FDATAL->Vf],FACE->QhatL_c[dim],Q_fIL[dim]);
				for (size_t i = 0; i < NfnI*Nvar; i++)
					Q_fIL[dim][i] += GradWL_fIL[dim][i];
			}

			double complex *const W_fIL = malloc(NfnI*Nvar * sizeof *W_fIL); // free
			for (size_t i = 0; i < NfnI*Nvar; i++)
				W_fIL[i] = 0.5*(WL_fIL[i]+WR_fIL[i]);

			flux_viscous_c(NfnI,1,W_fIL,(double complex const *const *const) Q_fIL,FluxViscNum_fIL);
			free(W_fIL);
		} else {
			// Evaluate from which side scaling should be computed based on area switch (Brdar(2012), eq. (4.5))
			struct S_VOLUME const *const VL = FDATAL->VOLUME,
			                      *const VR = FDATAR->VOLUME;

			double const *const detJVL_vIL = VL->detJV_vI,
			             *const detJVR_vIR = VR->detJV_vI;

			// Potentially scale by volume of reference element (necessary for mixed meshes?) (ToBeDeleted)
			double VolL = 0.0;
			for (size_t n = 0; n < NvnSL; n++)
				VolL = max(VolL,detJVL_vIL[n]);

			double VolR = 0.0;
			for (size_t n = 0; n < NvnSR; n++)
				VolR = max(VolR,detJVR_vIR[n]);

			if (EvalType == 1) {
				// Left side
				for (size_t dim = 0; dim < d; dim++) {
					memset(Q_fIL[dim],0.0,NfnI*Nvar * sizeof *Q_fIL[dim]);
					if (VolL <= VolR)
						mm_dcc(CBCM,CBT,CBNT,NfnI,Nvar,NvnSL,2.0*chi,0.0, OPSL[0]->ChiS_fI[FDATAL->Vf],FACE->QhatL_c[dim],Q_fIL[dim]);
					for (size_t i = 0; i < NfnI*Nvar; i++)
						Q_fIL[dim][i] += GradWL_fIL[dim][i];
				}

				double complex *const FluxViscL_fIL = malloc(NfnI*d*Neq * sizeof *FluxViscL_fIL); // free
				flux_viscous_c(NfnI,1,WL_fIL,(double complex const *const *const) Q_fIL,FluxViscL_fIL);

				// Right side
				for (size_t dim = 0; dim < d; dim++) {
					memset(Q_fIL[dim],0.0,NfnI*Nvar * sizeof *Q_fIL[dim]);
					if (!(VolL <= VolR)) {
						mm_dcc(CBCM,CBT,CBNT,NfnI,Nvar,NvnSR,2.0*chi,0.0, OPSR[0]->ChiS_fI[FDATAR->Vf],FACE->QhatR_c[dim],Q_fIL[dim]);
						array_rearrange_cmplx(NfnI,Nvar,nOrdRL,'C',Q_fIL[dim]);
					}
					for (size_t i = 0; i < NfnI*Nvar; i++)
						Q_fIL[dim][i] += GradWR_fIL[dim][i];
				}

				double complex *const FluxViscR_fIL = malloc(NfnI*d*Neq * sizeof *FluxViscR_fIL); // free
				flux_viscous_c(NfnI,1,WR_fIL,(const complex double *const *const) Q_fIL,FluxViscR_fIL);

				// Compute numerical flux
				for (size_t i = 0; i < NfnI*d*Neq; i++)
					FluxViscNum_fIL[i] = 0.5*(FluxViscL_fIL[i]+FluxViscR_fIL[i]);

				free(FluxViscL_fIL);
				free(FluxViscR_fIL);
			} else if (EvalType == 2) {
				for (size_t dim = 0; dim < d; dim++) {
					if (VolL <= VolR) {
						mm_dcc(CBCM,CBT,CBNT,NfnI,Nvar,NvnSL,2.0*chi,0.0, OPSL[0]->ChiS_fI[FDATAL->Vf],FACE->QhatL_c[dim],Q_fIL[dim]);
					} else {
						mm_dcc(CBCM,CBT,CBNT,NfnI,Nvar,NvnSR,2.0*chi,0.0, OPSR[0]->ChiS_fI[FDATAR->Vf],FACE->QhatR_c[dim],Q_fIL[dim]);
						array_rearrange_cmplx(NfnI,Nvar,nOrdRL,'C',Q_fIL[dim]);
					}
					for (size_t i = 0; i < NfnI*Nvar; i++) {
						Q_fIL[dim][i] += GradWL_fIL[dim][i];
						Q_fIL[dim][i] += GradWR_fIL[dim][i];
					}
				}

				double complex *const W_fIL = malloc(NfnI*Nvar * sizeof *W_fIL); // free
				for (size_t i = 0; i < NfnI*Nvar; i++)
					W_fIL[i] = 0.5*(WL_fIL[i]+WR_fIL[i]);

				flux_viscous_c(NfnI,1,W_fIL,(double complex const *const *const) Q_fIL,FluxViscNum_fIL);
				free(W_fIL);
			} else {
				EXIT_UNSUPPORTED;
			}
		}
	} else if (DB.ViscousFluxType == FLUX_BR2) {
		double const chi = 1.0*(DB.NfMax);

		if (Boundary) {
			for (size_t dim = 0; dim < d; dim++) {
				mm_dcc(CBCM,CBT,CBNT,NfnI,Nvar,NvnSL,chi,0.0,OPSL[0]->ChiS_fI[FDATAL->Vf],FACE->QhatL_c[dim],Q_fIL[dim]);
				for (size_t i = 0; i < NfnI*Nvar; i++)
					Q_fIL[dim][i] += GradWL_fIL[dim][i];
			}

			double complex *const W_fIL = malloc(NfnI*Nvar * sizeof *W_fIL); // free
			for (size_t i = 0; i < NfnI*Nvar; i++)
				W_fIL[i] = 0.5*(WL_fIL[i]+WR_fIL[i]);

			flux_viscous_c(NfnI,1,W_fIL,(double complex const *const *const) Q_fIL,FluxViscNum_fIL);
			free(W_fIL);
		} else {
			if (EvalType == 1) {
				// Left side
				for (size_t dim = 0; dim < d; dim++) {
					mm_dcc(CBCM,CBT,CBNT,NfnI,Nvar,NvnSL,chi,0.0,OPSL[0]->ChiS_fI[FDATAL->Vf],FACE->QhatL_c[dim],Q_fIL[dim]);
					for (size_t i = 0; i < NfnI*Nvar; i++)
						Q_fIL[dim][i] += GradWL_fIL[dim][i];
				}

				double complex *const FluxViscL_fIL = malloc(NfnI*d*Neq * sizeof *FluxViscL_fIL); // free
				flux_viscous_c(NfnI,1,WL_fIL,(double complex const *const *const) Q_fIL,FluxViscL_fIL);

				// Right side
				for (size_t dim = 0; dim < d; dim++) {
					mm_dcc(CBCM,CBT,CBNT,NfnI,Nvar,NvnSR,chi,0.0,OPSR[0]->ChiS_fI[FDATAR->Vf],FACE->QhatR_c[dim],Q_fIL[dim]);

					array_rearrange_cmplx(NfnI,Nvar,nOrdRL,'C',Q_fIL[dim]);
					for (size_t i = 0; i < NfnI*Nvar; i++)
						Q_fIL[dim][i] += GradWR_fIL[dim][i];
				}

				double complex *const FluxViscR_fIL = malloc(NfnI*d*Neq * sizeof *FluxViscR_fIL); // free
				flux_viscous_c(NfnI,1,WR_fIL,(double complex const *const *const) Q_fIL,FluxViscR_fIL);

				// Compute numerical flux
				for (size_t i = 0; i < NfnI*d*Neq; i++)
					FluxViscNum_fIL[i] = 0.5*(FluxViscL_fIL[i]+FluxViscR_fIL[i]);

				free(FluxViscL_fIL);
				free(FluxViscR_fIL);
			} else if (EvalType == 2) {
				for (size_t dim = 0; dim < d; dim++) {
					mm_dcc(CBCM,CBT,CBNT,NfnI,Nvar,NvnSR,chi,0.0,OPSR[0]->ChiS_fI[FDATAR->Vf],FACE->QhatR_c[dim],Q_fIL[dim]);
					array_rearrange_cmplx(NfnI,Nvar,nOrdRL,'C',Q_fIL[dim]);

					mm_dcc(CBCM,CBT,CBNT,NfnI,Nvar,NvnSL,chi,1.0,OPSL[0]->ChiS_fI[FDATAL->Vf],FACE->QhatL_c[dim],Q_fIL[dim]);
					for (size_t i = 0; i < NfnI*Nvar; i++) {
						Q_fIL[dim][i] += GradWL_fIL[dim][i];
						Q_fIL[dim][i] += GradWR_fIL[dim][i];
					}
				}

				double complex *const W_fIL = malloc(NfnI*Nvar * sizeof *W_fIL); // free
				for (size_t i = 0; i < NfnI*Nvar; i++)
					W_fIL[i] = 0.5*(WL_fIL[i]+WR_fIL[i]);

				flux_viscous_c(NfnI,1,W_fIL,(double complex const *const *const) Q_fIL,FluxViscNum_fIL);
				free(W_fIL);
			} else {
				EXIT_UNSUPPORTED;
			}
		}
	} else {
		EXIT_UNSUPPORTED;
	}
	array_free2_cmplx(d,Q_fIL);


	// Take dot product with normal vector
	for (size_t i = 0; i < NfnI*Neq; i++)
		nFluxViscNum_fIL[i] = 0.0;

	for (size_t eq = 0; eq < Neq; eq++) {
		double complex const *const FluxViscNum_ptr = &FluxViscNum_fIL[NfnI*d*eq];
		for (size_t dim = 0; dim < d; dim++) {
		for (size_t n = 0; n < NfnI; n++) {
			nFluxViscNum_fIL[n+NfnI*eq] += n_fIL[n*d+dim]*FluxViscNum_ptr[NfnI*dim+n];
		}}
	}
	free(FluxViscNum_fIL);

	// Modify nFluxViscNum_fIL to account for boundary conditions (if necessary)
	unsigned int const BC = FACE->BC;

	if (BC % BC_STEP_SC == BC_NOSLIP_ADIABATIC) {
		// Set last component of nFluxViscNum to zero
		size_t const eq = Neq-1;
		for (size_t n = 0; n < NfnI; n++)
			nFluxViscNum_fIL[n+NfnI*eq] = 0.0;
	}
}

void add_Jacobian_scaling_FACE_c(struct S_FDATA const *const FDATA, char const imex_type, char const coef_type)
{
	if (imex_type != 'E')
		EXIT_UNSUPPORTED;

	struct S_OPERATORS_F const *const *const OPS  = (struct S_OPERATORS_F const *const *const) FDATA->OPS;
	struct S_FACE        const *const        FACE = FDATA->FACE;

	unsigned int const d        = DB.d,
	                   Neq      = DB.Neq,
	                   IndFType = FDATA->IndFType,
	                   NfnI = OPS[IndFType]->NfnI;

	double const *const detJF_fIL = FACE->detJF_fI;

	if (coef_type == 'W' || coef_type == 'V') {
		double complex *nFluxNum_fIL;
		if (coef_type == 'W') {
			nFluxNum_fIL = (double complex *const) FDATA->NFluxData->nFluxNum_fI_c;
		} else if (coef_type == 'V') {
			nFluxNum_fIL = (double complex *const) FDATA->NFluxData->nFluxViscNum_fI_c;
		}

		for (size_t eq = 0; eq < Neq; eq++) {
			size_t const IndnF = eq*NfnI;
			for (size_t n = 0; n < NfnI; n++)
				nFluxNum_fIL[IndnF+n] *= detJF_fIL[n];
		}
	} else if (coef_type == 'Q') {
		double complex *const *const nSolNum_fIL = FDATA->NFluxData->nSolNum_fI_c;

		for (size_t dim = 0; dim < d; dim++) {
			for (size_t eq = 0; eq < Neq; eq++) {
				size_t const IndnF = eq*NfnI;
				for (size_t n = 0; n < NfnI; n++)
					nSolNum_fIL[dim][IndnF+n] *= detJF_fIL[n];
			}
		}
	} else {
		EXIT_UNSUPPORTED;
	}

}

static void swap_FACE_orientation_c(struct S_FDATA const *const FDATA, char const imex_type, char const coef_type)
{
	if (!(imex_type == 'E'))
		EXIT_UNSUPPORTED;

	struct S_OPERATORS_F const *const *const OPS  = (struct S_OPERATORS_F const *const *const) FDATA->OPS;

	unsigned int const d              = DB.d,
	                   Neq            = d+2,
	                   IndFType       = FDATA->IndFType,
	                   NfnI           = OPS[IndFType]->NfnI,
	                   *const nOrdLR  = OPS[IndFType]->nOrdLR;


	if (coef_type == 'W' || coef_type == 'V') {
		double complex *nFluxNum_fI;
		if (coef_type == 'W') {
			nFluxNum_fI = (double complex *const) FDATA->NFluxData->nFluxNum_fI_c;
		} else if (coef_type == 'V') {
			nFluxNum_fI = (double complex *const) FDATA->NFluxData->nFluxViscNum_fI_c;
		}

		for (size_t i = 0, iMax = Neq*NfnI; i < iMax; i++)
			nFluxNum_fI[i] *= -1.0;

		array_rearrange_cmplx(NfnI,Neq,nOrdLR,'C',nFluxNum_fI);
	} else if (coef_type == 'Q') {
		double complex *const *const nSolNum_fI = FDATA->NFluxData->nSolNum_fI_c;

		for (size_t dim = 0; dim < d; dim++) {
			for (size_t i = 0, iMax = Neq*NfnI; i < iMax; i++)
				nSolNum_fI[dim][i] *= -1.0;

			array_rearrange_cmplx(NfnI,Neq,nOrdLR,'C',nSolNum_fI[dim]);
		}
	} else {
		EXIT_UNSUPPORTED;
	}
}

void finalize_FACE_Inviscid_Weak_c(struct S_FDATA const *const FDATAL, struct S_FDATA const *const FDATAR,
                                   double complex const *const nANumL_fI, char const side, char const imex_type,
                                   char const coef_type)
{
	if (!(imex_type == 'E'))
		EXIT_UNSUPPORTED;

	struct S_FDATA const *FDATA = NULL;
	double complex       *RHS   = NULL;

	if (side == 'L') {
		FDATA = FDATAL;
		RHS   = FDATA->FACE->RHSIn_c;
	} else if (side == 'R') {
		FDATA = FDATAR;
		RHS   = FDATA->FACE->RHSOut_c;

		swap_FACE_orientation_c(FDATAR,imex_type,coef_type);
	} else {
		EXIT_UNSUPPORTED;
	}

	struct S_OPERATORS_F const *const *const OPS  = (struct S_OPERATORS_F const *const *const) FDATA->OPS;

	unsigned int const Neq      = DB.Neq,
	                   Vf       = FDATA->Vf,
	                   IndFType = FDATA->IndFType,
	                   NfnI     = OPS[IndFType]->NfnI;

	mm_dcc(CBCM,CBT,CBNT,OPS[0]->NvnS,Neq,NfnI,1.0,0.0,OPS[0]->I_Weak_FF[Vf],nANumL_fI,RHS);
}

void finalize_QhatF_Weak_c(struct S_FDATA const *const FDATAL, struct S_FDATA const *const FDATAR, char const side,
                           char const imex_type, char const FORM_MF1)
{
	if (!(imex_type == 'E'))
		EXIT_UNSUPPORTED;

	struct S_FACE const *const FACE = FDATAL->FACE;

	if (FORM_MF1 == 'S')
		correct_numerical_solution_strong_c(FDATAL,'E',side,FORM_MF1);

	struct S_FDATA const *FDATA = NULL;

	double complex *const *QhatF;
	if (side == 'L') {
		FDATA = FDATAL;
		QhatF = (double complex *const *const) FACE->QhatL_c;
	} else if (side == 'R') {
		FDATA = FDATAR;
		QhatF = (double complex *const *const) FACE->QhatR_c;
		swap_FACE_orientation_c(FDATAR,imex_type,'Q');
	} else {
		EXIT_UNSUPPORTED;
	}

	struct S_OPERATORS_F const *const *const OPS  = (struct S_OPERATORS_F const *const *const) FDATA->OPS;

	unsigned int const d        = DB.d,
	                   Neq      = d+2,
	                   Nvar     = d+2,
	                   Vf       = FDATA->Vf,
	                   IndFType = FDATA->IndFType,
	                   NfnI     = OPS[IndFType]->NfnI,
	                   NvnS     = OPS[0]->NvnS;

	double complex const *const *const nSolNum_fI = (double complex const *const *const) FDATA->NFluxData->nSolNum_fI_c;

	double complex *const InSNum = malloc(NvnS*Nvar * sizeof *InSNum); // free
	for (size_t dim = 0; dim < d; dim++) {
		// Note that there is a minus sign included in the definition of I_Weak_FF.
		mm_dcc(CBCM,CBT,CBNT,OPS[0]->NvnS,Neq,NfnI,1.0,0.0,OPS[0]->I_Weak_FF[Vf],nSolNum_fI[dim],InSNum);
		for (size_t i = 0, iMax = NvnS*Nvar; i < iMax; i++)
			QhatF[dim][i] = -InSNum[i];
	}
	free(InSNum);
}

void finalize_FACE_Viscous_Weak_c(struct S_FDATA const *const FDATAL, struct S_FDATA const *const FDATAR,
                                  double complex const *const nANumL_fI, char const side, char const imex_type,
                                  char const coef_type)
{
	finalize_FACE_Inviscid_Weak_c(FDATAL,FDATAR,nANumL_fI,side,imex_type,coef_type);
}
