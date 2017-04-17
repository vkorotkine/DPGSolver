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
		} else {
			EXIT_UNSUPPORTED;
		}
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

void add_Jacobian_scaling_FACE_c(struct S_FDATA const *const FDATA, char const imex_type)
{
	if (imex_type != 'E')
		EXIT_UNSUPPORTED;

	struct S_OPERATORS_F const *const *const OPS  = (struct S_OPERATORS_F const *const *const) FDATA->OPS;
	struct S_FACE        const *const        FACE = FDATA->FACE;

	unsigned int const Neq      = DB.Neq,
	                   IndFType = FDATA->IndFType,
	                   NfnI = OPS[IndFType]->NfnI;

	double const *const detJF_fIL = FACE->detJF_fI;

	double complex *const nFluxNum_fIL = FDATA->NFluxData->nFluxNum_fI_c;

	for (size_t eq = 0; eq < Neq; eq++) {
		size_t const IndnF = eq*NfnI;
		for (size_t n = 0; n < NfnI; n++)
			nFluxNum_fIL[IndnF+n] *= detJF_fIL[n];
	}
}

static void swap_FACE_orientation_c(struct S_FDATA const *const FDATA, char const imex_type)
{
	if (!(imex_type == 'E'))
		EXIT_UNSUPPORTED;

	struct S_OPERATORS_F const *const *const OPS  = (struct S_OPERATORS_F const *const *const) FDATA->OPS;

	unsigned int const Neq            = DB.Neq,
	                   IndFType       = FDATA->IndFType,
	                   NfnI           = OPS[IndFType]->NfnI,
	                   *const nOrdLR  = OPS[IndFType]->nOrdLR;

	double complex *const nFluxNum_fI = FDATA->NFluxData->nFluxNum_fI_c;

	for (size_t i = 0, iMax = Neq*NfnI; i < iMax; i++)
		nFluxNum_fI[i] *= -1.0;

	array_rearrange_cmplx(NfnI,Neq,nOrdLR,'C',nFluxNum_fI);
}

void finalize_FACE_Inviscid_Weak_c(struct S_FDATA const *const FDATAL, struct S_FDATA const *const FDATAR,
                                   char const side, char const imex_type)
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

		swap_FACE_orientation_c(FDATAR,'E');
	} else {
		EXIT_UNSUPPORTED;
	}

	struct S_OPERATORS_F const *const *const OPS  = (struct S_OPERATORS_F const *const *const) FDATA->OPS;

	unsigned int const Neq      = DB.Neq,
	                   Vf       = FDATA->Vf,
					   IndFType = FDATA->IndFType,
	                   NfnI     = OPS[IndFType]->NfnI;

	double complex const *nFluxNum_fI = FDATA->NFluxData->nFluxNum_fI_c;

	mm_dcc(CBCM,CBT,CBNT,OPS[0]->NvnS,Neq,NfnI,1.0,0.0,OPS[0]->I_Weak_FF[Vf],nFluxNum_fI,RHS);
}
