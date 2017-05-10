// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "explicit_GradW_c.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#include "Macros.h"
#include "Parameters.h"
#include "S_DB.h"
#include "S_VOLUME.h"
#include "S_FACE.h"

#include "solver_functions.h"
#include "solver_functions_c.h"
#include "matrix_functions.h"
#include "array_free.h"

/*
 *	Purpose:
 *		Identical to explicit_GradW using complex variables (for complex step verification).
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

static void explicit_GradW_VOLUME_c   (void);
static void explicit_GradW_FACE_c     (void);
static void explicit_GradW_finalize_c (void);

void explicit_GradW_c(void)
{
	if (!DB.Viscous)
		return;

	explicit_GradW_VOLUME_c();
	explicit_GradW_FACE_c();
	explicit_GradW_finalize_c();
}

static void explicit_GradW_VOLUME_c(void)
{
	unsigned int const d    = DB.d,
	                   Nvar = DB.Nvar;

	struct S_OPERATORS_V *OPS[2];

	struct S_VDATA *const VDATA = malloc(sizeof *VDATA); // free
	VDATA->OPS = (struct S_OPERATORS_V const *const *) OPS;

	for (size_t i = 0; i < 2; i++)
		OPS[i] = malloc(sizeof *OPS[i]); // free

	struct S_Dxyz *const DxyzInfo = malloc(sizeof *DxyzInfo); // free

	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		init_VDATA(VDATA,VOLUME);

		unsigned int const NvnS = VDATA->OPS[0]->NvnS,
		                   NvnI = VDATA->OPS[0]->NvnI;

		DxyzInfo->Nbf = VDATA->OPS[0]->NvnS;
		DxyzInfo->Nn  = VDATA->OPS[0]->NvnI;
		DxyzInfo->D   = (double const *const *const) VDATA->OPS[0]->D_Weak;
		DxyzInfo->C   = VOLUME->C_vI;

		double const *const        ChiS_vI  = VDATA->OPS[0]->ChiS_vI;
		for (size_t dim = 0; dim < d; dim++) {
			DxyzInfo->dim = dim;
			double *const Dxyz = compute_Dxyz(DxyzInfo,d); // free

			// Note: The detJ_vI term cancels with the gradient operator (Zwanenburg(2016), eq. (B.2))
			double *DxyzChiS = NULL;
			if (DB.Collocated) { // ChiS_vI == I
				DxyzChiS = Dxyz;
			} else {
				DxyzChiS = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS,NvnS,NvnI,1.0,Dxyz,ChiS_vI); // free
				free(Dxyz);
			}

			// Compute intermediate Qhat contribution
			if (VOLUME->QhatV_c[dim] != NULL)
				free(VOLUME->QhatV_c[dim]);
			VOLUME->QhatV_c[dim] = malloc(NvnS*Nvar * sizeof *(VOLUME->QhatV_c[dim])); // keep

			// Note: Using CBCM with CBNT for DxyzChiS (stored in row-major ordering) gives DxyzChiS transposed in the
			//       operation below.
			mm_dcc(CBCM,CBNT,CBNT,NvnS,Nvar,NvnS,1.0,0.0,DxyzChiS,VOLUME->What_c,VOLUME->QhatV_c[dim]);
			free(DxyzChiS);

			if (VOLUME->Qhat_c[dim] != NULL)
				free(VOLUME->Qhat_c[dim]);
			VOLUME->Qhat_c[dim] = malloc(NvnS*Nvar * sizeof *(VOLUME->Qhat_c[dim])); // keep

			for (size_t i = 0; i < NvnS*Nvar; i++)
				VOLUME->Qhat_c[dim][i] = VOLUME->QhatV_c[dim][i];
		}
	}

	free(VDATA);
	for (size_t i = 0; i < 2; i++)
		free(OPS[i]);
	free(DxyzInfo);
}

static void explicit_GradW_FACE_c(void)
{
	// Initialize DB Parameters
	unsigned int const d    = DB.d,
	                   Nvar = DB.Nvar,
	                   Neq  = DB.Neq;

	// Standard datatypes
	struct S_OPERATORS_F *OPSL[2], *OPSR[2];
	struct S_FDATA       *const FDATAL = malloc(sizeof *FDATAL), // free
	                     *const FDATAR = malloc(sizeof *FDATAR); // free
	FDATAL->OPS = (struct S_OPERATORS_F const *const *) OPSL;
	FDATAR->OPS = (struct S_OPERATORS_F const *const *) OPSR;

	struct S_NumericalFlux *const NFluxData = malloc(sizeof *NFluxData); // free
	FDATAL->NFluxData = NFluxData;
	FDATAR->NFluxData = NFluxData;

	for (size_t i = 0; i < 2; i++) {
		OPSL[i]  = malloc(sizeof *OPSL[i]);  // free
		OPSR[i]  = malloc(sizeof *OPSR[i]);  // free
	}

	for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
		init_FDATA(FDATAL,FACE,'L');
		init_FDATA(FDATAR,FACE,'R');

		// Compute WL_fIL and WR_fIL (i.e. as seen from the (L)eft VOLUME)
		unsigned int const IndFType = FDATAL->IndFType,
		                   NfnI     = OPSL[IndFType]->NfnI;

		FDATAL->W_fIL_c = malloc(NfnI*Nvar * sizeof *(FDATAL->W_fIL_c)), // free
		FDATAR->W_fIL_c = malloc(NfnI*Nvar * sizeof *(FDATAR->W_fIL_c)); // free

		coef_to_values_fI_c(FDATAL,'W','E');
		compute_WR_fIL_c(FDATAR,FDATAL->W_fIL_c,FDATAR->W_fIL_c);

		// Compute numerical flux as seen from the left VOLUME
		NFluxData->WL_fIL_c     = FDATAL->W_fIL_c;
		NFluxData->WR_fIL_c     = FDATAR->W_fIL_c;
		NFluxData->nSolNum_fI_c = malloc(d * sizeof *(NFluxData->nSolNum_fI_c)); // free
		for (size_t dim = 0; dim < d; dim++)
			NFluxData->nSolNum_fI_c[dim] = malloc(NfnI*Neq * sizeof *(NFluxData->nSolNum_fI_c[dim])); // free

		compute_numerical_solution_c(FDATAL,'E');
		add_Jacobian_scaling_FACE_c(FDATAL,'E','Q');

		for (size_t dim = 0; dim < d; dim++) {
			if (FACE->QhatL_c[dim] != NULL)
				free(FACE->QhatL_c[dim]);
			FACE->QhatL_c[dim] = malloc(OPSL[0]->NvnS*Nvar * sizeof *(FACE->QhatL_c[dim])); // keep

			if (!FACE->Boundary) {
				if (FACE->QhatR_c[dim] != NULL)
					free(FACE->QhatR_c[dim]);
				FACE->QhatR_c[dim] = malloc(OPSR[0]->NvnS*Nvar * sizeof *(FACE->QhatR_c[dim])); // keep
			}
		}

		finalize_QhatF_Weak_c(FDATAL,FDATAR,'L','E');
		if (!FACE->Boundary)
			finalize_QhatF_Weak_c(FDATAL,FDATAR,'R','E');

		free(FDATAL->W_fIL_c);
		free(FDATAR->W_fIL_c);
		array_free2_cmplx(d,NFluxData->nSolNum_fI_c);
	}

	for (size_t i = 0; i < 2; i++) {
		free(OPSL[i]);
		free(OPSR[i]);
	}
	free(NFluxData);
	free(FDATAL);
	free(FDATAR);
}

static void finalize_Qhat_c(struct S_VOLUME const *const VOLUME, unsigned int const NvnS,
                          double complex *const *const Qhat)
{
	unsigned int const d    = DB.d,
	                   Nvar = DB.Nvar;

	if (DB.Collocated) {
		double const *const detJV_vI = VOLUME->detJV_vI;
		for (size_t dim = 0; dim < d; dim++) {
			for (size_t var = 0; var < Nvar; var++) {
			for (size_t n = 0; n < NvnS; n++) {
				Qhat[dim][var*NvnS+n] /= detJV_vI[n];
			}}
		}
	} else {
		double complex *const Qhat_tmp = malloc(NvnS*Nvar * sizeof *Qhat_tmp); // free
		for (size_t dim = 0; dim < d; dim++) {
			mm_dcc(CBCM,CBT,CBNT,NvnS,Nvar,NvnS,1.0,0.0,VOLUME->MInv,Qhat[dim],Qhat_tmp);
			for (size_t var = 0; var < Nvar; var++) {
			for (size_t n = 0; n < NvnS; n++) {
				Qhat[dim][var*NvnS+n] = Qhat_tmp[var*NvnS+n];
			}}
		}
		free(Qhat_tmp);
	}
}

static void explicit_GradW_finalize_c(void)
{
	unsigned int const d    = DB.d,
	                   Nvar = DB.Nvar;

	// Add FACE contributions to VOLUME->Qhat then multiply by MInv
	for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
		struct S_VOLUME const *const VL = FACE->VL,
		                      *const VR = FACE->VR;

		unsigned int const NvnSL = VL->NvnS,
		                   NvnSR = VR->NvnS;

		for (size_t dim = 0; dim < d; dim++) {
			for (size_t i = 0; i < NvnSL*Nvar; i++)
				VL->Qhat_c[dim][i] += FACE->QhatL_c[dim][i];

			if (!FACE->Boundary) {
				for (size_t i = 0; i < NvnSR*Nvar; i++)
					VR->Qhat_c[dim][i] += FACE->QhatR_c[dim][i];
			}
		}

		finalize_Qhat_c(VL,NvnSL,FACE->QhatL_c);
		if (!FACE->Boundary)
			finalize_Qhat_c(VR,NvnSR,FACE->QhatR_c);
	}

	// Multiply VOLUME Qhat terms by MInv
	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		unsigned int const NvnS = VOLUME->NvnS;

		finalize_Qhat_c(VOLUME,NvnS,VOLUME->Qhat_c);
		finalize_Qhat_c(VOLUME,NvnS,VOLUME->QhatV_c);
	}
}
