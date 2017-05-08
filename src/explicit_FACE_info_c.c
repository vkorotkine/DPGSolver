// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "explicit_FACE_info_c.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_FACE.h"

#include "solver_functions.h"
#include "solver_functions_c.h"
#include "array_free.h"

/*
 *	Purpose:
 *		Identical to explicit_FACE_info using complex variables (for complex step verification).
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

static void compute_Inviscid_FACE_RHS_EFE (void);
static void compute_Viscous_FACE_RHS_EFE  (void);

void explicit_FACE_info_c(void)
{
	compute_Inviscid_FACE_RHS_EFE();
	compute_Viscous_FACE_RHS_EFE();
}

static void compute_Inviscid_FACE_RHS_EFE(void)
{
	unsigned int const d    = DB.d,
	                   Neq  = d+2,
	                   Nvar = d+2;

	struct S_OPERATORS_F *OPSL[2], *OPSR[2];

	struct S_FDATA *const FDATAL = malloc(sizeof *FDATAL), // free
	               *const FDATAR = malloc(sizeof *FDATAR); // free
	FDATAL->OPS = (struct S_OPERATORS_F const *const *) OPSL;
	FDATAR->OPS = (struct S_OPERATORS_F const *const *) OPSR;

	struct S_NumericalFlux *const NFluxData = malloc(sizeof *NFluxData); // free
	FDATAL->NFluxData = NFluxData;
	FDATAR->NFluxData = NFluxData;

	for (size_t i = 0; i < 2; i++) {
		OPSL[i] = malloc(sizeof *OPSL[i]); // free
		OPSR[i] = malloc(sizeof *OPSR[i]); // free
	}

	if (strstr(DB.Form,"Weak")) {
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
			NFluxData->WL_fIL_c      = FDATAL->W_fIL_c;
			NFluxData->WR_fIL_c      = FDATAR->W_fIL_c;
			NFluxData->nFluxNum_fI_c = malloc(NfnI*Neq * sizeof *(NFluxData->nFluxNum_fI_c)); // free

			compute_numerical_flux_c(FDATAL,'E');
			add_Jacobian_scaling_FACE_c(FDATAL,'E','W');

			free(FDATAL->W_fIL_c);
			free(FDATAR->W_fIL_c);


			// Compute FACE RHS terms
			unsigned int const NvnSL = OPSL[0]->NvnS,
			                   NvnSR = OPSR[0]->NvnS;

			if (FACE->RHSIn_c != NULL)
				free(FACE->RHSIn_c);
			FACE->RHSIn_c = calloc(NvnSL*Neq , sizeof *(FACE->RHSIn_c)); // keep
			finalize_FACE_Inviscid_Weak_c(FDATAL,FDATAR,NFluxData->nFluxNum_fI_c,'L','E','W');

			if (!FACE->Boundary) {
				if (FACE->RHSOut_c != NULL)
					free(FACE->RHSOut_c);
				FACE->RHSOut_c = calloc(NvnSR*Neq , sizeof *(FACE->RHSOut_c)); // keep

				finalize_FACE_Inviscid_Weak_c(FDATAL,FDATAR,NFluxData->nFluxNum_fI_c,'R','E','W');
			}

			free(NFluxData->nFluxNum_fI_c);
		}
	} else if (strstr(DB.Form,"Strong")) {
		EXIT_UNSUPPORTED;
	}

	for (size_t i = 0; i < 2; i++) {
		free(OPSL[i]);
		free(OPSR[i]);
	}
	free(NFluxData);
	free(FDATAL);
	free(FDATAR);
}

static void compute_Viscous_FACE_RHS_EFE(void)
{
	if (!DB.Viscous)
		return;

	unsigned int const d    = DB.d,
	                   Nvar = DB.Nvar,
	                   Neq  = DB.Neq;

	struct S_OPERATORS_F *OPSL[2], *OPSR[2];

	struct S_FDATA *const FDATAL = malloc(sizeof *FDATAL), // free
	               *const FDATAR = malloc(sizeof *FDATAR); // free
	FDATAL->OPS = (struct S_OPERATORS_F const *const *) OPSL;
	FDATAR->OPS = (struct S_OPERATORS_F const *const *) OPSR;

	struct S_NumericalFlux *const NFluxData = malloc(sizeof *NFluxData); // free
	FDATAL->NFluxData = NFluxData;
	FDATAR->NFluxData = NFluxData;

	for (size_t i = 0; i < 2; i++) {
		OPSL[i] = malloc(sizeof *OPSL[i]); // free
		OPSR[i] = malloc(sizeof *OPSR[i]); // free
	}

	if (strstr(DB.Form,"Weak")) {
		for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
			init_FDATA(FDATAL,FACE,'L');
			init_FDATA(FDATAR,FACE,'R');

			// Compute WL_fIL and WR_fIL and their gradients (as seen from the (L)eft VOLUME)
			unsigned int const IndFType = FDATAL->IndFType,
			                   NfnI     = OPSL[IndFType]->NfnI;

			FDATAL->W_fIL_c = malloc(NfnI*Nvar * sizeof *(FDATAL->W_fIL_c)), // free
			FDATAR->W_fIL_c = malloc(NfnI*Nvar * sizeof *(FDATAR->W_fIL_c)); // free

			double complex **const QpL_fIL_c = malloc(d * sizeof *QpL_fIL_c), // free
			               **const QpR_fIL_c = malloc(d * sizeof *QpR_fIL_c); // free

			for (size_t dim = 0; dim < d; dim++) {
				QpL_fIL_c[dim] = malloc(NfnI*Nvar * sizeof *QpL_fIL_c[dim]); // free
				QpR_fIL_c[dim] = malloc(NfnI*Nvar * sizeof *QpR_fIL_c[dim]); // free
			}

			FDATAL->Qp_fIL_c = QpL_fIL_c;
			FDATAR->Qp_fIL_c = QpR_fIL_c;

			coef_to_values_fI_c(FDATAL,'W','E');
			coef_to_values_fI_c(FDATAL,'Q','E');
			compute_WR_QpR_fIL_c(FDATAR,FDATAL->W_fIL_c,FDATAR->W_fIL_c,
			                     (double complex const *const *const) QpL_fIL_c,QpR_fIL_c,'E');


			// Compute numerical flux as seen from the left VOLUME
			NFluxData->WL_fIL_c          = FDATAL->W_fIL_c;
			NFluxData->WR_fIL_c          = FDATAR->W_fIL_c;
			NFluxData->nFluxViscNum_fI_c = malloc(NfnI*Neq * sizeof *(NFluxData->nFluxViscNum_fI_c)); // free

			compute_numerical_flux_viscous_c(FDATAL,FDATAR,'E');
			add_Jacobian_scaling_FACE_c(FDATAL,'E','V');

			free(FDATAL->W_fIL_c);
			free(FDATAR->W_fIL_c);
			array_free2_cmplx(d,QpL_fIL_c);
			array_free2_cmplx(d,QpR_fIL_c);

			// Compute FACE RHS terms
			finalize_FACE_Viscous_Weak_c(FDATAL,FDATAR,NFluxData->nFluxViscNum_fI_c,'L','E','V');
			if (!FACE->Boundary) {
				finalize_FACE_Viscous_Weak_c(FDATAL,FDATAR,NFluxData->nFluxViscNum_fI_c,'R','E','V');
			}

			free(NFluxData->nFluxViscNum_fI_c);
		}
	} else if (strstr(DB.Form,"Strong")) {
		// Note that the viscous flux is negated.
		EXIT_UNSUPPORTED;
	}

	for (size_t i = 0; i < 2; i++) {
		free(OPSL[i]);
		free(OPSR[i]);
	}
	free(NFluxData);
	free(FDATAL);
	free(FDATAR);
}
