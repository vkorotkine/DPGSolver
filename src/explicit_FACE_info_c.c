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

void explicit_FACE_info_c(void)
{
	compute_Inviscid_FACE_RHS_EFE();
}

static void compute_Inviscid_FACE_RHS_EFE(void)
{
	unsigned int Nvar = DB.Nvar,
	             Neq  = DB.Neq;

	struct S_OPERATORS_F *OPSL[2], *OPSR[2];
	struct S_FACE     *FACE;

	struct S_FDATA *FDATAL = malloc(sizeof *FDATAL), // free
	               *FDATAR = malloc(sizeof *FDATAR); // free
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
		for (FACE = DB.FACE; FACE; FACE = FACE->next) {
			init_FDATA(FDATAL,FACE,'L');
			init_FDATA(FDATAR,FACE,'R');

			// Compute WL_fIL and WR_fIL (i.e. as seen from the (L)eft VOLUME)
			unsigned int IndFType = FDATAL->IndFType;
			unsigned int NfnI     = OPSL[IndFType]->NfnI;

			FDATAL->W_fIL_c = malloc(NfnI*Nvar * sizeof *(FDATAL->W_fIL_c)), // free
			FDATAR->W_fIL_c = malloc(NfnI*Nvar * sizeof *(FDATAR->W_fIL_c)); // free

			coef_to_values_fI_c(FDATAL,'W');
			compute_WR_fIL_c(FDATAR,FDATAL->W_fIL_c,FDATAR->W_fIL_c);


			// Compute numerical flux as seen from the left VOLUME
			double complex *nFluxNum_fI_c = malloc(NfnI*Neq * sizeof *nFluxNum_fI_c); // free

			NFluxData->WL_fIL_c      = FDATAL->W_fIL_c;
			NFluxData->WR_fIL_c      = FDATAR->W_fIL_c;
			NFluxData->nFluxNum_fI_c = nFluxNum_fI_c;

			compute_numerical_flux_c(FDATAL,'E');
			add_Jacobian_scaling_FACE_c(FDATAL,'E');

			free(FDATAL->W_fIL_c);
			free(FDATAR->W_fIL_c);


			// Compute FACE RHS terms
			unsigned int NvnSL = OPSL[0]->NvnS,
			             NvnSR = OPSR[0]->NvnS;

			if (FACE->RHSIn_c)
				free(FACE->RHSIn_c);
			FACE->RHSIn_c  = calloc(NvnSL*Neq , sizeof *(FACE->RHSIn_c)); // keep

			if (FACE->RHSOut_c)
				free(FACE->RHSOut_c);
			FACE->RHSOut_c = calloc(NvnSR*Neq , sizeof *(FACE->RHSOut_c)); // keep

			finalize_FACE_Inviscid_Weak_c(FDATAL,FDATAR,'L','E');
			if (!FACE->Boundary)
				finalize_FACE_Inviscid_Weak_c(FDATAL,FDATAR,'R','E');

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
