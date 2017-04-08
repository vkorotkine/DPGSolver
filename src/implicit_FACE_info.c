// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "implicit_FACE_info.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_FACE.h"

#include "solver_functions.h"


/*
 *	Purpose:
 *		Evaluate the FACE contributions to the RHS and LHS terms.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

static void compute_FACE_EFE (void);

void implicit_FACE_info(void)
{
	compute_FACE_EFE();
}

static void compute_FACE_EFE(void)
{
	const char         *Form = DB.Form;
	const unsigned int Nvar  = DB.Nvar,
	                   Neq   = DB.Neq;

	struct S_OPERATORS_F *OPSL[2], *OPSR[2];
	struct S_FACE     *FACE;

	struct S_FDATA *FDATAL = malloc(sizeof *FDATAL), // free
	               *FDATAR = malloc(sizeof *FDATAR); // free
	FDATAL->OPS = OPSL;
	FDATAR->OPS = OPSR;

	struct S_NumericalFlux *NFluxData = malloc(sizeof *NFluxData); // free
	FDATAL->NFluxData = NFluxData;
	FDATAR->NFluxData = NFluxData;

	for (size_t i = 0; i < 2; i++) {
		OPSL[i] = malloc(sizeof *OPSL[i]); // free
		OPSR[i] = malloc(sizeof *OPSR[i]); // free
	}

	if (strstr(Form,"Weak")) {
		for (FACE = DB.FACE; FACE; FACE = FACE->next) {
			init_FDATA(FDATAL,FACE,'L');
			init_FDATA(FDATAR,FACE,'R');

			// Compute WL_fIL and WR_fIL (i.e. as seen from the (L)eft VOLUME)
			unsigned int IndFType = FDATAL->IndFType;
			unsigned int NfnI     = OPSL[IndFType]->NfnI;

			FDATAL->W_fIL = malloc(NfnI*Nvar * sizeof *(FDATAL->W_fIL)), // free
			FDATAR->W_fIL = malloc(NfnI*Nvar * sizeof *(FDATAR->W_fIL)); // free

			coef_to_values_fI(FDATAL,'W');
			compute_WR_fIL(FDATAR,FDATAL->W_fIL,FDATAR->W_fIL);

			// Compute numerical flux and its Jacobian as seen from the left VOLUME
			NFluxData->WL_fIL          = FDATAL->W_fIL;
			NFluxData->WR_fIL          = FDATAR->W_fIL;
			NFluxData->nFluxNum_fI     = malloc(NfnI*Neq      * sizeof *(NFluxData->nFluxNum_fI));     // free
			NFluxData->dnFluxNumdWL_fI = malloc(NfnI*Neq*Nvar * sizeof *(NFluxData->dnFluxNumdWL_fI)); // free
			NFluxData->dnFluxNumdWR_fI = malloc(NfnI*Neq*Nvar * sizeof *(NFluxData->dnFluxNumdWR_fI)); // free

			compute_numerical_flux(FDATAL,'I');
			add_Jacobian_scaling_FACE(FDATAL,'I');

			free(FDATAL->W_fIL);
			free(FDATAR->W_fIL);


			// Compute FACE RHS and LHS terms
			unsigned int NvnSL = OPSL[0]->NvnS,
			             NvnSR = OPSR[0]->NvnS;

			if (FACE->RHSIn)
				free(FACE->RHSIn);
			FACE->RHSIn  = calloc(NvnSL*Neq , sizeof *(FACE->RHSIn)); // keep

			if (FACE->RHSOut)
				free(FACE->RHSOut);
			FACE->RHSOut = calloc(NvnSR*Neq , sizeof *(FACE->RHSOut)); // keep

			FACE->LHSInIn   = calloc(NvnSL*NvnSL*Neq*Nvar , sizeof *(FACE->LHSInIn));   // keep
			FACE->LHSOutIn  = calloc(NvnSL*NvnSR*Neq*Nvar , sizeof *(FACE->LHSOutIn));  // keep
			FACE->LHSInOut  = calloc(NvnSR*NvnSL*Neq*Nvar , sizeof *(FACE->LHSInOut));  // keep
			FACE->LHSOutOut = calloc(NvnSR*NvnSR*Neq*Nvar , sizeof *(FACE->LHSOutOut)); // keep

			finalize_FACE_Inviscid_Weak(FDATAL,FDATAR,'L','E');
			finalize_FACE_Inviscid_Weak(FDATAL,FDATAR,'L','I');

			if (!FACE->Boundary) {
				// LHS
				finalize_FACE_Inviscid_Weak(FDATAL,FDATAR,'R','I');

				// RHS
				finalize_FACE_Inviscid_Weak(FDATAL,FDATAR,'R','E');
			}
			free(NFluxData->nFluxNum_fI);
			free(NFluxData->dnFluxNumdWL_fI);
			free(NFluxData->dnFluxNumdWR_fI);
		}
	} else if (strstr(Form,"Strong")) {
		EXIT_UNSUPPORTED;
	}

	for (size_t i = 0; i < 2; i++) {
		free(OPSL[i]);
		free(OPSR[i]);
	}

	free(FDATAL);
	free(FDATAR);
	free(NFluxData);
}
