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
#include "array_free.h"


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

static void compute_Inviscid_FACE_EFE (void);
static void compute_Viscous_FACE_EFE (void);

void implicit_FACE_info(void)
{
	compute_Inviscid_FACE_EFE();
	compute_Viscous_FACE_EFE();
}

static void compute_Inviscid_FACE_EFE(void)
{
	unsigned int const d    = DB.d,
	                   Nvar = d+2,
	                   Neq  = d+2;

	struct S_OPERATORS_F *OPSL[2], *OPSR[2];

	struct S_FDATA *FDATAL = malloc(sizeof *FDATAL), // free
	               *FDATAR = malloc(sizeof *FDATAR); // free
	FDATAL->OPS = (struct S_OPERATORS_F const *const *) OPSL;
	FDATAR->OPS = (struct S_OPERATORS_F const *const *) OPSR;

	struct S_NumericalFlux *NFluxData = malloc(sizeof *NFluxData); // free
	FDATAL->NFluxData = NFluxData;
	FDATAR->NFluxData = NFluxData;

	for (size_t i = 0; i < 2; i++) {
		OPSL[i] = malloc(sizeof *OPSL[i]); // free
		OPSR[i] = malloc(sizeof *OPSR[i]); // free
	}

	struct S_DATA *const DATA = malloc(sizeof *DATA); // free
	DATA->FDATAL    = FDATAL;
	DATA->FDATAR    = FDATAR;
	DATA->NFluxData = NFluxData;
	DATA->feature   = 'F';
	DATA->imex_type = 'I';

	if (strstr(DB.Form,"Weak")) {
		for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
			init_FDATA(FDATAL,FACE,'L');
			init_FDATA(FDATAR,FACE,'R');

			// Compute WL_fIL and WR_fIL (i.e. as seen from the (L)eft VOLUME)
			manage_solver_memory(DATA,'A','W'); // free

			coef_to_values_fI(FDATAL,'W','I');
			compute_WR_fIL(FDATAR,FDATAL->W_fIL,FDATAR->W_fIL);


			// Compute numerical flux and its Jacobian as seen from the left VOLUME
			NFluxData->WL_fIL = FDATAL->W_fIL;
			NFluxData->WR_fIL = FDATAR->W_fIL;
			manage_solver_memory(DATA,'A','I'); // free

			compute_numerical_flux(FDATAL,'I');
			add_Jacobian_scaling_FACE(FDATAL,'I','W');

			manage_solver_memory(DATA,'F','W');


			// Compute FACE RHS and LHS terms
			unsigned int const NvnSL = OPSL[0]->NvnS,
			                   NvnSR = OPSR[0]->NvnS;

			memset(FACE->RHSL, 0.0,NvnSL*Neq            * sizeof *(FACE->RHSL));
			memset(FACE->LHSLL,0.0,NvnSL*NvnSL*Neq*Nvar * sizeof *(FACE->LHSLL));

			finalize_FACE_Inviscid_Weak(FDATAL,FDATAR,NFluxData->nFluxNum_fI,NULL,'L','E','W');
			finalize_FACE_Inviscid_Weak(FDATAL,FDATAR,NFluxData->dnFluxNumdWL_fI,NFluxData->dnFluxNumdWR_fI,'L','I','W');

			if (!FACE->Boundary) {
				memset(FACE->RHSR, 0.0,NvnSR*Neq            * sizeof *(FACE->RHSR));
				memset(FACE->LHSRL,0.0,NvnSL*NvnSR*Neq*Nvar * sizeof *(FACE->LHSRL));
				memset(FACE->LHSLR,0.0,NvnSR*NvnSL*Neq*Nvar * sizeof *(FACE->LHSLR));
				memset(FACE->LHSRR,0.0,NvnSR*NvnSR*Neq*Nvar * sizeof *(FACE->LHSRR));

				finalize_FACE_Inviscid_Weak(FDATAL,FDATAR,NFluxData->nFluxNum_fI,NULL,'R','E','W');
				finalize_FACE_Inviscid_Weak(FDATAL,FDATAR,NFluxData->dnFluxNumdWL_fI,NFluxData->dnFluxNumdWR_fI,'R','I','W');
			}
			manage_solver_memory(DATA,'F','I');
		}
	} else if (strstr(DB.Form,"Strong")) {
		EXIT_UNSUPPORTED;
	}

	for (size_t i = 0; i < 2; i++) {
		free(OPSL[i]);
		free(OPSR[i]);
	}

	free(FDATAL);
	free(FDATAR);
	free(NFluxData);
	free(DATA);
}

static void compute_Viscous_FACE_EFE(void)
{
	if (!DB.Viscous)
		return;

	struct S_OPERATORS_F *OPSL[2], *OPSR[2];

	struct S_FDATA *FDATAL = malloc(sizeof *FDATAL), // free
	               *FDATAR = malloc(sizeof *FDATAR); // free
	FDATAL->OPS = (struct S_OPERATORS_F const *const *) OPSL;
	FDATAR->OPS = (struct S_OPERATORS_F const *const *) OPSR;

	struct S_NumericalFlux *NFluxData = malloc(sizeof *NFluxData); // free
	FDATAL->NFluxData = NFluxData;
	FDATAR->NFluxData = NFluxData;

	struct S_DATA *const DATA = malloc(sizeof *DATA); // free
	DATA->FDATAL    = FDATAL;
	DATA->FDATAR    = FDATAR;
	DATA->NFluxData = NFluxData;
	DATA->feature   = 'F';
	DATA->imex_type = 'I';

	for (size_t i = 0; i < 2; i++) {
		OPSL[i] = malloc(sizeof *OPSL[i]); // free
		OPSR[i] = malloc(sizeof *OPSR[i]); // free
	}

	if (strstr(DB.Form,"Weak")) {
		for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
			// FACE contribution to V(L/R)->LHS and related off-diagonal contributions from the VOLUME term
			finalize_VOLUME_LHSQF_Weak(FACE);

			init_FDATA(FDATAL,FACE,'L');
			init_FDATA(FDATAR,FACE,'R');

			// Compute WL_fIL, WR_fIL, QpL_fIL, and QpR_fIL (i.e. as seen from the (L)eft VOLUME)
			manage_solver_memory(DATA,'A','W'); // free
			manage_solver_memory(DATA,'A','Q'); // free

			coef_to_values_fI(FDATAL,'W','I');
			coef_to_values_fI(FDATAL,'Q','I');
			compute_WR_QpR_fIL(FDATAR,FDATAL->W_fIL,FDATAR->W_fIL,
			                   (double const *const *const) FDATAL->Qp_fIL,FDATAR->Qp_fIL,'I');


			// Compute numerical flux and its Jacobian as seen from the left VOLUME
			NFluxData->WL_fIL = FDATAL->W_fIL;
			NFluxData->WR_fIL = FDATAR->W_fIL;
			manage_solver_memory(DATA,'A','V'); // free

			compute_numerical_flux_viscous(FDATAL,FDATAR,'I');
			add_Jacobian_scaling_FACE(FDATAL,'I','V');
			add_Jacobian_scaling_FACE(FDATAL,'I','P');

			manage_solver_memory(DATA,'F','W');
			manage_solver_memory(DATA,'F','Q');


			// Compute FACE RHS and LHS terms
			finalize_FACE_Viscous_Weak(FDATAL,FDATAR,NFluxData->nFluxViscNum_fI,NULL,'L','E','V');
			finalize_FACE_Viscous_Weak(FDATAL,FDATAR,NFluxData->dnFluxViscNumdWL_fI,NFluxData->dnFluxViscNumdWR_fI,'L','I','V');
			finalize_implicit_FACE_Q_Weak(FDATAL,FDATAR,'L');

			if (!FACE->Boundary) {
				finalize_FACE_Viscous_Weak(FDATAL,FDATAR,NFluxData->nFluxViscNum_fI,NULL,'R','E','V');
				finalize_FACE_Viscous_Weak(FDATAL,FDATAR,NFluxData->dnFluxViscNumdWL_fI,NFluxData->dnFluxViscNumdWR_fI,'R','I','V');
				finalize_implicit_FACE_Q_Weak(FDATAL,FDATAR,'R');
			}
			manage_solver_memory(DATA,'F','V');
		}
	} else if (strstr(DB.Form,"Strong")) {
		EXIT_UNSUPPORTED;
	}

	for (size_t i = 0; i < 2; i++) {
		free(OPSL[i]);
		free(OPSR[i]);
	}

	free(FDATAL);
	free(FDATAR);
	free(NFluxData);
	free(DATA);
}
