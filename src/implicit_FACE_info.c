// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "compute_FACE_RLHS_DG.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_FACE.h"

#include "solver.h"
#include "solver_functions.h"
#include "array_free.h"
#include "support.h"


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

static void set_memory_to_zero_FACEs  (const char imex_type);
static void compute_Inviscid_FACE_EFE (const char imex_type);
static void compute_Viscous_FACE_EFE  (const char imex_type);

void compute_FACE_RLHS_DG (const struct S_solver_info*const solver_info)
{
	if (solver_info->display)
		printf("F");

	set_memory_to_zero_FACEs(solver_info->imex_type);
	compute_Inviscid_FACE_EFE(solver_info->imex_type);
	compute_Viscous_FACE_EFE(solver_info->imex_type);
}

static void set_memory_to_zero_FACEs (const char imex_type)
{
	unsigned int const Nvar = DB.Nvar,
	                   Neq  = DB.Neq;

	for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
		unsigned int const NvnSL = FACE->VL->NvnS;
		set_to_zero_d(NvnSL*Neq,FACE->RHSL);
		if (imex_type == 'I')
			set_to_zero_d(NvnSL*NvnSL*Neq*Nvar,FACE->LHSLL);

		if (!FACE->Boundary) {
			unsigned int const NvnSR = FACE->VR->NvnS;
			set_to_zero_d(NvnSR*Neq,FACE->RHSR);
			if (imex_type == 'I') {
				set_to_zero_d(NvnSL*NvnSR*Neq*Nvar,FACE->LHSRL);
				set_to_zero_d(NvnSR*NvnSL*Neq*Nvar,FACE->LHSLR);
				set_to_zero_d(NvnSR*NvnSR*Neq*Nvar,FACE->LHSRR);
			}
		}
	}
}

static void compute_Inviscid_FACE_EFE (const char imex_type)
{
	if (!DB.Inviscid)
		return;

	struct S_OPERATORS_F *OPSL[2], *OPSR[2];

	struct S_FDATA *FDATAL = malloc(sizeof *FDATAL), // free
	               *FDATAR = malloc(sizeof *FDATAR); // free
	FDATAL->OPS = (struct S_OPERATORS_F const *const *) OPSL;
	FDATAR->OPS = (struct S_OPERATORS_F const *const *) OPSR;

	struct S_NUMERICALFLUX *NFLUXDATA = malloc(sizeof *NFLUXDATA); // free
	FDATAL->NFLUXDATA = NFLUXDATA;
	FDATAR->NFLUXDATA = NFLUXDATA;

	for (size_t i = 0; i < 2; i++) {
		OPSL[i] = malloc(sizeof *OPSL[i]); // free
		OPSR[i] = malloc(sizeof *OPSR[i]); // free
	}

	struct S_DATA *const DATA = malloc(sizeof *DATA); // free
	DATA->FDATAL    = FDATAL;
	DATA->FDATAR    = FDATAR;
	DATA->NFLUXDATA = NFLUXDATA;
	DATA->feature   = 'F';
	DATA->imex_type = imex_type;

	if (strstr(DB.Form,"Weak")) {
		for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
			init_FDATA(FDATAL,FACE,'L');
			init_FDATA(FDATAR,FACE,'R');

			// Compute WL_fIL and WR_fIL (i.e. as seen from the (L)eft VOLUME)
			manage_solver_memory(DATA,'A','W'); // free

			coef_to_values_fI(FDATAL,'W',imex_type);
			compute_WR_fIL(FDATAR,FDATAL->W_fIL,FDATAR->W_fIL);


			// Compute numerical flux and its Jacobian as seen from the left VOLUME
			NFLUXDATA->WL = FDATAL->W_fIL;
			NFLUXDATA->WR = FDATAR->W_fIL;
			manage_solver_memory(DATA,'A','I'); // free

			compute_numerical_flux(FDATAL,imex_type);
			add_Jacobian_scaling_FACE(FDATAL,imex_type,'W');

			manage_solver_memory(DATA,'F','W');


			// Compute FACE RHS and LHS terms
			finalize_FACE_Inviscid_Weak(FDATAL,FDATAR,NFLUXDATA->nFluxNum,NULL,'L','E','W');
			if (imex_type == 'I')
				finalize_FACE_Inviscid_Weak(FDATAL,FDATAR,NFLUXDATA->dnFluxNumdWL,NFLUXDATA->dnFluxNumdWR,'L','I','W');

			if (!FACE->Boundary) {
				finalize_FACE_Inviscid_Weak(FDATAL,FDATAR,NFLUXDATA->nFluxNum,NULL,'R','E','W');
				if (imex_type == 'I')
					finalize_FACE_Inviscid_Weak(FDATAL,FDATAR,NFLUXDATA->dnFluxNumdWL,NFLUXDATA->dnFluxNumdWR,'R','I','W');
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
	free(NFLUXDATA);
	free(DATA);
}

static void compute_Viscous_FACE_EFE (const char imex_type)
{
	if (!DB.Viscous)
		return;

	struct S_OPERATORS_F *OPSL[2], *OPSR[2];

	struct S_FDATA *FDATAL = malloc(sizeof *FDATAL), // free
	               *FDATAR = malloc(sizeof *FDATAR); // free
	FDATAL->OPS = (struct S_OPERATORS_F const *const *) OPSL;
	FDATAR->OPS = (struct S_OPERATORS_F const *const *) OPSR;

	struct S_NUMERICALFLUX *NFLUXDATA = malloc(sizeof *NFLUXDATA); // free
	FDATAL->NFLUXDATA = NFLUXDATA;
	FDATAR->NFLUXDATA = NFLUXDATA;

	struct S_DATA *const DATA = malloc(sizeof *DATA); // free
	DATA->FDATAL    = FDATAL;
	DATA->FDATAR    = FDATAR;
	DATA->NFLUXDATA = NFLUXDATA;
	DATA->feature   = 'F';
	DATA->imex_type = imex_type;

	for (size_t i = 0; i < 2; i++) {
		OPSL[i] = malloc(sizeof *OPSL[i]); // free
		OPSR[i] = malloc(sizeof *OPSR[i]); // free
	}

	if (strstr(DB.Form,"Weak")) {
		for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
			init_FDATA(FDATAL,FACE,'L');
			init_FDATA(FDATAR,FACE,'R');

			// Compute WL_fIL, WR_fIL, QpL_fIL, and QpR_fIL (i.e. as seen from the (L)eft VOLUME)
			manage_solver_memory(DATA,'A','W'); // free
			manage_solver_memory(DATA,'A','Q'); // free

			coef_to_values_fI(FDATAL,'W',imex_type);
			coef_to_values_fI(FDATAL,'Q',imex_type);
			compute_WR_QpR_fIL(FDATAR,FDATAL->W_fIL,FDATAR->W_fIL,
			                   (double const *const *const) FDATAL->Qp_fIL,FDATAR->Qp_fIL,'I');


			// Compute numerical flux and its Jacobian as seen from the left VOLUME
			NFLUXDATA->WL = FDATAL->W_fIL;
			NFLUXDATA->WR = FDATAR->W_fIL;
			manage_solver_memory(DATA,'A','V'); // free

			compute_numerical_flux_viscous(FDATAL,FDATAR,imex_type);
			manage_solver_memory(DATA,'F','W');
			manage_solver_memory(DATA,'F','Q');

			add_Jacobian_scaling_FACE(FDATAL,imex_type,'V');
			if (imex_type == 'I')
				add_Jacobian_scaling_FACE(FDATAL,'I','P');

			if (imex_type == 'I') {
				// FACE contribution to V(L/R)->LHS and related off-diagonal contributions from the VOLUME term
				finalize_VOLUME_LHSQF_Weak(FACE);
			}

			// Compute FACE RHS and LHS terms
			finalize_FACE_Viscous_Weak(FDATAL,FDATAR,NFLUXDATA->nFluxNum,NULL,'L','E','V');
			if (imex_type == 'I') {
				finalize_FACE_Viscous_Weak(FDATAL,FDATAR,NFLUXDATA->dnFluxNumdWL,NFLUXDATA->dnFluxNumdWR,'L','I','V');
				finalize_implicit_FACE_Q_Weak(FDATAL,FDATAR,'L');
			}

			if (!FACE->Boundary) {
				finalize_FACE_Viscous_Weak(FDATAL,FDATAR,NFLUXDATA->nFluxNum,NULL,'R','E','V');
				if (imex_type == 'I') {
					finalize_FACE_Viscous_Weak(FDATAL,FDATAR,NFLUXDATA->dnFluxNumdWL,NFLUXDATA->dnFluxNumdWR,'R','I','V');
					finalize_implicit_FACE_Q_Weak(FDATAL,FDATAR,'R');
				}
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
	free(NFLUXDATA);
	free(DATA);
}
