// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "solver_Poisson_c.h"

#include <stdlib.h>
#include <stdio.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_VOLUME.h"
#include "S_FACE.h"

#include "solver_functions.h"
#include "solver_functions_c.h"
#include "matrix_functions.h"

/*
 *	Purpose:
 *		Compute RHS for Poisson solver using complex variables (for linearization testing).
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void compute_What_VOLUME_c(void)
{
	unsigned int const d   = DB.d,
	                   Neq = 1;

	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		unsigned int const NvnS = VOLUME->NvnS;

		// Compute RHS
		if (VOLUME->RHS_c)
			free(VOLUME->RHS_c);
		VOLUME->RHS_c = calloc(NvnS*Neq , sizeof *(VOLUME->RHS_c)); // keep

		for (size_t dim = 0; dim < d; dim++)
			mm_dcc(CBCM,CBT,CBNT,NvnS,1,NvnS,1.0,1.0,VOLUME->LHSQ[dim],VOLUME->Qhat_c[dim],VOLUME->RHS_c);
	}
}

void compute_What_FACE_c()
{
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
	DATA->imex_type = 'E';

	for (size_t i = 0; i < 2; i++) {
		OPSL[i] = malloc(sizeof *OPSL[i]); // free
		OPSR[i] = malloc(sizeof *OPSR[i]); // free
	}

	for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
		init_FDATA(FDATAL,FACE,'L');
		init_FDATA(FDATAR,FACE,'R');

		// Compute WL_fIL, WR_fIL, QpL_fIL, and QpR_fIL (i.e. as seen from the (L)eft VOLUME)
		manage_solver_memory_c(DATA,'A','W'); // free
		manage_solver_memory_c(DATA,'A','Q'); // free

		coef_to_values_fI_c(FDATAL,'W','E');
		coef_to_values_fI_c(FDATAL,'Q','E');
		compute_WR_QpR_fIL_c(FDATAR,FDATAL->W_fIL_c,FDATAR->W_fIL_c,
		                     (double complex const *const *const) FDATAL->Qp_fIL_c,FDATAR->Qp_fIL_c,'E');

		// Compute numerical flux and its Jacobian as seen from the left VOLUME
		NFLUXDATA->WL_c = FDATAL->W_fIL_c;
		NFLUXDATA->WR_c = FDATAR->W_fIL_c;
		manage_solver_memory_c(DATA,'A','P'); // free

		compute_numerical_flux_viscous_c(FDATAL,FDATAR,'E');
		manage_solver_memory_c(DATA,'F','W');
		manage_solver_memory_c(DATA,'F','Q');

		add_Jacobian_scaling_FACE_c(FDATAL,'E','V');

		// Negate the flux and flux Jacobian terms
		unsigned int const IndFType = FDATAL->IndFType,
		                   NfnI     = OPSL[IndFType]->NfnI;

		for (size_t n = 0; n < NfnI; n++)
			NFLUXDATA->nFluxNum_c[n] *= -1.0;

		// Finalize FACE RHS terms
		unsigned int const NvnSL = OPSL[0]->NvnS,
		                   NvnSR = OPSR[0]->NvnS;

		// VOLUME contribution taken care of in compute_What_VOLUME_c.
		if (FACE->RHSL_c)
			free(FACE->RHSL_c);
		FACE->RHSL_c = calloc(NvnSL , sizeof *(FACE->RHSL_c)); // keep

		if (!FACE->Boundary) {
			if (FACE->RHSR_c)
				free(FACE->RHSR_c);
			FACE->RHSR_c = calloc(NvnSR , sizeof *(FACE->RHSR_c)); // keep
		}

		// Add FACE contributions to RHS

		// Interior FACE
		finalize_FACE_Viscous_Weak_c(FDATAL,FDATAR,NFLUXDATA->nFluxNum_c,'L','E','V');

		// Exterior FACE
		if (!FACE->Boundary) {
			finalize_FACE_Viscous_Weak_c(FDATAL,FDATAR,NFLUXDATA->nFluxNum_c,'R','E','V');
		}
		manage_solver_memory_c(DATA,'F','P');
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
