// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "explicit_FACE_info.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_FACE.h"

#include "solver_functions.h"
#include "fluxes_structs.h"
#include "array_free.h"


/*
 *	Purpose:
 *		Evaluate the FACE contributions to the RHS term.
 *
 *	Comments:
 *
 *		When adaptation is enabled, the method used here to lift the FACE information to the VOLUME remains quite
 *		similar to that of the conforming case. The normal numerical flux is computed on each of the FACEs of the mesh
 *		and is used directly to form RHSL/RHSR terms to be added to the appropriate VOLUME. This is in contrast to the
 *		traditional mortar element method (based on my current understanding), which first uses an L2 projection of the
 *		normal numerical flux of all non-conforming FACEs to standard VOLUME FACEs and then computes RHSL/RHSR exactly
 *		as if the mesh were conforming. As no L2 projection is required for the current approach,
 *		compute_Inviscid_FACE_RHS_EFE can in fact be called even for the case of non-conforming discretizations allowing
 *		for:
 *			1) Reduced cost because the interpolation from FACE solution to FACE cubature nodes is not required;
 *			2) Reduced aliasing because the normal numerical flux can be computed at the cubature nodes.
 *
 *		Based on this discussion, it is unclear why this approach is not adopted instead of the traditional mortar
 *		method (Kopriva(1996)) as this alternative seems to satisfy both the conservation and outflow condition
 *		requirements which motivated the use of the mortar element method.
 *
 *		Vectorization is more involved for FACE terms as there are many more possible combinations than for VOLUMEs.
 *		Given that the VOLUME vectorization seems not to have a significant impact on performance, the FACE
 *		vectorization may not be pursued. (ToBeModified)
 *
 *	Notation:
 *		Qp : (p)artially corrected weak solution gradients (Q).
 *
 *	References:
 *		Kopriva(1996)-A_Conservative_Staggered-Grid_Chebyshev_Multidomain_Method_for_Compressible_Flows_II._A_Semi-Structured_Method
 */

static void compute_Inviscid_FACE_RHS_EFE (void);
static void compute_Viscous_FACE_RHS_EFE  (void);

void explicit_FACE_info(void)
{
	compute_Inviscid_FACE_RHS_EFE();
	compute_Viscous_FACE_RHS_EFE();
}

static void compute_Inviscid_FACE_RHS_EFE(void)
{
	unsigned int const Neq = DB.Neq;

	struct S_OPERATORS_F *OPSL[2], *OPSR[2];

	struct S_FDATA *const FDATAL = malloc(sizeof *FDATAL), // free
	               *const FDATAR = malloc(sizeof *FDATAR); // free
	FDATAL->OPS = (struct S_OPERATORS_F const *const *) OPSL;
	FDATAR->OPS = (struct S_OPERATORS_F const *const *) OPSR;

	struct S_NUMERICALFLUX *const NFLUXDATA = malloc(sizeof *NFLUXDATA); // free
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
	DATA->imex_type = 'E';

	if (strstr(DB.Form,"Weak")) {
		for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
			init_FDATA(FDATAL,FACE,'L');
			init_FDATA(FDATAR,FACE,'R');

			// Compute WL_fIL and WR_fIL (i.e. as seen from the (L)eft VOLUME)
			manage_solver_memory(DATA,'A','W'); // free

			coef_to_values_fI(FDATAL,'W','E');
			compute_WR_fIL(FDATAR,FDATAL->W_fIL,FDATAR->W_fIL);


			// Compute numerical flux as seen from the left VOLUME
			NFLUXDATA->WL = FDATAL->W_fIL;
			NFLUXDATA->WR = FDATAR->W_fIL;
			manage_solver_memory(DATA,'A','I'); // free

			compute_numerical_flux(FDATAL,'E');
			add_Jacobian_scaling_FACE(FDATAL,'E','W');

			manage_solver_memory(DATA,'F','W');


			// Compute FACE RHS terms
			unsigned int const NvnSL = OPSL[0]->NvnS,
			                   NvnSR = OPSR[0]->NvnS;

			memset(FACE->RHSL,0.0,NvnSL*Neq * sizeof *(FACE->RHSL));
			finalize_FACE_Inviscid_Weak(FDATAL,FDATAR,NFLUXDATA->nFluxNum,NULL,'L','E','W');

			if (!FACE->Boundary) {
				memset(FACE->RHSR,0.0,NvnSR*Neq * sizeof *(FACE->RHSR));
				finalize_FACE_Inviscid_Weak(FDATAL,FDATAR,NFLUXDATA->nFluxNum,NULL,'R','E','W');
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

static void compute_Viscous_FACE_RHS_EFE(void)
{
	if (!DB.Viscous)
		return;

	struct S_OPERATORS_F *OPSL[2], *OPSR[2];

	struct S_FDATA *const FDATAL = malloc(sizeof *FDATAL), // free
	               *const FDATAR = malloc(sizeof *FDATAR); // free
	FDATAL->OPS = (struct S_OPERATORS_F const *const *) OPSL;
	FDATAR->OPS = (struct S_OPERATORS_F const *const *) OPSR;

	struct S_NUMERICALFLUX *const NFLUXDATA = malloc(sizeof *NFLUXDATA); // free
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

	if (strstr(DB.Form,"Weak")) {
		for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
			init_FDATA(FDATAL,FACE,'L');
			init_FDATA(FDATAR,FACE,'R');

			// Compute WL_fIL and WR_fIL and their gradients (as seen from the (L)eft VOLUME)
			manage_solver_memory(DATA,'A','W'); // free
			manage_solver_memory(DATA,'A','Q'); // free

			coef_to_values_fI(FDATAL,'W','E');
			coef_to_values_fI(FDATAL,'Q','E');
			compute_WR_QpR_fIL(FDATAR,FDATAL->W_fIL,FDATAR->W_fIL,
			                   (double const *const *const) FDATAL->Qp_fIL,FDATAR->Qp_fIL,'E');


			// Compute numerical flux as seen from the left VOLUME
			NFLUXDATA->WL = FDATAL->W_fIL;
			NFLUXDATA->WR = FDATAR->W_fIL;
			manage_solver_memory(DATA,'A','V'); // free

			compute_numerical_flux_viscous(FDATAL,FDATAR,'E');
			add_Jacobian_scaling_FACE(FDATAL,'E','V');

			manage_solver_memory(DATA,'F','W');
			manage_solver_memory(DATA,'F','Q');


			// Compute FACE RHS terms
			finalize_FACE_Viscous_Weak(FDATAL,FDATAR,NFLUXDATA->nFluxNum,NULL,'L','E','V');
			if (!FACE->Boundary)
				finalize_FACE_Viscous_Weak(FDATAL,FDATAR,NFLUXDATA->nFluxNum,NULL,'R','E','V');

			manage_solver_memory(DATA,'F','V');
		}
	} else if (strstr(DB.Form,"Strong")) {
		// Note that the viscous flux is negated.
		EXIT_UNSUPPORTED;
	}

	for (size_t i = 0; i < 2; i++) {
		free(OPSL[i]);
		free(OPSR[i]);
	}
	free(NFLUXDATA);
	free(FDATAL);
	free(FDATAR);
	free(DATA);
}
