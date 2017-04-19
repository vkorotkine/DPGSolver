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
	unsigned int const Nvar = DB.Nvar,
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

			// Compute WL_fIL and WR_fIL (i.e. as seen from the (L)eft VOLUME)
			unsigned int const IndFType = FDATAL->IndFType,
			                   NfnI     = OPSL[IndFType]->NfnI;

			FDATAL->W_fIL = malloc(NfnI*Nvar * sizeof *(FDATAL->W_fIL)), // free
			FDATAR->W_fIL = malloc(NfnI*Nvar * sizeof *(FDATAR->W_fIL)); // free

			coef_to_values_fI(FDATAL,'W');
			compute_WR_fIL(FDATAR,FDATAL->W_fIL,FDATAR->W_fIL);


			// Compute numerical flux as seen from the left VOLUME
			NFluxData->WL_fIL      = FDATAL->W_fIL;
			NFluxData->WR_fIL      = FDATAR->W_fIL;
			NFluxData->nFluxNum_fI = malloc(NfnI*Neq * sizeof *(NFluxData->nFluxNum_fI)); // free

			compute_numerical_flux(FDATAL,'E');
			add_Jacobian_scaling_FACE(FDATAL,'E','W');

			free(FDATAL->W_fIL);
			free(FDATAR->W_fIL);


			// Compute FACE RHS terms
			unsigned int const NvnSL = OPSL[0]->NvnS,
			                   NvnSR = OPSR[0]->NvnS;

			memset(FACE->RHSIn,0.0,NvnSL*Neq * sizeof *(FACE->RHSIn));
			finalize_FACE_Inviscid_Weak(FDATAL,FDATAR,NFluxData->nFluxNum_fI,NULL,'L','E','W');

			if (!FACE->Boundary) {
				memset(FACE->RHSOut,0.0,NvnSR*Neq * sizeof *(FACE->RHSOut));
				finalize_FACE_Inviscid_Weak(FDATAL,FDATAR,NFluxData->nFluxNum_fI,NULL,'R','E','W');
			}

			free(NFluxData->nFluxNum_fI);
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
	/*
	 *	Comments:
	 *		It is currently hard-coded that the viscous flux has no dependence on the weak gradient. This is done so
	 *		that element coupling is restricted to VOLUMEs and their immediate neighbours, maintaining the previous data
	 *		structure when the linearization is computed.
	 */

	if (!DB.Viscous)
		return;

	unsigned int const d    = DB.d,
	                   Nvar = DB.Nvar,
	                   Neq  = DB.Neq;

	struct S_OPERATORS_F *OPSL[2], *OPSR[2];
	struct S_FACE     *FACE;

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
		for (FACE = DB.FACE; FACE; FACE = FACE->next) {
			init_FDATA(FDATAL,FACE,'L');
			init_FDATA(FDATAR,FACE,'R');

			// Compute WL_fIL and WR_fIL and their gradients (as seen from the (L)eft VOLUME)
			unsigned int const IndFType = FDATAL->IndFType,
			                   NfnI     = OPSL[IndFType]->NfnI;

			FDATAL->W_fIL = malloc(NfnI*Nvar * sizeof *(FDATAL->W_fIL)), // free
			FDATAR->W_fIL = malloc(NfnI*Nvar * sizeof *(FDATAR->W_fIL)); // free

			double **const GradWL_fIL = malloc(d * sizeof *GradWL_fIL), // free
			       **const GradWR_fIL = malloc(d * sizeof *GradWR_fIL); // free

			for (size_t dim = 0; dim < d; dim++) {
				GradWL_fIL[dim] = malloc(NfnI*Nvar * sizeof *GradWL_fIL[dim]); // free
				GradWR_fIL[dim] = malloc(NfnI*Nvar * sizeof *GradWR_fIL[dim]); // free
			}

			FDATAL->GradW_fIL = GradWL_fIL;
			FDATAR->GradW_fIL = GradWR_fIL;

			coef_to_values_fI(FDATAL,'W');
			coef_to_values_fI(FDATAL,'Q');
			compute_WR_GradWR_fIL(FDATAR,FDATAL->W_fIL,FDATAR->W_fIL,(double const *const *const) GradWL_fIL,GradWR_fIL);


			// Compute numerical flux as seen from the left VOLUME
			NFluxData->WL_fIL          = FDATAL->W_fIL;
			NFluxData->WR_fIL          = FDATAR->W_fIL;
			NFluxData->nFluxViscNum_fI = malloc(NfnI*Neq * sizeof *(NFluxData->nFluxViscNum_fI)); // free

			compute_numerical_flux_viscous(FDATAL,FDATAR,'E');
			add_Jacobian_scaling_FACE(FDATAL,'E','V');

			free(FDATAL->W_fIL);
			free(FDATAR->W_fIL);
			array_free2_d(d,GradWL_fIL);
			array_free2_d(d,GradWR_fIL);


			// Compute FACE RHS terms
			finalize_FACE_Viscous_Weak(FDATAL,FDATAR,NFluxData->nFluxViscNum_fI,NULL,'L','E','V');
			if (!FACE->Boundary)
				finalize_FACE_Viscous_Weak(FDATAL,FDATAR,NFluxData->nFluxViscNum_fI,NULL,'R','E','V');

			free(NFluxData->nFluxViscNum_fI);
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
