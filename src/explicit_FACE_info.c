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
#include "array_swap.h"

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
 *		as if the mesh were conforming. As no L2 projection is required for the current approach, compute_FACE_RHS_EFE
 *		can in fact be called even for the case of non-conforming discretizations allowing for:
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

static void compute_FACE_RHS_EFE    (void);

void explicit_FACE_info(void)
{
	compute_FACE_RHS_EFE();
}

static void compute_FACE_RHS_EFE(void)
{
	// Initialize DB Parameters
	char         *Form            = DB.Form;
	unsigned int Nvar             = DB.Nvar,
	             Neq              = DB.Neq;

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

			// Compute WL_fI and WR_fIL
			unsigned int IndFType = FDATAL->IndFType;
			unsigned int NfnI     = OPSL[IndFType]->NfnI;

			double *WIn_fI = malloc(NfnI*Nvar * sizeof *WIn_fI); // free
			compute_W_fI(FDATAL,WIn_fI);

			double *WOut_fIIn = malloc(NfnI*Nvar * sizeof *WOut_fIIn); // free
			compute_WR_fIL(FDATAR,WIn_fI,WOut_fIIn);


			// Compute numerical flux as seen from the left VOLUME
			double *nFluxNum_fI = malloc(NfnI*Neq * sizeof *nFluxNum_fI); // free

			NFluxData->WL_fIL      = WIn_fI;
			NFluxData->WR_fIL      = WOut_fIIn;
			NFluxData->nFluxNum_fI = nFluxNum_fI;

			compute_numerical_flux(FDATAL,'E');
			add_Jacobian_scaling_FACE(FDATAL,'E');

			free(WIn_fI);
			free(WOut_fIIn);


			// Compute FACE RHS terms
			unsigned int NvnSL = OPSL[0]->NvnS,
						 NvnSR = OPSR[0]->NvnS;

			double *RHSL = calloc(NvnSL*Neq , sizeof *RHSL), // keep
				   *RHSR = calloc(NvnSR*Neq , sizeof *RHSR); // keep
			FACE->RHSIn  = RHSL;
			FACE->RHSOut = RHSR;

			finalize_FACE_Inviscid_Weak(FDATAL,'L','E');
			if (!FACE->Boundary) {
				swap_FACE_orientation(FDATAR,'E');
				finalize_FACE_Inviscid_Weak(FDATAR,'R','E');
			}
			free(nFluxNum_fI);
		}
	} else if (strstr(Form,"Strong")) {
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
