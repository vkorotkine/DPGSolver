// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "finalize_RHS.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "S_FACET.h"

#include "array_norm.h"
#include "matrix_functions.h"
#include "element_functions.h"
#include "exact_solutions.h"

/*
 *	Purpose:
 *		Finalize RHS term by summing VOLUME and FACET contributions and multiplying by the inverse mass matrix for
 *		explicit runs.
 *
 *	Comments:
 *		When the 'Collocated' option is enabled, the inverse VOLUME cubature weights are already included in the VV and
 *		FF operators used in the VOLUME and FACET info functions. Thus, only the inverse VOLUME element must be applied
 *		to RHS terms. When the option is disabled, the standard procedure of multiplication with the inverse mass matrix
 *		is carried out. The motivation for this is that a particular sum factorized operator becomes identity for the
 *		collocated scheme when this approach is adopted; see cases where diag = 2 in the VOLUME/FACET info routines.
 *
 *	Notation:
 *
 *	References:
 */

struct S_OPERATORS {
	unsigned int NvnI;
	double       *I_vG_vI, *I_Weak;
};

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME)
{
	// Standard datatypes
	unsigned int P, type, curved;
	struct S_ELEMENT *ELEMENT;

	P      = VOLUME->P;
	type   = VOLUME->type;
	curved = VOLUME->curved;

	ELEMENT = get_ELEMENT_type(type);

	if (!curved) {
		OPS->NvnI = ELEMENT->NvnIs[P];

		OPS->I_vG_vI = ELEMENT->I_vGs_vIs[1][P][0];
		OPS->I_Weak  = ELEMENT->Is_Weak_VV[P][P][0];
	} else {
		OPS->NvnI = ELEMENT->NvnIc[P];

		OPS->I_vG_vI = ELEMENT->I_vGc_vIc[P][P][0];
		OPS->I_Weak  = ELEMENT->Ic_Weak_VV[P][P][0];
	}
}

static void add_source(const struct S_VOLUME *VOLUME)
{
	// Initialize DB Parameters
	unsigned int d   = DB.d,
	             Neq = DB.Neq;

	// Standard datatypes
	unsigned int eq, n, NvnI;
	double    *XYZ_vI, *f_vI, *detJV_vI_ptr;

	struct S_OPERATORS *OPS;

	OPS = malloc(sizeof *OPS); // free

	init_ops(OPS,VOLUME);

	NvnI = OPS->NvnI;

	XYZ_vI = malloc(NvnI*d * sizeof *XYZ_vI); // free
	mm_CTN_d(NvnI,d,VOLUME->NvnG,OPS->I_vG_vI,VOLUME->XYZ,XYZ_vI);

	f_vI = malloc(NvnI*Neq * sizeof *f_vI); // free
	compute_source(NvnI,XYZ_vI,f_vI);
	free(XYZ_vI);

	for (eq = 0; eq < Neq; eq++) {
		detJV_vI_ptr = VOLUME->detJV_vI;
		for (n = 0; n < NvnI; n++)
			f_vI[eq*NvnI+n] *= *detJV_vI_ptr++;
	}

	mm_d(CBCM,CBT,CBNT,VOLUME->NvnS,Neq,NvnI,-1.0,1.0,OPS->I_Weak,f_vI,VOLUME->RHS);

	free(f_vI);

	free(OPS);
}

double finalize_RHS(void)
{
	// Initialize DB Parameters
	char         *SolverType   = DB.SolverType;
	unsigned int Neq           = DB.Neq,
	             SourcePresent = DB.SourcePresent,
	             Collocated    = DB.Collocated;

	// Standard datatypes
	unsigned int iMax, jMax,
	             NvnSIn, NvnSOut, Boundary;
	double       maxRHS, maxRHSV, *VRHSIn_ptr, *VRHSOut_ptr, *FRHSIn_ptr, *FRHSOut_ptr, *detJV_vI_ptr,
	             *RHS_Final;

	struct S_VOLUME    *VIn, *VOut, *VOLUME;
	struct S_FACET     *FACET;

	// silence
	NvnSIn = 0;

	for (FACET = DB.FACET; FACET; FACET = FACET->next) {
		VIn    = FACET->VIn;
		NvnSIn = VIn->NvnS;

		VOut    = FACET->VOut;
		NvnSOut = VOut->NvnS;

		VRHSIn_ptr  = VIn->RHS;
		VRHSOut_ptr = VOut->RHS;

		FRHSIn_ptr  = FACET->RHSIn;
		FRHSOut_ptr = FACET->RHSOut;

		Boundary = FACET->Boundary;
		for (iMax = Neq; iMax--; ) {
			for (jMax = NvnSIn; jMax--; )
				*VRHSIn_ptr++ += *FRHSIn_ptr++;
			if (!Boundary) {
				for (jMax = NvnSOut; jMax--; )
					*VRHSOut_ptr++ += *FRHSOut_ptr++;
			}
		}

		free(FACET->RHSIn);
		free(FACET->RHSOut);
	}

	// Add source contribution
	if (SourcePresent) {
		for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next)
			add_source(VOLUME);
	}

	// Compute maxRHS for convergence monitoring
	maxRHS = 0.0;
	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		// Compute maxRHS for convergence monitoring
		NvnSIn = VOLUME->NvnS;
		maxRHSV = array_norm_d(NvnSIn,VOLUME->RHS,"Inf");
		if (maxRHSV > maxRHS)
			maxRHS = maxRHSV;
	}

	// Add MInv contribution to RHS for explicit runs
	if (strstr(SolverType,"Explicit")) {
		if (SourcePresent)
			printf("Warning: Ensure that sources are being treated properly.\n");

		for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			// Add (remaining) MInv contribution to RHS
			if (Collocated) {
				VRHSIn_ptr = VOLUME->RHS;
				for (iMax = Neq; iMax--; ) {
					detJV_vI_ptr = VOLUME->detJV_vI;
					for (jMax = NvnSIn; jMax--; )
						*VRHSIn_ptr++ /= *detJV_vI_ptr++;
				}
			} else {
				RHS_Final = malloc(NvnSIn*Neq * sizeof *RHS_Final);

				mm_CTN_d(NvnSIn,Neq,NvnSIn,VOLUME->MInv,VOLUME->RHS,RHS_Final);
				free(VOLUME->RHS);
				VOLUME->RHS = RHS_Final;
			}

		}
	}
	return maxRHS;
}
