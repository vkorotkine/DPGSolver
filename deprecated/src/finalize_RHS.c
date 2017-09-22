// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "finalize_RHS.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "S_FACE.h"

#include "array_norm.h"
#include "matrix_functions.h"
#include "element_functions.h"
#include "exact_solutions.h"

/*
 *	Purpose:
 *		Finalize RHS term by summing VOLUME and FACE contributions and multiplying by the inverse mass matrix for
 *		explicit runs.
 *
 *	Comments:
 *		When the 'Collocated' option is enabled, the inverse VOLUME cubature weights are already included in the VV and
 *		FV operators used in the VOLUME and FACE info functions. Thus, only the inverse VOLUME element must be applied
 *		to RHS terms. When the option is disabled, the standard procedure of multiplication with the inverse mass matrix
 *		is carried out. The motivation for this is that a particular sum factorized operator becomes identity for the
 *		collocated scheme when this approach is adopted; see cases where diag = 2 in the VOLUME/FACE info routines.
 *
 *		Likely delete EFE = 0 option for source computation; implemented only for comparison with Hesthaven.
 *		(ToBeDeleted)
 *
 *	Notation:
 *
 *	References:
 */

struct S_OPERATORS {
	unsigned int NvnI;
	double       *w_vI, *I_vG_vI, *ChiS_vI, *I_Weak_VV;
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

		OPS->w_vI = ELEMENT->w_vIs[P];

		OPS->I_vG_vI = ELEMENT->I_vGs_vIs[1][P][0];
		OPS->ChiS_vI = ELEMENT->ChiS_vIs[P][P][0];

		OPS->I_Weak_VV = ELEMENT->Is_Weak_VV[P][P][0];
	} else {
		OPS->NvnI = ELEMENT->NvnIc[P];

		OPS->w_vI = ELEMENT->w_vIc[P];

		OPS->I_vG_vI = ELEMENT->I_vGc_vIc[P][P][0];
		OPS->ChiS_vI = ELEMENT->ChiS_vIc[P][P][0];

		OPS->I_Weak_VV = ELEMENT->Ic_Weak_VV[P][P][0];
	}
}

static void add_source(const struct S_VOLUME *VOLUME)
{
// Potential different treatment for Poisson and Euler equations -> Investigate (ToBeDeleted).

	// Initialize DB Parameters
	unsigned int d   = DB.d,
	             Neq = DB.Neq;

	// Standard datatypes
	unsigned int eq, n, NvnI;
	double       *XYZ_vI, *f_vI, *detJV_vI;

	struct S_OPERATORS *OPS;

	OPS = malloc(sizeof *OPS); // free

	init_ops(OPS,VOLUME);

	NvnI = OPS->NvnI;

	XYZ_vI = malloc(NvnI*d * sizeof *XYZ_vI); // free
	mm_CTN_d(NvnI,d,VOLUME->NvnG,OPS->I_vG_vI,VOLUME->XYZ,XYZ_vI);

	f_vI = malloc(NvnI*Neq * sizeof *f_vI); // free
	compute_source(NvnI,XYZ_vI,f_vI);
	free(XYZ_vI);

	detJV_vI = VOLUME->detJV_vI;
	for (eq = 0; eq < Neq; eq++) {
		for (n = 0; n < NvnI; n++)
			f_vI[eq*NvnI+n] *= detJV_vI[n];
	}

	mm_d(CBCM,CBT,CBNT,VOLUME->NvnS,Neq,NvnI,1.0,1.0,OPS->I_Weak_VV,f_vI,VOLUME->RHS);
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
	             NvnSL, NvnSR, Boundary;
	double       maxRHS, maxRHSV, *VRHSL_ptr, *VRHSR_ptr, *FRHSL_ptr, *FRHSR_ptr, *detJV_vI_ptr,
	             *RHS_Final;

	struct S_VOLUME    *VL, *VR, *VOLUME;
	struct S_FACE     *FACE;

	// silence
	NvnSL = 0;

	for (FACE = DB.FACE; FACE; FACE = FACE->next) {
		VL    = FACE->VL;
		NvnSL = VL->NvnS;

		VR    = FACE->VR;
		NvnSR = VR->NvnS;

		VRHSL_ptr = VL->RHS;
		VRHSR_ptr = VR->RHS;

		FRHSL_ptr = FACE->RHSL;
		FRHSR_ptr = FACE->RHSR;

		Boundary = FACE->Boundary;
		for (iMax = Neq; iMax--; ) {
			for (jMax = NvnSL; jMax--; )
				*VRHSL_ptr++ += *FRHSL_ptr++;
			if (!Boundary) {
				for (jMax = NvnSR; jMax--; )
					*VRHSR_ptr++ += *FRHSR_ptr++;
			}
		}
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
		NvnSL = VOLUME->NvnS;
		maxRHSV = array_norm_d(NvnSL,VOLUME->RHS,"Inf");
		if (maxRHSV > maxRHS)
			maxRHS = maxRHSV;
	}

	// Add MInv contribution to RHS for explicit runs
	if (strstr(SolverType,"Explicit")) {
		for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			// Add (remaining) MInv contribution to RHS
			if (Collocated) {
				VRHSL_ptr = VOLUME->RHS;
				for (iMax = Neq; iMax--; ) {
					detJV_vI_ptr = VOLUME->detJV_vI;
					for (jMax = NvnSL; jMax--; )
						*VRHSL_ptr++ /= *detJV_vI_ptr++;
				}
			} else {
				RHS_Final = malloc(NvnSL*Neq * sizeof *RHS_Final);

				mm_CTN_d(NvnSL,Neq,NvnSL,VOLUME->MInv,VOLUME->RHS,RHS_Final);
				free(VOLUME->RHS);
				VOLUME->RHS = RHS_Final;
			}

		}
	}
	return maxRHS;
}
