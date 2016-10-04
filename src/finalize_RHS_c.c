// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "finalize_RHS_c.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>

#include "petscvec.h"

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
 *		Identical to finalize_RHS using complex variables (for complex step verification).
 *
 *	Comments:
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

static void compute_source_c(const unsigned int Nn, double *XYZ, double complex *source)
{
	unsigned int   n;
	double         *source_d;

	source_d = malloc(Nn * sizeof *source_d); // free
	compute_source(Nn,XYZ,source_d);

	for (n = 0; n < Nn; n++)
		source[n] = source_d[n];
	free(source_d);
}

static void add_source(const struct S_VOLUME *VOLUME)
{
	// Initialize DB Parameters
	unsigned int d   = DB.d,
	             Neq = DB.Neq;

	// Standard datatypes
	unsigned int   eq, n, NvnI;
	double         *XYZ_vI, *detJV_vI_ptr;
	double complex *f_vI;

	struct S_OPERATORS *OPS;

	OPS = malloc(sizeof *OPS); // free

	init_ops(OPS,VOLUME);

	NvnI = OPS->NvnI;

	XYZ_vI = malloc(NvnI*d * sizeof *XYZ_vI); // free
	mm_CTN_d(NvnI,d,VOLUME->NvnG,OPS->I_vG_vI,VOLUME->XYZ,XYZ_vI);

	f_vI = malloc(NvnI*Neq * sizeof *f_vI); // free
	compute_source_c(NvnI,XYZ_vI,f_vI);
	free(XYZ_vI);

	for (eq = 0; eq < Neq; eq++) {
		detJV_vI_ptr = VOLUME->detJV_vI;
		for (n = 0; n < NvnI; n++)
			f_vI[eq*NvnI+n] *= *detJV_vI_ptr++;
	}

	mm_dcc(CBCM,CBT,CBNT,VOLUME->NvnS,Neq,NvnI,-1.0,1.0,OPS->I_Weak,f_vI,VOLUME->RHS_c);

	free(f_vI);

	free(OPS);
}

void finalize_RHS_c(void)
{
	// Initialize DB Parameters
	char         *SolverType = DB.SolverType;
	unsigned int Neq           = DB.Neq,
	             SourcePresent = DB.SourcePresent;

	// Standard datatypes
	unsigned int   iMax, jMax, NvnSIn, NvnSOut, Boundary;
	double complex *VRHSIn_ptr, *VRHSOut_ptr, *FRHSIn_ptr, *FRHSOut_ptr;

	struct S_VOLUME  *VOLUME, *VIn, *VOut;
	struct S_FACET   *FACET;

	for (FACET = DB.FACET; FACET; FACET = FACET->next) {
		VIn    = FACET->VIn;
		NvnSIn = VIn->NvnS;

		VOut    = FACET->VOut;
		NvnSOut = VOut->NvnS;

		VRHSIn_ptr  = VIn->RHS_c;
		VRHSOut_ptr = VOut->RHS_c;

		FRHSIn_ptr  = FACET->RHSIn_c;
		FRHSOut_ptr = FACET->RHSOut_c;

		Boundary = FACET->Boundary;
		for (iMax = Neq; iMax--; ) {
			for (jMax = NvnSIn; jMax--; )
				*VRHSIn_ptr++ += *FRHSIn_ptr++;
			if (!Boundary) {
				for (jMax = NvnSOut; jMax--; )
					*VRHSOut_ptr++ += *FRHSOut_ptr++;
			}
		}

		free(FACET->RHSIn_c);  FACET->RHSIn_c  = NULL;
		free(FACET->RHSOut_c); FACET->RHSOut_c = NULL;
	}

	// Add source contribution
	if (SourcePresent) {
		for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next)
			add_source(VOLUME);
	}

	if (strstr(SolverType,"Explicit"))
		printf("Error: This functions should only be used for implicit linearization verification.\n"), EXIT_MSG;
}
