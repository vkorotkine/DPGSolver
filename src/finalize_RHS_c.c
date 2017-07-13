// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

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
#include "S_FACE.h"

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
	double       *w_vI, *I_vG_vI, *ChiS_vI;
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
	} else {
		OPS->NvnI = ELEMENT->NvnIc[P];

		OPS->w_vI = ELEMENT->w_vIc[P];

		OPS->I_vG_vI = ELEMENT->I_vGc_vIc[P][P][0];
		OPS->ChiS_vI = ELEMENT->ChiS_vIc[P][P][0];
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
	double         *XYZ_vI, *detJV_vI, *w_vI;
	double complex *f_vI;

	struct S_OPERATORS *OPS;

	OPS = malloc(sizeof *OPS); // free

	init_ops(OPS,VOLUME);

	NvnI = OPS->NvnI;

	w_vI = OPS->w_vI;

	XYZ_vI = malloc(NvnI*d * sizeof *XYZ_vI); // free
	mm_CTN_d(NvnI,d,VOLUME->NvnG,OPS->I_vG_vI,VOLUME->XYZ,XYZ_vI);

	f_vI = malloc(NvnI*Neq * sizeof *f_vI); // free
	compute_source_c(NvnI,XYZ_vI,f_vI);
	free(XYZ_vI);

	detJV_vI = VOLUME->detJV_vI;
	for (eq = 0; eq < Neq; eq++) {
		for (n = 0; n < NvnI; n++)
			f_vI[eq*NvnI+n] *= w_vI[n]*detJV_vI[n];
	}

	mm_dcc(CBCM,CBNT,CBNT,VOLUME->NvnS,Neq,NvnI,-1.0,1.0,OPS->ChiS_vI,f_vI,VOLUME->RHS_c);

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
	unsigned int   iMax, jMax, NvnSL, NvnSR, Boundary;
	double complex *VRHSL_ptr, *VRHSR_ptr, *FRHSL_ptr, *FRHSR_ptr;

	struct S_VOLUME  *VOLUME, *VL, *VR;
	struct S_FACE   *FACE;

	for (FACE = DB.FACE; FACE; FACE = FACE->next) {
		VL    = FACE->VL;
		NvnSL = VL->NvnS;

		VR    = FACE->VR;
		NvnSR = VR->NvnS;

		VRHSL_ptr = VL->RHS_c;
		VRHSR_ptr = VR->RHS_c;

		FRHSL_ptr = FACE->RHSL_c;
		FRHSR_ptr = FACE->RHSR_c;

		Boundary = FACE->Boundary;
		for (iMax = Neq; iMax--; ) {
			for (jMax = NvnSL; jMax--; )
				*VRHSL_ptr++ += *FRHSL_ptr++;
			if (!Boundary) {
				for (jMax = NvnSR; jMax--; )
					*VRHSR_ptr++ += *FRHSR_ptr++;
			}
		}

		// ToBeDeleted (2 lines below)
		free(FACE->RHSL_c); FACE->RHSL_c = NULL;
		free(FACE->RHSR_c); FACE->RHSR_c = NULL;
	}

	// Add source contribution
	if (SourcePresent) {
		for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next)
			add_source(VOLUME);
	}

	if (strstr(SolverType,"Explicit"))
		printf("Error: This functions should only be used for implicit linearization verification.\n"), EXIT_MSG;
}
