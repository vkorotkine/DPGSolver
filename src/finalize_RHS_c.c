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
#include "S_VOLUME.h"
#include "S_FACET.h"

#include "array_norm.h"
#include "matrix_functions.h"

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

void finalize_RHS_c(void)
{
	// Initialize DB Parameters
	char         *SolverType = DB.SolverType;
	unsigned int Neq         = DB.Neq;

	// Standard datatypes
	unsigned int   iMax, jMax, NvnSIn, NvnSOut, Boundary;
	double complex *VRHSIn_ptr, *VRHSOut_ptr, *FRHSIn_ptr, *FRHSOut_ptr;

	struct S_VOLUME  *VIn, *VOut;
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

	if (strstr(SolverType,"Explicit"))
		printf("Error: This functions should only be used for implicit linearization verification.\n"), EXIT_MSG;
}
