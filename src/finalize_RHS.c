// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Finalize RHS term by summing VOLUME and FACET contributions and multiplying by the inverse mass matrix for
 *		explicit runs.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

double finalize_RHS(void)
{
	// Initialize DB Parameters
	char         *SolverType = DB.SolverType;
	unsigned int Neq         = DB.Neq,
	             Collocated  = DB.Collocated;

	// Standard datatypes
	unsigned int iMax, jMax,
	             VfIn, VfOut, fIn, fOut, NvnSIn, NvnSOut, Boundary;
	double       maxRHS, maxRHSV, *VRHSIn_ptr, *VRHSOut_ptr, *FRHSIn_ptr, *FRHSOut_ptr, *wdetJV_vI, *wdetJV_vI_ptr,
	             *RHS_Final;

	struct S_VOLUME  *VIn, *VOut, *VOLUME;
	struct S_FACET   *FACET;

	for (FACET = DB.FACET; FACET != NULL; FACET = FACET->next) {
		VIn    = FACET->VIn;
		VfIn   = FACET->VfIn;
		fIn    = VfIn/NFREFMAX;
		NvnSIn = VIn->NvnS;

		VOut    = FACET->VOut;
		VfOut   = FACET->VfOut;
		fOut    = VfOut/NFREFMAX;
		NvnSOut = VOut->NvnS;

		VRHSIn_ptr  = VIn->RHS;
		VRHSOut_ptr = VOut->RHS;

		FRHSIn_ptr  = FACET->RHSIn;
		FRHSOut_ptr = FACET->RHSOut;

/*
printf("%d\n",FACET->indexg);
array_print_d(NvnSIn,Neq,VRHSIn_ptr,'C');
array_print_d(NvnSOut,Neq,VRHSOut_ptr,'C');
//array_print_d(NvnSIn,Neq,FRHSIn_ptr,'C');
//array_print_d(NvnSOut,Neq,FRHSOut_ptr,'C');
*/

		Boundary = !((VIn->indexg != VOut->indexg) || (VIn->indexg == VOut->indexg && fIn != fOut));
		for (iMax = Neq; iMax--; ) {
			for (jMax = NvnSIn; jMax--; )
				*VRHSIn_ptr++ += *FRHSIn_ptr++;
			if (!Boundary) {
				for (jMax = NvnSOut; jMax--; )
					*VRHSOut_ptr++ += *FRHSOut_ptr++;
			}
		}
/*
printf("%d\n",FACET->indexg);
VRHSIn_ptr  = VIn->RHS;
VRHSOut_ptr = VOut->RHS;
array_print_d(NvnSIn,Neq,VRHSIn_ptr,'C');
array_print_d(NvnSOut,Neq,VRHSOut_ptr,'C');
*/

		free(FACET->RHSIn);
		free(FACET->RHSOut);
	}

	// Add MInv contribution to RHS for explicit runs
	maxRHS = 0.0;
	if (strstr(SolverType,"Explicit") != NULL) {
		for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
			// Add (remaining) MInv contribution to RHS
			NvnSIn = VOLUME->NvnS;

			if (Collocated) {
				wdetJV_vI = VOLUME->wdetJV_vI;

				VRHSIn_ptr = VOLUME->RHS;
				for (iMax = Neq; iMax--; ) {
					wdetJV_vI_ptr = wdetJV_vI;
					for (jMax = NvnSIn; jMax--; )
						*VRHSIn_ptr++ /= *wdetJV_vI_ptr++;
				}
			} else {
				RHS_Final = malloc(NvnSIn*Neq * sizeof *RHS_Final);

				mm_CTN_d(NvnSIn,Neq,NvnSIn,VOLUME->MInv,VOLUME->RHS,RHS_Final);
				free(VOLUME->RHS);
				VOLUME->RHS = RHS_Final;
/*
printf("%d\n",VOLUME->indexg);
array_print_d(NvnSIn,NvnSIn,VOLUME->MInv,'R');
array_print_d(NvnSIn,Neq,VOLUME->RHS,'C');
array_print_d(NvnSIn,Neq,RHS_Final,'C');
*/
			}

			// Compute maxRHS for convergence monitoring
			maxRHSV = array_norm_d(NvnSIn,VOLUME->RHS,"Inf");
			if (maxRHSV > maxRHS)
				maxRHS = maxRHSV;
		}
	}
	return maxRHS;
}
