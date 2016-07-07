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

double finalize_RHS(void)
{
	// Initialize DB Parameters
	char         *SolverType = DB.SolverType;
	unsigned int Neq         = DB.Neq,
	             Collocated  = DB.Collocated;

	// Standard datatypes
	unsigned int iMax, jMax,
	             VfIn, VfOut, fIn, fOut, NvnSIn, NvnSOut, Boundary;
	double       maxRHS, maxRHSV, *VRHSIn_ptr, *VRHSOut_ptr, *FRHSIn_ptr, *FRHSOut_ptr, *detJV_vI_ptr,
	             *RHS_Final;

	struct S_VOLUME  *VIn, *VOut, *VOLUME;
	struct S_FACET   *FACET;

	// silence
	NvnSIn = 0;

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

/*
printf("%d\n",FACET->indexg);
array_print_d(NvnSIn,Neq,VIn->RHS,'C');
array_print_d(NvnSOut,Neq,VOut->RHS,'C');
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
//if (VIn->indexg < 2) {
printf("%d\n",VIn->indexg);
array_print_d(NvnSIn,Neq,VIn->RHS,'C');
//}
//if (VOut->indexg < 2) {
printf("%d\n",VOut->indexg);
array_print_d(NvnSOut,Neq,VOut->RHS,'C');
//}
*/


		free(FACET->RHSIn);
		free(FACET->RHSOut);
	}
//exit(1);

	// Add MInv contribution to RHS for explicit runs
	maxRHS = 0.0;
	if (strstr(SolverType,"Explicit") != NULL) {
		for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
			// Compute maxRHS for convergence monitoring
			NvnSIn = VOLUME->NvnS;
			maxRHSV = array_norm_d(NvnSIn,VOLUME->RHS,"Inf");
			if (maxRHSV > maxRHS)
				maxRHS = maxRHSV;

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

/*
printf("%d\n",VOLUME->indexg);
array_print_d(NvnSIn,NvnSIn,VOLUME->MInv,'R');
array_print_d(NvnSIn,Neq,VOLUME->RHS,'C');
*/

				mm_CTN_d(NvnSIn,Neq,NvnSIn,VOLUME->MInv,VOLUME->RHS,RHS_Final);
				free(VOLUME->RHS);
				VOLUME->RHS = RHS_Final;

/*
//printf("%d\n",VOLUME->indexg);
array_print_d(NvnSIn,Neq,RHS_Final,'C');
//exit(1);
*/

			}

		}
	}
//exit(1);
	return maxRHS;
}
