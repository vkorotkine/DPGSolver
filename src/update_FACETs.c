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
 *		Update FACET information/operators in ELEMENTs which have undergone hp refinement.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void update_FACET_hp(void)
{
	/*
	 *	Comments:
	 *		Parameters to (potentially) modify:
	 *			P, typeInt, BC, indexg,
	 *			VIn/VOut, VfIn/VfOut, IndOrdInOut/IndOrdOutIn,
	 *			n_fS, XYZ_fS, detJF_fS
	 */
	// Initialize DB Parameters
	unsigned int Adapt = DB.Adapt;

	// Standard datatypes
	unsigned int dummy_ui;

	struct S_VOLUME    *VIn, *VOut, *dummyPtr_V;
	struct S_FACET     *FACET;

	switch (Adapt) {
	default: // ADAPT_H, ADAPT_HP
		break;
	case ADAPT_P:
		/*	No modifications required for:
		 *		indexg, typeInt, BC
		 */
		for (FACET = DB.FACET; FACET != NULL; FACET = FACET->next) {
			VIn  = FACET->VIn;
			VOut = FACET->VOut;

			if (VIn->update || VOut->update) {
				FACET->P = max(VIn->P,VOut->P);

				// Ensure that VIn is the highest order VOLUME
				if (VIn->P < VOut->P) {
					FACET->VIn  = VOut;
					FACET->VOut = VIn;

					dummy_ui     = FACET->VfIn;
					FACET->VfIn  = FACET->VfOut;
					FACET->VfOut = dummy_ui;

					dummy_ui           = FACET->IndOrdInOut;
					FACET->IndOrdInOut = FACET->IndOrdOutIn;
					FACET->IndOrdOutIn = dummy_ui;
				}

				// Recompute XYZ_fS, n_fS, and detJF_fS
				free(FACET->XYZ_fS);
				setup_FACET_XYZ(FACET);

				free(FACET->n_fS);
				free(FACET->detJF_fS);
				setup_normals(FACET);
			}
		}
		break;
	}
}
