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
 *		Update VOLUME operators in ELEMENTs which have undergone hp refinement.
 *
 *	Comments:
 *		Think whether it would be possible to add sum factorization capability here for MInv computation (ToBeDeleted).
 *
 *	Notation:
 *
 *	References:
 */

struct S_OPERATORS {
	unsigned int NvnS, NvnI;
	double       *w_vI, *ChiS_vI;
};

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const unsigned int IndClass)
{
	unsigned int P, type, curved;
	struct S_ELEMENT *ELEMENT, *ELEMENT_OPS;

	// silence
	ELEMENT_OPS = NULL;

	P      = VOLUME->P;
	type   = VOLUME->type;
	curved = VOLUME->curved;

	ELEMENT = get_ELEMENT_type(type);
	if (1 || type == TRI || type == TET || type == PYR)
		ELEMENT_OPS = ELEMENT;
	else if (type == LINE || type == QUAD || type == HEX || type == WEDGE)
		ELEMENT_OPS = ELEMENT->ELEMENTclass[IndClass];

	OPS->NvnS = ELEMENT_OPS->NvnS[P];
	if (!curved) {
		OPS->NvnI = ELEMENT_OPS->NvnIs[P];

		OPS->w_vI    = ELEMENT_OPS->w_vIs[P];
		OPS->ChiS_vI = ELEMENT_OPS->ChiS_vIs[P][P][0];
	} else {
		OPS->NvnI = ELEMENT_OPS->NvnIc[P];

		OPS->w_vI    = ELEMENT_OPS->w_vIc[P];
		OPS->ChiS_vI = ELEMENT_OPS->ChiS_vIc[P][P][0];
	}
}

void update_VOLUME_Ops(void)
{
	// Initialize DB Parameters
	char         *SolverType = DB.SolverType;
	unsigned int Collocated  = DB.Collocated;

	// Standard datatypes
	unsigned int iMax, jMax,
	             NvnS, NvnI;
	double       *detJV_vI, *detJV_vI_ptr, *w_vI, *w_vI_ptr, *wdetJV_vI, *wdetJV_vI_ptr, *ChiS_vI, *ChiS_vI_ptr,
	             *wdetJVChiS_vI, *wdetJVChiS_vI_ptr, *IS, *M, *MInv;

	struct S_OPERATORS *OPS;
	struct S_VOLUME    *VOLUME;

	OPS = malloc(sizeof *OPS); // free

	for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
		if (VOLUME->update) {
			VOLUME->update = 0;
			if (strstr(SolverType,"Explicit") != NULL) {
				init_ops(OPS,VOLUME,0);

				NvnS = OPS->NvnS;
				NvnI = OPS->NvnI;
				w_vI    = OPS->w_vI;
				ChiS_vI = OPS->ChiS_vI;

				detJV_vI = VOLUME->detJV_vI;

				// Compute required portion of MInv
				wdetJV_vI = malloc(NvnI * sizeof *wdetJV_vI); // keep
				MInv      = NULL;

				w_vI_ptr      = w_vI;
				detJV_vI_ptr  = detJV_vI;
				wdetJV_vI_ptr = wdetJV_vI;
				for (iMax = NvnI; iMax--; )
					*wdetJV_vI_ptr++ = (*w_vI_ptr++)*(*detJV_vI_ptr++);

				if (!Collocated) {
					wdetJVChiS_vI = malloc(NvnI*NvnS * sizeof *wdetJVChiS_vI); // free

					ChiS_vI_ptr       = ChiS_vI;
					wdetJVChiS_vI_ptr = wdetJVChiS_vI;
					for (iMax = NvnI*NvnS; iMax--; )
						*wdetJVChiS_vI_ptr++ = *ChiS_vI_ptr++;

					wdetJV_vI_ptr     = wdetJV_vI;
					wdetJVChiS_vI_ptr = wdetJVChiS_vI;
					for (iMax = NvnI; iMax--; ) {
						for (jMax = NvnS; jMax--; )
							*wdetJVChiS_vI_ptr++ *= *wdetJV_vI_ptr;
						wdetJV_vI_ptr++;
					}

					M    = mm_Alloc_d(CBRM,CBT,CBNT,NvnS,NvnS,NvnI,1.0,ChiS_vI,wdetJVChiS_vI); // free
					IS   = identity_d(NvnS); // free
					MInv = inverse_d(NvnS,NvnS,M,IS); // keep

					free(wdetJVChiS_vI);
					free(M);
					free(IS);
				}
				free(VOLUME->wdetJV_vI);
				free(VOLUME->MInv);

				VOLUME->wdetJV_vI = wdetJV_vI;
				VOLUME->MInv      = MInv;

			} else {
				; // Updates for implicit.
			}
		}
	}

	free(OPS);
}
