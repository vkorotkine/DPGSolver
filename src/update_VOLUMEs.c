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
 *		May be better not to allow P adaptivity to go down to P = 0, but have the limit at P = 1. Try/THINK (ToBeDeleted).
 *		Sum factorization not currently used here. Profile and re-evaluate. (ToBeDeleted)
 *
 *	Notation:
 *
 *	References:
 */

struct S_OPERATORS {
	unsigned int NvnGs, NvnGc, NvnS, NvnSP, NvnI;
	double       *I_vGs_vGc, *Ihat_vS_vS, *w_vI, *ChiS_vI;
};

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const unsigned int IndClass)
{
	unsigned int P, PNew, type, curved;
	struct S_ELEMENT *ELEMENT;

	P      = VOLUME->P;
	PNew   = VOLUME->PNew;
	type   = VOLUME->type;
	curved = VOLUME->curved;

	ELEMENT = get_ELEMENT_type(type);

	OPS->NvnS  = ELEMENT->NvnS[P];
	OPS->NvnSP = ELEMENT->NvnS[PNew];
	OPS->Ihat_vS_vS = ELEMENT->Ihat_vS_vS[P][PNew];
	if (!curved) {
		OPS->NvnI = ELEMENT->NvnIs[P];

		OPS->w_vI    = ELEMENT->w_vIs[P];
		OPS->ChiS_vI = ELEMENT->ChiS_vIs[P][P][0];
	} else {
		OPS->NvnGs = ELEMENT->NvnGs[0];
		OPS->NvnGc = ELEMENT->NvnGc[PNew];
		OPS->NvnI  = ELEMENT->NvnIc[P];

		OPS->I_vGs_vGc = ELEMENT->I_vGs_vGc[0][PNew][0];

		OPS->w_vI    = ELEMENT->w_vIc[P];
		OPS->ChiS_vI = ELEMENT->ChiS_vIc[P][P][0];
	}
}

void update_VOLUME_hp(void)
{
	// Initialize DB Parameters
	unsigned int d    = DB.d,
	             PMax = DB.PMax,
	             Nvar = DB.Nvar;

	char         *MeshType = DB.MeshType;

	// Standard datatypes
	unsigned int P, PNew, adapt_type;
	double       NvnGs, NvnGc, NvnS, NvnSP, NCols,
	             *I_vGs_vGc, *XYZ_vC, *XYZ_S,
	             *Ihat_vS_vS, *What, *RES, *WhatP, *RESP;

	struct S_OPERATORS *OPS;
	struct S_VOLUME    *VOLUME;

	OPS = malloc(sizeof *OPS); // free

	for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
		if (VOLUME->Vadapt) {
			VOLUME->update = 1;
			P = VOLUME->P;
			adapt_type = VOLUME->adapt_type;

			switch(adapt_type) {
			case PREFINE:
				if (P+1 <= PMax)
					PNew = P+1;
				else
					printf("Error: Should not be entering PREFINE in update_VOLUME_hp for P = %d.\n",P), exit(1);
				VOLUME->PNew = PNew;
				break;
			case PCOARSE:
				if (P >= 1)
					PNew = P-1;
				else
					printf("Error: Should not be entering PREFINE in update_VOLUME_hp for P = %d.\n",P), exit(1);
				VOLUME->PNew = PNew;
				break;
			case HREFINE:
				// Vh as flag for type of VOLUME h refinement
				// Also need another flag to make sure that no h refinement is performed on elements of the coarsest
				// mesh.
				VOLUME->PNew = P;
				break;
			case HCOARSE:
				VOLUME->PNew = P;
				break;
			default:
				printf("Error: Unsupported adapt_type = %d in update_VOLUME_hp.\n",adapt_type), exit(1);
				break;
			}

			init_ops(OPS,VOLUME,0);

			switch(adapt_type) {
			default: // PREFINE or PCOARSE
				// Update geometry
// need normals and detJF_fI
				if (VOLUME->curved) {
					NvnGs      = OPS->NvnGs;
					NvnGc      = OPS->NvnGc;
					I_vGs_vGc  = OPS->I_vGs_vGc;

					NCols = d;

					XYZ_vC = VOLUME->XYZ_vC;
					XYZ_S  = malloc(NvnGc*NCols * sizeof *XYZ_S); // keep
					mm_CTN_d(NvnGc,NCols,NvnGs,I_vGs_vGc,XYZ_vC,XYZ_S);
					free(VOLUME->XYZ_S);
					VOLUME->XYZ_S = XYZ_S;
				}

				if (strstr(MeshType,"ToBeCurved") != NULL)
					setup_ToBeCurved(VOLUME);
				else
					printf("Error: Add in support for MeshType != ToBeCurved\n"), exit(1);

				setup_geom_factors(VOLUME); // ToBeDeleted: Don't forget to free VOLUME->C_vC

				// Project What and RES
				NvnS       = OPS->NvnS;
				NvnSP      = OPS->NvnSP;
				Ihat_vS_vS = OPS->Ihat_vS_vS;

				VOLUME->NvnS = NvnSP;

				What  = VOLUME->What;
				RES   = VOLUME->RES;

				WhatP = malloc(NvnSP*Nvar * sizeof *WhatP); // keep
				RESP  = malloc(NvnSP*Nvar * sizeof *RESP);  // keep

				mm_CTN_d(NvnSP,Nvar,NvnS,Ihat_vS_vS,What,WhatP);
				mm_CTN_d(NvnSP,Nvar,NvnS,Ihat_vS_vS,RES,RESP);

				free(What);
				free(RES);

				VOLUME->What = WhatP;
				VOLUME->RES  = RESP;

				VOLUME->P = PNew;
				break;
			case HREFINE:
				// Interpolate to finer space
				break;
			case HCOARSE:
				// Galerkin projection to coarser space
				break;
			}
		}
	}

	free(OPS);
}

void update_Vgrp(void)
{
	// Initialize DB Parameters
	unsigned int NP     = DB.NP,
	             NTVgrp = DB.NTVgrp,
	             *NVgrp = DB.NVgrp;

	// Standard datatypes
	unsigned int i, P, Eclass, curved, IndVgrp;

	struct S_VOLUME *VOLUME, **VOLUME_prev;

	if (DB.update) {
		DB.update = 0;

		VOLUME_prev = malloc(NTVgrp * sizeof *VOLUME_prev); // free

		for (i = 0; i < NTVgrp; i++)
			NVgrp[i] = 0;

		for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
			Eclass = VOLUME->Eclass;
			P      = VOLUME->P;
			curved = VOLUME->curved;

			IndVgrp = Eclass*NP*2 + P*2 + curved;

			if (!NVgrp[IndVgrp])
				DB.Vgrp[IndVgrp] = VOLUME;
			else
				VOLUME_prev[IndVgrp]->grpnext = VOLUME;

			NVgrp[IndVgrp]++;
			VOLUME_prev[IndVgrp] = VOLUME;
		}
		free(VOLUME_prev);
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
				if (!Collocated) {
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
					free(VOLUME->wdetJV_vI);
					free(VOLUME->MInv);

					VOLUME->wdetJV_vI = wdetJV_vI;
					VOLUME->MInv      = MInv;
				}
			} else {
				; // Updates for implicit.
			}
		}
	}

	free(OPS);
}
