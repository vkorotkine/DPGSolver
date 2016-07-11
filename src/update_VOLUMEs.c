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
 *		Update VOLUME related information/operators in ELEMENTs which have undergone hp refinement.
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
	double       *I_vGs_vGc, **I_vGs_vGs, **Ihat_vS_vS, **Ghat_vS_vS, *w_vI, *ChiS_vI;
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

	OPS->NvnGs = ELEMENT->NvnGs[1];
	OPS->NvnGc = ELEMENT->NvnGc[PNew];
	OPS->NvnS  = ELEMENT->NvnS[P];
	OPS->NvnSP = ELEMENT->NvnS[PNew];
	OPS->I_vGs_vGs  = ELEMENT->I_vGs_vGs[1][1];
	OPS->I_vGs_vGc  = ELEMENT->I_vGs_vGc[1][PNew][0];
	OPS->Ihat_vS_vS = ELEMENT->Ihat_vS_vS[P][PNew]; // ToBeDeleted: Remove all instances of Ihat_vS_vS from the code if
	                                                //              not used here.
	OPS->Ghat_vS_vS = ELEMENT->Ghat_vS_vS[P][PNew];
	if (!curved) {
		OPS->NvnI = ELEMENT->NvnIs[P];

		OPS->w_vI    = ELEMENT->w_vIs[P];
		OPS->ChiS_vI = ELEMENT->ChiS_vIs[P][P][0];
	} else {
		OPS->NvnI  = ELEMENT->NvnIc[P];

		OPS->w_vI    = ELEMENT->w_vIc[P];
		OPS->ChiS_vI = ELEMENT->ChiS_vIc[P][P][0];
	}
}

static unsigned int get_VOLUMEh_type(const unsigned int VType, const unsigned int vh)
{
	switch (VType) {
	case TET:
		if (vh < 4)
			return TET;
		else
			return PYR;
		break;
	case PYR:
		if (vh < 4 || vh > 7)
			return PYR;
		else
			return TET;
		break;
	default:
		printf("Error: Unsupported VType in get_VOLUMEh_type.\n"), exit(1);
		break;
	}
}

void update_VOLUME_hp(void)
{
	// Initialize DB Parameters
	unsigned int d         = DB.d,
	             NV        = DB.NV,
	             AC        = DB.AC,
	             PMax      = DB.PMax,
				 LevelsMax = DB.LevelsMax,
	             Nvar      = DB.Nvar;

	char         *MeshType = DB.MeshType;

	// Standard datatypes
	unsigned int P, PNew, f, sf, level, adapt_type, vh, vhMin, vhMax, sfMax, href_type, VType, Nf,
	             fInd,
	             NvnGs, NvnGc, NvnG, NvnS, NvnSP, NCols, *NsubF;
	double       *I_vGs_vGc, *XYZ_vC, *XYZ_S,
	             **Ihat_vS_vS, **I_vGs_vGs, **Ghat_vS_vS, *What, *RES, *WhatP, *RESP;

	struct S_OPERATORS *OPS;
	struct S_ELEMENT   *ELEMENT;
	struct S_VOLUME    *VOLUME, *VOLUMEh, **VOLUME_ptr;
	struct S_FACET     *FACET;

	// silence
	I_vGs_vGs = NULL;
	VOLUMEh   = NULL;

	OPS = malloc(sizeof *OPS); // free

	VOLUME_ptr = &(DB.VOLUME);
	for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
		if (VOLUME->Vadapt) {
			P     = VOLUME->P;
			level = VOLUME->level;
			adapt_type = VOLUME->adapt_type;

			switch(adapt_type) {
			case PREFINE:
				if (P < PMax)
					PNew = P+1;
				else
					printf("Error: Should not be entering PREFINE in update_VOLUME_hp for P = %d.\n",P), exit(1);
				VOLUME->PNew = PNew;
				break;
			case PCOARSE:
				if (P >= 1)
					PNew = P-1;
				else
					printf("Error: Should not be entering PCOARSE in update_VOLUME_hp for P = %d.\n",P), exit(1);
				VOLUME->PNew = PNew;
				break;
			case HREFINE:
				if (level >= LevelsMax)
					printf("Error: Should not be entering HREFINE in update_VOLUME_hp for level = %d.\n",level), exit(1);
				VOLUME->PNew = P;
				break;
			case HCOARSE:
				if (level == 0)
					printf("Error: Should not be entering HCOARSE in update_VOLUME_hp for level = %d.\n",level), exit(1);
				VOLUME->PNew = P;
				break;
			default:
				printf("Error: Unsupported adapt_type = %d in update_VOLUME_hp.\n",adapt_type), exit(1);
				break;
			}

			init_ops(OPS,VOLUME,0);

			VOLUME->update = 1;
			switch(adapt_type) {
			default: // PREFINE or PCOARSE
				VOLUME->P = PNew;

				// Update geometry
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
					VOLUME->NvnG  = NvnGc;
				}

				free(VOLUME->XYZ);
				if (strstr(MeshType,"ToBeCurved") != NULL)
					setup_ToBeCurved(VOLUME);
				else
					printf("Error: Add in support for MeshType != ToBeCurved\n"), exit(1);

				free(VOLUME->detJV_vI);
				free(VOLUME->C_vI);
				setup_geom_factors(VOLUME);

				// Project What and RES
				NvnS       = OPS->NvnS;
				NvnSP      = OPS->NvnSP;
				Ihat_vS_vS = OPS->Ihat_vS_vS;
				Ghat_vS_vS = OPS->Ghat_vS_vS;

				VOLUME->NvnS = NvnSP;

				What  = VOLUME->What;
				RES   = VOLUME->RES;

				WhatP = malloc(NvnSP*Nvar * sizeof *WhatP); // keep
				RESP  = malloc(NvnSP*Nvar * sizeof *RESP);  // keep

//				mm_CTN_d(NvnSP,Nvar,NvnS,Ihat_vS_vS[0],What,WhatP);
//				mm_CTN_d(NvnSP,Nvar,NvnS,Ihat_vS_vS[0],RES,RESP);
				mm_CTN_d(NvnSP,Nvar,NvnS,Ghat_vS_vS[0],What,WhatP);
				mm_CTN_d(NvnSP,Nvar,NvnS,Ghat_vS_vS[0],RES,RESP);

				free(What);
				free(RES);

				VOLUME->What = WhatP;
				VOLUME->RES  = RESP;

				VOLUME_ptr = &VOLUME->next;
				break;
			case HREFINE:
				VType = VOLUME->type;
				NsubF = VOLUME->NsubF;


				NvnGs     = OPS->NvnGs;
				NvnGc     = OPS->NvnGc;
				I_vGs_vGs = OPS->I_vGs_vGs;
				I_vGs_vGc = OPS->I_vGs_vGc;

				NCols = d;


				ELEMENT = get_ELEMENT_type(VType);
				Nf = ELEMENT->Nf;

				href_type = VOLUME->hrefine_type;

				get_vh_range(VType,href_type,&vhMin,&vhMax);
				for (vh = vhMin; vh <= vhMax; vh++) {
					if (vh == vhMin) {
						VOLUMEh = New_VOLUME();
						VOLUME->child0 = VOLUMEh;
					} else {
						VOLUMEh->next = New_VOLUME();
						VOLUMEh = VOLUMEh->next;
					}
					VOLUMEh->update = 1;
					VOLUMEh->parent = VOLUME;
					if (vh == vhMin)
						*VOLUME_ptr = VOLUMEh;
					else if (vh == vhMax)
						VOLUMEh->next = VOLUME->next;
					VOLUMEh->indexg = NV++;

					VOLUMEh->P = VOLUME->P;
					VOLUMEh->level = (VOLUME->level)+1;
					switch (VType) {
					default: // LINE, TRI, QUAD, HEX, WEDGE
						VOLUMEh->type = VType;
						break;
					case TET:
					case PYR:
						VOLUMEh->type = get_VOLUMEh_type(VType,vh);
						break;
					}
					VOLUMEh->Eclass = get_Eclass(VOLUMEh->type);

					if (AC) {
						VOLUMEh->curved = 1;
					} else if (VOLUME->curved) {
						printf("Error: Add support for h-refinement VOLUMEh->curved.\n"), exit(1);
						// Use VToBC and knowledge of whether the new VOLUME shares the BC.
					} else {
						VOLUMEh->curved = 0;
					}

					// Update geometry
// When updating XYZ_vC, ensure that corners on curved boundaries are placed on the boundary.
					VOLUMEh->XYZ_vC = malloc(NvnGs*d * sizeof *XYZ_vC); // keep
					XYZ_vC = VOLUMEh->XYZ_vC;

					mm_CTN_d(NvnGs,NCols,NvnGs,I_vGs_vGs[vh],VOLUME->XYZ_vC,VOLUMEh->XYZ_vC);
					if (!VOLUMEh->curved) {
						double *XYZ;

						VOLUMEh->NvnG = NvnGs;

						VOLUMEh->XYZ_S = malloc(NvnGs*NCols * sizeof *XYZ_S); // keep
						VOLUMEh->XYZ   = malloc(NvnGs*NCols * sizeof *XYZ);   // keep
						XYZ_S = VOLUMEh->XYZ_S;
						XYZ   = VOLUMEh->XYZ;
						for (unsigned int i = 0, iMax = NCols*NvnGs; i < iMax; i++) {
							XYZ_S[i] = XYZ_vC[i];
							XYZ[i]   = XYZ_S[i];
						}
					} else {
						VOLUMEh->NvnG = NvnGc;

						VOLUMEh->XYZ_S = malloc(NvnGc*NCols * sizeof *XYZ_S); // keep
						mm_CTN_d(NvnGc,NCols,NvnGs,I_vGs_vGc,XYZ_vC,VOLUMEh->XYZ_S);

						if (strstr(MeshType,"ToBeCurved") != NULL) {
							setup_ToBeCurved(VOLUMEh);
						} else {
							printf("Add in support for MeshType != ToBeCurved");
							exit(1);
						}
					}
					setup_geom_factors(VOLUMEh);
				}

				// Fix VOLUME linked list and Vgrp linked list
				// Also update indexg and indexl

				// When updating connectivity, start with groups of elements of the lowest level and move to those with
				// the highest level to avoid having more than 1 level of non-conformity.


				// Project What and RES
				// free VOLUME->What and VOLUME->RES after projection.




				break;
			case HCOARSE:
			// Don't forget VOLUME_ptr
				// Galerkin projection to coarser space
				break;
			}
		} else {
			VOLUME_ptr = &VOLUME->next;
		}
	}
	free(OPS);

	DB.NV = NV;
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
//printf("update_VOLUME_Ops: %d\n",VOLUME->indexg);
			if (strstr(SolverType,"Explicit") != NULL) {
				if (!Collocated) {
					init_ops(OPS,VOLUME,0);

					NvnS = OPS->NvnS;
					NvnI = OPS->NvnI;
					w_vI    = OPS->w_vI;
					ChiS_vI = OPS->ChiS_vI;

					detJV_vI = VOLUME->detJV_vI;

					// Compute required portion of MInv
					wdetJV_vI = malloc(NvnI * sizeof *wdetJV_vI); // free
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
//array_print_d(NvnS,NvnS,M,'R');
//array_print_d(NvnI,NvnS,wdetJVChiS_vI,'R');
//exit(1);

					free(wdetJVChiS_vI);
					free(M);
					free(IS);
					free(wdetJV_vI);

					free(VOLUME->MInv);
					VOLUME->MInv = MInv;
				}
			} else {
				; // Updates for implicit.
			}
		}
	}

	free(OPS);
}

void update_VOLUME_finalize(void)
{
	struct S_VOLUME *VOLUME;

	for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
		if (VOLUME->Vadapt) {
			VOLUME->Vadapt = 0;
		}
	}
}
