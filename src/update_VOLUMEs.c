// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "update_VOLUMEs.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <limits.h>

#include "mkl.h"

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "S_FACET.h"

#include "adaptation.h"
#include "element_functions.h"
#include "matrix_functions.h"
#include "setup_ToBeCurved.h"
#include "setup_geom_factors.h"
#include "memory_constructors.h"

#include "array_print.h"

/*
 *	Purpose:
 *		Update VOLUME related information/operators in ELEMENTs which have undergone hp refinement.
 *
 *	Comments:
 *		May be better not to allow P adaptivity to go down to P = 0, but have the limit at P = 1. Try/THINK (ToBeDeleted).
 *		Sum factorization not currently used here. Profile and re-evaluate. (ToBeDeleted)
 *		If the HCOARSE projection is slow, consider allowing the BLAS call with beta = 1.0 (ToBeDeleted).
 *
 *	Notation:
 *
 *	References:
 */

struct S_OPERATORS {
	unsigned int NvnGs, NvnGc, NvnS, *Nvve, NvnSP, NvnI;
	double       *I_vGs_vGc, **I_vGs_vGs, **Ihat_vS_vS, **L2hat_vS_vS, *w_vI, *ChiS_vI;
};

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const unsigned int IndEType_h)
{
	unsigned int P, PNew, type, curved;
	struct S_ELEMENT *ELEMENT;

	P      = VOLUME->P;
	PNew   = VOLUME->PNew;
	type   = VOLUME->type;
	curved = VOLUME->curved;

	ELEMENT = get_ELEMENT_type(type);
	if (IndEType_h)
		ELEMENT = get_ELEMENT_type(ELEMENT->type_h[IndEType_h]);

	OPS->NvnGs = ELEMENT->NvnGs[1];
	OPS->NvnGc = ELEMENT->NvnGc[PNew];
	OPS->NvnS  = ELEMENT->NvnS[P];
	OPS->Nvve  = ELEMENT->Nvve;
	OPS->NvnSP = ELEMENT->NvnS[PNew];
	OPS->I_vGs_vGs  = ELEMENT->I_vGs_vGs[1][1];
	OPS->I_vGs_vGc  = ELEMENT->I_vGs_vGc[1][PNew][0];
	OPS->Ihat_vS_vS = ELEMENT->Ihat_vS_vS[P][PNew]; // ToBeDeleted: Remove all instances of Ihat_vS_vS from the code if
	                                                //              not used here.
	OPS->L2hat_vS_vS = ELEMENT->L2hat_vS_vS[P][PNew];
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

void update_VOLUME_hp(void)
{
	// Initialize DB Parameters
	unsigned int d         = DB.d,
	             NV        = DB.NV,
	             AC        = DB.AC,
	             PMax      = DB.PMax,
				 LevelsMax = DB.LevelsMax,
	             Nvar      = DB.Nvar;

	char         *MeshType = DB.MeshType,
	             *TestCase = DB.TestCase;

	// Standard datatypes
	unsigned int i, iMax, P, PNew, f, level, adapt_type, vh, vhMin, vhMax, VType, Nf,
	             IndEhref, NvnGs[2], NvnGc[2], NvnS[2], NvnSP, NCols, update, maxP;
	double       *I_vGs_vGc[2], *XYZ_vC, *XYZ_S,
	             **Ihat_vS_vS, **I_vGs_vGs, **L2hat_vS_vS, *What, *RES, *WhatP, *WhatH, *RESP, *RESH, *dummyPtr_d,
	             *uhat, *uhatP, *uhatH;

	struct S_OPERATORS *OPS;
	struct S_ELEMENT   *ELEMENT;
	struct S_VOLUME    *VOLUME, *VOLUMEc, *VOLUMEp;

	// silence
	I_vGs_vGs = NULL;
	VOLUMEc   = NULL;

	OPS = malloc(sizeof *OPS); // free

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		if (VOLUME->Vadapt) {
			P     = VOLUME->P;
			level = VOLUME->level;
			adapt_type = VOLUME->adapt_type;

			switch(adapt_type) {
			case PREFINE:
				if (P < PMax)
					PNew = P+1;
				else
					printf("Error: Should not be entering PREFINE for P = %d.\n",P), EXIT_MSG;
				VOLUME->PNew = PNew;
				break;
			case PCOARSE:
				if (P >= 1)
					PNew = P-1;
				else
					printf("Error: Should not be entering PCOARSE for P = %d.\n",P), EXIT_MSG;
				VOLUME->PNew = PNew;
				break;
			case HREFINE:
				if (level == LevelsMax)
					printf("Error: Should not be entering HREFINE for level = %d.\n",level), EXIT_MSG;
				VOLUME->PNew = P;
				break;
			case HCOARSE:
				if (level == 0)
					printf("Error: Should not be entering HCOARSE for level = %d.\n",level), EXIT_MSG;
				VOLUME->PNew = P;
				break;
			default:
				printf("Error: Unsupported adapt_type = %d.\n",adapt_type), EXIT_MSG;
				break;
			}

			if (adapt_type == HREFINE) {
				if (VOLUME->type == PYR || VOLUME->type == TET)
					init_ops(OPS,VOLUME,1);
				else
					init_ops(OPS,VOLUME,0);

				NvnGs[1]     = OPS->NvnGs;
				NvnGc[1]     = OPS->NvnGc;
				NvnS[1]      = OPS->NvnS;
				I_vGs_vGc[1] = OPS->I_vGs_vGc;
			} else if (adapt_type == HCOARSE) {
				if ((VOLUME->type == TET && VOLUME->parent->type == TET) ||
				    (VOLUME->type == PYR && VOLUME->parent->type == PYR))
					init_ops(OPS,VOLUME,1);
				else
					init_ops(OPS,VOLUME,0);

				NvnS[1]      = OPS->NvnS;
			}

			init_ops(OPS,VOLUME,0);

			VOLUME->update = 1;
			switch (adapt_type) {
			default: // PREFINE or PCOARSE
				VOLUME->P = PNew;

				// Update geometry
				if (VOLUME->curved) {
					NvnGs[0]     = OPS->NvnGs;
					NvnGc[0]     = OPS->NvnGc;
					I_vGs_vGc[0] = OPS->I_vGs_vGc;

					NCols = d;

					XYZ_vC = VOLUME->XYZ_vC;
					XYZ_S  = malloc(NvnGc[0]*NCols * sizeof *XYZ_S); // keep
					mm_CTN_d(NvnGc[0],NCols,NvnGs[0],I_vGs_vGc[0],XYZ_vC,XYZ_S);

					free(VOLUME->XYZ_S);
					VOLUME->XYZ_S = XYZ_S;
					VOLUME->NvnG  = NvnGc[0];

					free(VOLUME->XYZ);
					if (strstr(MeshType,"ToBeCurved"))
						setup_ToBeCurved(VOLUME);
					else if (strstr(MeshType,"Curved"))
						printf("Add in support for MeshType == Curved.\n"), EXIT_MSG;
				}

				free(VOLUME->detJV_vI);
				free(VOLUME->C_vI);
				setup_geom_factors(VOLUME);

				// Project solution coefficients (and RES if applicable)
				NvnS[0]    = OPS->NvnS;
				NvnSP      = OPS->NvnSP;

				VOLUME->NvnS = NvnSP;

				if (strstr(TestCase,"Poisson")) {
					uhat = VOLUME->uhat;

					uhatP = malloc(NvnSP*Nvar * sizeof *uhatP); // keep

					if (adapt_type == PREFINE) {
						Ihat_vS_vS = OPS->Ihat_vS_vS;
						mm_CTN_d(NvnSP,Nvar,NvnS[0],Ihat_vS_vS[0],uhat,uhatP);
					} else {
						L2hat_vS_vS = OPS->L2hat_vS_vS;
						mm_CTN_d(NvnSP,Nvar,NvnS[0],L2hat_vS_vS[0],uhat,uhatP);
					}
					free(uhat);
					VOLUME->uhat = uhatP;
				} else {
					What = VOLUME->What;
					RES  = VOLUME->RES;

					WhatP = malloc(NvnSP*Nvar * sizeof *WhatP); // keep
					RESP  = malloc(NvnSP*Nvar * sizeof *RESP);  // keep

					if (adapt_type == PREFINE) {
						Ihat_vS_vS = OPS->Ihat_vS_vS;
						mm_CTN_d(NvnSP,Nvar,NvnS[0],Ihat_vS_vS[0],What,WhatP);
						mm_CTN_d(NvnSP,Nvar,NvnS[0],Ihat_vS_vS[0],RES,RESP);
					} else {
						L2hat_vS_vS = OPS->L2hat_vS_vS;
						mm_CTN_d(NvnSP,Nvar,NvnS[0],L2hat_vS_vS[0],What,WhatP);
						mm_CTN_d(NvnSP,Nvar,NvnS[0],L2hat_vS_vS[0],RES,RESP);
					}

					free(What);
					free(RES);

					VOLUME->What = WhatP;
					VOLUME->RES  = RESP;
				}
				break;
			case HREFINE:
				VType = VOLUME->type;
				uhat  = VOLUME->uhat;
				What  = VOLUME->What;
				RES   = VOLUME->RES;


				NvnGs[0]     = OPS->NvnGs;
				NvnGc[0]     = OPS->NvnGc;
				NvnS[0]      = OPS->NvnS;
				I_vGs_vGs    = OPS->I_vGs_vGs;
				I_vGs_vGc[0] = OPS->I_vGs_vGc;

				NCols = d;

				VOLUME->hrefine_type = 0;

				get_vh_range(VOLUME,&vhMin,&vhMax);
				for (vh = vhMin; vh <= vhMax; vh++) {
					if (vh == vhMin) {
						VOLUMEc = New_VOLUME();
						VOLUME->child0 = VOLUMEc;
					} else {
						VOLUMEc->next = New_VOLUME();
						VOLUMEc = VOLUMEc->next;
					}
					VOLUMEc->update = 1;
					VOLUMEc->parent = VOLUME;
					VOLUMEc->indexg = NV++;

					VOLUMEc->P    = VOLUME->P;
					VOLUMEc->PNew = VOLUME->P;
					VOLUMEc->level = (VOLUME->level)+1;
					switch (VType) {
					default: // LINE, TRI, QUAD, HEX, WEDGE
						VOLUMEc->type = VType;
						break;
					case TET:
					case PYR:
						VOLUMEc->type = get_VOLUMEc_type(VType,vh);
						break;
					}
					ELEMENT = get_ELEMENT_type(VOLUMEc->type);
					Nf = ELEMENT->Nf;
					for (f = 0; f < Nf; f++)
						VOLUMEc->NsubF[f] = 1;

					VOLUMEc->Eclass = get_Eclass(VOLUMEc->type);

					if (AC) {
						VOLUMEc->curved = 1;
					} else if (VOLUME->curved) {
						printf("Error: Add support for h-refinement VOLUMEc->curved.\n"), EXIT_MSG;
						// Use VToBC and knowledge of whether the new VOLUME shares the BC.
					} else {
						VOLUMEc->curved = 0;
					}

					// Update geometry
					IndEhref = get_IndEhref(VType,vh);

// When updating XYZ_vC, ensure that corners on curved boundaries are placed on the boundary. (ToBeDeleted)
					VOLUMEc->XYZ_vC = malloc(NvnGs[IndEhref]*d * sizeof *XYZ_vC); // keep
					XYZ_vC = VOLUMEc->XYZ_vC;

					mm_CTN_d(NvnGs[IndEhref],NCols,NvnGs[0],I_vGs_vGs[vh],VOLUME->XYZ_vC,VOLUMEc->XYZ_vC);
					if (!VOLUMEc->curved) {
						double *XYZ;

						VOLUMEc->NvnG = NvnGs[IndEhref];

						VOLUMEc->XYZ_S = malloc(NvnGs[IndEhref]*NCols * sizeof *XYZ_S); // keep
						VOLUMEc->XYZ   = malloc(NvnGs[IndEhref]*NCols * sizeof *XYZ);   // keep
						XYZ_S = VOLUMEc->XYZ_S;
						XYZ   = VOLUMEc->XYZ;
						for (unsigned int i = 0, iMax = NCols*NvnGs[IndEhref]; i < iMax; i++) {
							XYZ_S[i] = XYZ_vC[i];
							XYZ[i]   = XYZ_S[i];
						}
					} else {
						VOLUMEc->NvnG = NvnGc[IndEhref];

						VOLUMEc->XYZ_S = malloc(NvnGc[IndEhref]*NCols * sizeof *XYZ_S); // keep
						mm_CTN_d(NvnGc[IndEhref],NCols,NvnGs[IndEhref],I_vGs_vGc[IndEhref],XYZ_vC,VOLUMEc->XYZ_S);

						if (strstr(MeshType,"ToBeCurved"))
							setup_ToBeCurved(VOLUMEc);
						else if (strstr(MeshType,"Curved"))
							printf("Add in support for MeshType == Curved.\n"), EXIT_MSG;
					}
					setup_geom_factors(VOLUMEc);

// Fix Vgrp linked list (ToBeDeleted)

					// Project solution coefficients (and RES if applicable)
					Ihat_vS_vS = OPS->Ihat_vS_vS;
					VOLUMEc->NvnS = NvnS[IndEhref];
					if (strstr(TestCase,"Poisson")) {
						uhatH = malloc(NvnS[IndEhref]*Nvar * sizeof *uhatH); // keep

						mm_CTN_d(NvnS[IndEhref],Nvar,NvnS[0],Ihat_vS_vS[vh],uhat,uhatH);

						VOLUMEc->uhat = uhatH;
					} else {
						WhatH = malloc(NvnS[IndEhref]*Nvar * sizeof *WhatH); // keep
						RESH  = malloc(NvnS[IndEhref]*Nvar * sizeof *RESH);  // keep

						mm_CTN_d(NvnS[IndEhref],Nvar,NvnS[0],Ihat_vS_vS[vh],What,WhatH);
						mm_CTN_d(NvnS[IndEhref],Nvar,NvnS[0],Ihat_vS_vS[vh],RES,RESH);

						VOLUMEc->What = WhatH;
						VOLUMEc->RES  = RESH;
					}
				}
				if (strstr(TestCase,"Poisson")) {
					free(VOLUME->uhat);
				} else {
					free(VOLUME->What);
					free(VOLUME->RES);
				}
				break;
			case HCOARSE:
				VOLUMEp = VOLUME->parent;

				if (VOLUME == VOLUMEp->child0) {
					level = VOLUME->level;

					update = 1;
					maxP = VOLUME->P;

					get_vh_range(VOLUMEp,&vhMin,&vhMax);
					VOLUMEc = VOLUME;
					for (vh = vhMin+1; vh <= vhMax; vh++) {
						VOLUMEc = VOLUMEc->next;
						maxP = max(maxP,VOLUMEc->P);

						if (VOLUMEc->level != level || !VOLUMEc->Vadapt || VOLUMEc->adapt_type != HCOARSE) {
							update = 0;
							break;
						}
					}
					VOLUMEp->update = update;

					if (update) {
						// Most of the information of VOLUMEp has been stored.
						VOLUMEp->update = 1;
						VOLUMEp->indexg = NV++;

						VOLUMEp->P    = maxP;
						VOLUMEp->PNew = maxP;
						VOLUMEp->adapt_type = HCOARSE;

						// Project solution coefficients (and RES if applicable)
						NvnS[0] = OPS->NvnS;
						L2hat_vS_vS = OPS->L2hat_vS_vS;

						NCols = d;

						dummyPtr_d = malloc(NvnS[0]*Nvar * sizeof *dummyPtr_d); // free

						uhat = What = RES = NULL;
						if (strstr(TestCase,"Poisson")) {
							uhat = calloc(NvnS[0]*Nvar , sizeof *uhat); // keep
						} else {
							What = calloc(NvnS[0]*Nvar , sizeof *What); // keep
							RES  = calloc(NvnS[0]*Nvar , sizeof *RES);  // keep
						}

						VOLUMEc = VOLUME;
						for (vh = vhMin; vh <= vhMax; vh++) {
							IndEhref = get_IndEhref(VOLUMEp->type,vh);
							if (vh > vhMin)
								VOLUMEc = VOLUMEc->next;

							if (strstr(TestCase,"Poisson")) {
								uhatH = VOLUMEc->uhat;
								mm_CTN_d(NvnS[0],Nvar,NvnS[IndEhref],L2hat_vS_vS[vh],uhatH,dummyPtr_d);
								for (i = 0, iMax = NvnS[0]*Nvar; i < iMax; i++)
									uhat[i] += dummyPtr_d[i];
							} else {
								WhatH = VOLUMEc->What;
								RESH  = VOLUMEc->RES;

								mm_CTN_d(NvnS[0],Nvar,NvnS[IndEhref],L2hat_vS_vS[vh],WhatH,dummyPtr_d);
								for (i = 0, iMax = NvnS[0]*Nvar; i < iMax; i++)
									What[i] += dummyPtr_d[i];
								mm_CTN_d(NvnS[0],Nvar,NvnS[IndEhref],L2hat_vS_vS[vh],RESH,dummyPtr_d);
								for (i = 0, iMax = NvnS[0]*Nvar; i < iMax; i++)
									RES[i] += dummyPtr_d[i];
							}
						}
						free(dummyPtr_d);

						if (strstr(TestCase,"Poisson")) {
							VOLUMEp->uhat = uhat;
						} else {
							VOLUMEp->What = What;
							VOLUMEp->RES  = RES;
						}
					} else {
						// Ensure that all children are marked as not to be coarsened.
						VOLUMEc = VOLUME;
						for (vh = vhMin; vh <= vhMax; vh++) {
							if (VOLUMEc->adapt_type == HCOARSE) {
								VOLUMEc->Vadapt = 0;
								VOLUMEc->update = 0;
							}
							VOLUMEc = VOLUMEc->next;
						}
					}
				} else {
					VOLUMEc = VOLUMEp->child0;
					if (!(VOLUMEc->adapt_type == HCOARSE && VOLUMEc->Vadapt)) {
						VOLUME->Vadapt = 0;
						VOLUME->update = 0;
					}
				}
				break;
			}
		}
	}
	free(OPS);
}

void update_VOLUME_list(void)
{
	/*
	 *	Comments:
	 *		This is done inelegantly because it was required to keep the pointer to the parent VOLUME.
	 */

	unsigned int adapt_type;

	struct S_VOLUME *VOLUME, *VOLUMEc, *VOLUMEnext;

	// Fix list head if necessary
	VOLUME = DB.VOLUME;

	if (VOLUME->update) {
		adapt_type = VOLUME->adapt_type;
		if (adapt_type == HREFINE) {
			DB.VOLUME = VOLUME->child0;
			for (VOLUMEc = DB.VOLUME; VOLUMEc->next; VOLUMEc = VOLUMEc->next)
				;
			VOLUMEc->next = VOLUME->next;
		} else if (adapt_type == HCOARSE) {
			DB.VOLUME = VOLUME->parent;
			for (VOLUMEc = VOLUME; VOLUMEc->next && VOLUMEc->next->parent == DB.VOLUME; VOLUMEc = VOLUMEc->next)
				;
			DB.VOLUME->next = VOLUMEc->next;
			VOLUMEc->next = NULL;
		}
	}

	// Fix remainder of list
	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		VOLUMEnext = VOLUME->next;
		if (VOLUMEnext && VOLUMEnext->update) {
			adapt_type = VOLUMEnext->adapt_type;
			if (adapt_type == HREFINE) {
				VOLUME->next = VOLUMEnext->child0;
				for (VOLUMEc = VOLUME->next; VOLUMEc->next; VOLUMEc = VOLUMEc->next)
					;
				VOLUMEc->next = VOLUMEnext->next;
			} else if (adapt_type == HCOARSE) {
				VOLUME->next = VOLUMEnext->parent;
				for (VOLUMEc = VOLUMEnext; VOLUMEc->next && VOLUMEc->next->parent == VOLUME->next; VOLUMEc = VOLUMEc->next)
					;
				VOLUME->next->next = VOLUMEc->next;
				VOLUMEc->next = NULL;
			}
		}
	}
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

		for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			Eclass = VOLUME->Eclass;
			P      = VOLUME->P;
			curved = VOLUME->curved;

			IndVgrp = Eclass*NP*2 + P*2;
			if (curved)
				IndVgrp += 1;

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

void compute_inverse_mass(struct S_VOLUME *VOLUME)
{
	// Standard datatypes
	unsigned int iMax, jMax,
	             NvnS, NvnI;
	double       *detJV_vI, *detJV_vI_ptr, *w_vI, *w_vI_ptr, *wdetJV_vI, *wdetJV_vI_ptr, *ChiS_vI, *ChiS_vI_ptr,
	             *wdetJVChiS_vI, *wdetJVChiS_vI_ptr, *IS, *M, *MInv;

	struct S_OPERATORS *OPS;

	OPS = malloc(sizeof *OPS); // free

	init_ops(OPS,VOLUME,0);

	NvnS = OPS->NvnS;
	NvnI = OPS->NvnI;
	w_vI    = OPS->w_vI;
	ChiS_vI = OPS->ChiS_vI;

	detJV_vI = VOLUME->detJV_vI;

	// Compute required portion of MInv
	// ToBeModified (Extra work being done here for collocated)
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
	IS   = identity_d(NvnS);                                                   // free
	MInv = inverse_d(NvnS,NvnS,M,IS);                                          // keep

	free(wdetJVChiS_vI);
	free(M);
	free(IS);
	free(wdetJV_vI);

	free(VOLUME->MInv);
	VOLUME->MInv = MInv;

	free(OPS);
}

void update_VOLUME_Ops(void)
{
	// Initialize DB Parameters
	char         *SolverType = DB.SolverType;
	unsigned int Collocated  = DB.Collocated;

	// Standard datatypes
	struct S_VOLUME    *VOLUME;

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		if (VOLUME->update) {
			VOLUME->update = 0;
			if (strstr(SolverType,"Explicit")) {
				if (!Collocated)
					compute_inverse_mass(VOLUME);
			} else {
				// Do nothing
//				printf("Error: Unsupported SolverType.\n"), EXIT_MSG;
			}
		}
	}
}

void update_VOLUME_finalize(void)
{
	unsigned int NV = 0;
	unsigned int VfIn, VfOut, fIn, fOut;

	struct S_VOLUME *VOLUME, *VIn, *VOut;
	struct S_FACET  *FACET;

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		VOLUME->indexg = NV++;
		VOLUME->Vadapt = 0;
		VOLUME->adapt_type = UINT_MAX; // Done for debugging only
		VOLUME->update = 0;
	}

	DB.NV = NV;
	DB.NVglobal = NV;

	for (FACET = DB.FACET; FACET; FACET = FACET->next) {
		VIn   = FACET->VIn;
		VfIn  = FACET->VfIn;
		fIn   = VfIn/NFREFMAX;

		VOut  = FACET->VOut;
		VfOut = FACET->VfOut;
		fOut  = VfOut/NFREFMAX;

		FACET->Boundary = !((VIn->indexg != VOut->indexg) || (VIn->indexg == VOut->indexg && fIn != fOut));

		VIn->neigh[VfIn]   = VOut->indexg;
		VOut->neigh[VfOut] = VIn->indexg;

		if (fabs((int) VIn->level - (int) VOut->level) > 1.0) {
			printf("%d %d %d\n",VIn->indexg,VOut->indexg,VIn->parent->indexg);
			printf("Error: Adjacent VOLUMEs are more than 1-irregular.\n"), EXIT_MSG;
		}
	}
}
