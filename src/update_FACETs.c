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
 *		H-adaptation FACET updating is done from the coarsest to finest for HREFINE and in the opposite sense for
 *		HCOARSE in order to ensure that neighbouring elements are no more than 1-irregular.
 *
 *	Notation:
 *
 *	References:
 */

static void get_FACET_IndVIn(const unsigned int Vf, const unsigned int fh, const unsigned int VType, unsigned int
*IndVInh, unsigned int *Vfh)
{
	/*
	 *	Comments:
	 *		Vfh != f for TET/PYR refinement.
	 */

	// Standard datatypes
	unsigned int f;

	f = Vf/NFREFMAX;

	switch (VType) {
	case TRI:
		switch (f) {
		default: // FACE 0
			switch (fh) {
			default: // fh == 0
				*IndVInh = 1; break;
			case 1:
				*IndVInh = 2; break;
			}
			*Vfh = 0;
			break;
		case 1:
			switch (fh) {
			default: // fh == 0
				*IndVInh = 0; break;
			case 1:
				*IndVInh = 2; break;
			}
			*Vfh = 1;
			break;
		case 2:
			switch (fh) {
			default: // fh == 0
				*IndVInh = 0; break;
			case 1:
				*IndVInh = 1; break;
			}
			*Vfh = 2;
			break;
		}
		break;
	default:
		printf("Error: Unsupported VType in get_FACET_IndVIn.\n");
		break;
	}
}

static unsigned int get_FACET_VfOut(const unsigned int fh, const unsigned int IndOrd, const unsigned int neigh_f,
                                    const unsigned int VType)
{
	switch (VType) {
	case TRI:
		// Isotropic refinement only
		switch (IndOrd) {
		default: // case 0
			return neigh_f*NFREFMAX+fh+1;
			break;
		case 1:
			return neigh_f*NFREFMAX+((fh+1)%2)+1;
			break;
		}
		break;
	default:
		printf("Error: Unsupported VType in get_VfOut.\n"), exit(1);
		break;
	}
}

static unsigned int get_FACET_type(struct S_FACET *FACET)
{
	// Initialize DB Parameters
	const unsigned int d = DB.d;

	// Standard datatypesc
	unsigned int VType, VfIn, fIn;

	struct S_VOLUME *VIn = FACET->VIn;

	VType = VIn->type;
	VfIn  = FACET->VfIn;
	fIn   = VfIn/NFREFMAX;

	switch (d) {
	default: // d = 3
		if (VType == HEX || (VType == WEDGE && fIn < 3) || (VType == PYR && fIn > 3))
			return QUAD;
		else
			return TRI;
		break;
	case 2:
		return LINE;
		break;
	case 1:
		return POINT;
		break;
	}
}

static unsigned int get_fhMax(const unsigned int VType, const unsigned int href_type)
{
	switch (VType) {
	case TRI:
		// Supported href_type: 0 (Isotropic)
		return 2;
		break;
	default:
		printf("Error: Unsupported VType in get_fhMax.\n"), exit(1);
		return 0;
		break;
	}
}

static void get_Indsf(struct S_FACET *FACET, unsigned int *sfIn, unsigned int *sfOut);

static void set_FACET_Out(const unsigned int vh, struct S_FACET *FACETh, struct S_VOLUME *VOLUME)
{
	unsigned int i, VType, IndVhOut, f, IndOrdInOut, IndOrdOutIn, sfIn, sfOut;

	struct S_VOLUME *VOLUMEh;

printf("%d %d\n",VOLUME->indexg,vh);
	VType = VOLUME->type;
	switch (VType) {
	case TRI:
		// Isotropic refinement only.
		IndVhOut = 3;
		f = vh;
		IndOrdInOut = 1; // Reversed
		IndOrdOutIn = 1; // Reversed
		if (vh == 3)
			printf("Error: Should not be entering set_FACET_Out for vh %d for VType %d.\n",vh,VType), exit(1);
			// Should already have found all FACETs
		break;
	default:
		printf("Error: Unsupported VType in set_FACET_Out.\n"), exit(1);
		break;
	}

	VOLUMEh = VOLUME->child0;
	for (i = 0; i < IndVhOut; i++)
		VOLUMEh = VOLUMEh->next;

	FACETh->VOut  = VOLUMEh;
	FACETh->VfOut = f*NFREFMAX;

	FACETh->IndOrdInOut = IndOrdInOut;
	FACETh->IndOrdOutIn = IndOrdOutIn;

	get_Indsf(FACETh,&sfIn,&sfOut);
//printf("sf: %d %d\n",sfIn,sfOut);
	FACETh->VIn->FACET[sfIn] = FACETh;
	FACETh->VOut->FACET[sfOut] = FACETh;
	FACETh->type = get_FACET_type(FACETh);
}

static void set_FACET_Out_External(struct S_FACET *FACETh, struct S_VOLUME *VOLUME)
{
	unsigned int i, VType, IndVhOut, f;

	struct S_VOLUME *VOLUMEh;

	VType = VOLUME->type;
	switch (VType) {
	case TRI:
		// Isotropic refinement only.
		switch (FACETh->VfOut) {
		case 2:
		case NFREFMAX+2:
			IndVhOut = 2;
			break;
		case NFREFMAX+1:
		case 2*NFREFMAX+1:
			IndVhOut = 0;
			break;
		case 1:
		case 2*NFREFMAX+2:
			IndVhOut = 1;
			break;
		default:
			printf("Error: Unsupported VfOut in set_FACET_Out_External.\n"), exit(1);
			break;
		}
		f = (FACETh->VfOut)/NFREFMAX;
		break;
	default:
		printf("Error: Unsupported VType in set_FACET_Out_External.\n"), exit(1);
		break;
	}

	VOLUMEh = VOLUME->child0;
	for (i = 0; i < IndVhOut; i++)
		VOLUMEh = VOLUMEh->next;

	FACETh->VOut  = VOLUMEh;
	FACETh->VfOut = f*NFREFMAX;

}

static void get_Indsf(struct S_FACET *FACET, unsigned int *sfIn, unsigned int *sfOut)
{
	unsigned int VType, Vf, f, Vfl;

	Vf   = FACET->VfIn;
	f    = Vf/NFREFMAX;
	Vfl  = Vf%NFREFMAX;

	VType = FACET->VIn->type;
	switch (VType) {
	case TRI:
		// Isotropic refinement only.
		switch (Vfl) {
		case 0:
		case 1:
			*sfIn = 0;
			break;
		case 2:
			*sfIn = 1;
			break;
		default:
			printf("Error: Unsupported VfInl in case %d of get_Indsf.\n",VType), exit(1);
			break;
		}
		*sfIn += f*NSUBFMAX;
		break;
	default:
		printf("Error: Unsupported VIn->type in get_Indsf.\n"), exit(1);
		break;
	}

	Vf   = FACET->VfOut;
	f    = Vf/NFREFMAX;
	Vfl  = Vf%NFREFMAX;

	VType = FACET->VOut->type;
	switch (VType) {
	case TRI:
		// Isotropic refinement only.
		switch (Vfl) {
		case 0:
		case 1:
			*sfOut = 0;
			break;
		case 2:
			*sfOut = 1;
			break;
		default:
			printf("Error: Unsupported VfOutl in case %d of get_Indsf.\n",VType), exit(1);
			break;
		}
		*sfOut += f*NSUBFMAX;
		break;
	default:
		printf("Error: Unsupported VOut->type in get_Indsf.\n"), exit(1);
		break;
	}
}

static unsigned int is_VOLUME_VIn(const unsigned int indexgVOLUME, const unsigned int indexgVIn)
{
	if (indexgVOLUME == indexgVIn)
		return 1;
	return 0;
}

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
	unsigned int Adapt     = DB.Adapt,
	             NGF       = DB.NGF,
	             LevelsMax = DB.LevelsMax;

	// Standard datatypes
	unsigned int l, Nf, f, fh, Vfh, vh, sf, sfIn, sfOut, fhMax, sfMax, dummy_ui,
	             Indf, IndVInh, adapt_type,
	             *NsubF;

	struct S_ELEMENT *ELEMENT;
	struct S_VOLUME  *VOLUME, *VOLUMEh, *VIn, *VOut;
	struct S_FACET   *FACET, *FACETh, *FACETnext;

	// silence
	IndVInh = Vfh = 0;
	FACETh = NULL;

	switch (Adapt) {
	default: // ADAPT_HP
		break;
	case ADAPT_H:
// Different for HREFINE and HCOARSE
		for (l = 0; l < LevelsMax; l++) {
			for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
				if (VOLUME->update && VOLUME->adapt_type == HREFINE && l == VOLUME->level) {
					NsubF = VOLUME->NsubF;

					// External FACETs
					ELEMENT = get_ELEMENT_type(VOLUME->type);
					Nf = ELEMENT->Nf;

					for (f = 0; f < Nf; f++) {
						Indf = f*NSUBFMAX;
						sfMax = NsubF[f];
						if (sfMax == 1) { // Create new FACETs between two VOLUMEs which were previously on the same level.
							FACET = VOLUME->FACET[Indf];
//printf("208 %d\n",FACET->indexg);
							FACET->update = 1;
							FACET->adapt_type = HREFINE;

							VIn   = FACET->VIn;
							VOut  = FACET->VOut;

							fhMax = get_fhMax(VOLUME->type,VOLUME->hrefine_type);
							for (fh = 0; fh < fhMax; fh++) {
								if (fh) {
									FACETh->next = New_FACET();
									FACETh = FACETh->next;
								} else {
									FACETh = New_FACET();
									FACET->child0 = FACETh;
								}
								FACETh->update = 1;
								FACETh->parent = FACET;

								FACETh->indexg = NGF++;
								FACETh->level  = (VOLUME->level)+1;
								FACETh->BC     = FACET->BC;

								// Find out if VOLUME == VIn or VOut
								if (is_VOLUME_VIn(VOLUME->indexg,VIn->indexg)) {
// If condition can be outside of fh loop
									get_FACET_IndVIn(FACET->VfIn,fh,VIn->type,&IndVInh,&Vfh);

									VOLUMEh = VOLUME->child0;
									for (vh = 0; vh < IndVInh; vh++)
										VOLUMEh = VOLUMEh->next;

									FACETh->VIn  = VOLUMEh;
									FACETh->VfIn = Vfh*NFREFMAX;

									FACETh->VOut  = VOut;
									FACETh->VfOut = get_FACET_VfOut(fh,FACET->IndOrdOutIn,VIn->neigh_f[f*NFREFMAX],VIn->type);

									// Valid for TET/PYR? (ToBeDeleted)
									FACETh->IndOrdInOut = FACET->IndOrdInOut;
									FACETh->IndOrdOutIn = FACET->IndOrdOutIn;

									// Update VOLUME related parameters
// Add to else below.
// Potentially add VOLUME->(neigh,neigh_f,FACET)
									VOut->NsubF[(FACETh->VfOut)/NFREFMAX] = fhMax;
								} else {
									get_FACET_IndVIn(FACET->VfOut,fh,VOut->type,&IndVInh,&Vfh);

									VOLUMEh = VOLUME->child0;
									for (vh = 0; vh < IndVInh; vh++)
										VOLUMEh = VOLUMEh->next;

									FACETh->VIn  = VOLUMEh;
									FACETh->VfIn = Vfh*NFREFMAX;

									FACETh->VOut  = VIn;
									FACETh->VfOut = get_FACET_VfOut(fh,FACET->IndOrdInOut,VOut->neigh_f[f*NFREFMAX],VOut->type);

									FACETh->IndOrdInOut = FACET->IndOrdOutIn;
									FACETh->IndOrdOutIn = FACET->IndOrdInOut;
								}
								get_Indsf(FACETh,&sfIn,&sfOut);
//printf("sf: %d %d %d %d\n",FACETh->VfIn,FACETh->VfOut,sfIn,sfOut);
								FACETh->VIn->FACET[sfIn] = FACETh;
								FACETh->VOut->FACET[sfOut] = FACETh;
								FACETh->type = get_FACET_type(FACETh);
							}
						} else { // Connect to existing FACETs
							// VOLUME = VOut
							for (sf = 0; sf < sfMax; sf++) {
								FACETh = VOLUME->FACET[Indf+sf];
								printf("%d %d %d\n",FACETh->indexg,VOLUME->indexg,f);
								printf("%d\n",FACETh->VfOut);

								// Only connectivity needs updating in FACETh->VOut
								set_FACET_Out_External(FACETh,VOLUME);
							}
//printf("Exiting\n");
//exit(1);

						}
					}

					// Internal FACETs (Always created)
					vh = 0;
					for (VOLUMEh = VOLUME->child0; VOLUMEh != NULL; VOLUMEh = VOLUMEh->next) {
						ELEMENT = get_ELEMENT_type(VOLUMEh->type);
						Nf = ELEMENT->Nf;

						for (f = 0; f < Nf; f++) {
							if (!VOLUMEh->FACET[f*NSUBFMAX]) {
								FACETh->next = New_FACET();
								FACETh = FACETh->next;
								FACETh->update = 1;
								FACETh->parent = NULL;
								FACETh->indexg = NGF++;
								FACETh->level  = (VOLUME->level)+1;
								FACETh->BC     = 0; // Internal (Note: May have a curved edge in 3D)

								VOLUMEh->NsubF[f] = 1;

								FACETh->VIn = VOLUMEh;
								FACETh->VfIn = f*NFREFMAX;

								set_FACET_Out(vh,FACETh,VOLUME);
							}
						}
						vh++;
					}
//exit(1);
				}
			}
		}
		DB.NGF = NGF;

		// Update FACET linked list (For HREFINE)
// This is done inelegantly because it was required to keep the pointer to the parent FACET.

		// Fix list head if necessary
		FACET = DB.FACET;

		if (FACET->update) {
			adapt_type = FACET->adapt_type;
			if (adapt_type == HREFINE) {
				DB.FACET = FACET->child0;
				for (FACETh = DB.FACET; FACETh->next != NULL; FACETh = FACETh->next)
					;
				FACETh->next = FACET->next;
			} else if (adapt_type == HCOARSE) {
				DB.FACET = FACET->parent;
				for (FACETh = FACET; FACETh->parent == DB.FACET; FACETh = FACETh->next)
					FACETh->update = 0;
				DB.FACET->next = FACETh;
			}
		}

		// Fix remainder of list
		for (FACET = DB.FACET; FACET != NULL; FACET = FACET->next) {
			FACETnext = FACET->next;
			if (FACETnext && FACETnext->update) {
				adapt_type = FACETnext->adapt_type;
				if (adapt_type == HREFINE) {
					FACET->next = FACETnext->child0;
					for (FACETh = FACET->next; FACETh->next != NULL; FACETh = FACETh->next)
						;
					FACETh->next = FACETnext->next;
				} else if (adapt_type == HCOARSE) {
					FACET->next = FACETnext->parent;
					for (FACETh = FACETnext; FACETh->next->parent == FACET->next; FACETh = FACETh->next)
						FACETh->update = 0;
					FACET->next->next = FACETh->next;
				}
			}
		}

		for (FACET = DB.FACET; FACET != NULL; FACET = FACET->next) {
			if (FACET->update) {
				FACET->update = 0;
				if (FACET->parent && FACET->parent->update)
					FACET->parent->update = 0;
				FACET->P = max(FACET->VIn->P,FACET->VOut->P);

				// Compute XYZ_fS, n_fS, and detJF_fS
				setup_FACET_XYZ(FACET);
				setup_normals(FACET);
			}
		}
//exit(1);

// START HCOARSE HERE
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
