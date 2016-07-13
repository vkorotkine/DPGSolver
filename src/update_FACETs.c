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
//	if (BC == 0 || (BC % BC_STEP_SC > 50)) {
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
//	} else {
//	}
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

	struct S_VOLUME *VOLUMEc;

//printf("%d %d\n",VOLUME->indexg,vh);
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

	VOLUMEc = VOLUME->child0;
	for (i = 0; i < IndVhOut; i++)
		VOLUMEc = VOLUMEc->next;

	FACETh->VOut  = VOLUMEc;
	FACETh->VfOut = f*NFREFMAX;

	FACETh->IndOrdInOut = IndOrdInOut;
	FACETh->IndOrdOutIn = IndOrdOutIn;

	get_Indsf(FACETh,&sfIn,&sfOut);
//printf("sf: %d %d\n",sfIn,sfOut);
	FACETh->VIn->FACET[sfIn] = FACETh;
	FACETh->VOut->FACET[sfOut] = FACETh;
	FACETh->type = get_FACET_type(FACETh);

	FACETh->VIn->neigh_f[FACETh->VfIn] = (FACETh->VfOut)/NFREFMAX;
	FACETh->VOut->neigh_f[FACETh->VfOut] = (FACETh->VfIn)/NFREFMAX;
}

static void set_FACET_Out_External(struct S_FACET *FACETh, struct S_VOLUME *VOLUME)
{
	unsigned int i, VType, IndVhOut, f, sfIn, sfOut;

	struct S_VOLUME *VOLUMEc;

	VType = VOLUME->type;
	switch (VType) {
	case TRI:
		// Isotropic refinement only.
		switch (FACETh->VfOut) {
		case NFREFMAX+1:
		case 2*NFREFMAX+1:
			IndVhOut = 0;
			break;
		case 1:
		case 2*NFREFMAX+2:
			IndVhOut = 1;
			break;
		case 2:
		case NFREFMAX+2:
			IndVhOut = 2;
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

	VOLUMEc = VOLUME->child0;
	for (i = 0; i < IndVhOut; i++)
		VOLUMEc = VOLUMEc->next;

	FACETh->VOut  = VOLUMEc;
	FACETh->VfOut = f*NFREFMAX;
//printf("\t\t\t\tset_FACET_Ex: %d %d %d\n",VOLUME->indexg,VOLUMEc->indexg-8,FACETh->VfOut);
	FACETh->VOut->neigh_f[FACETh->VfOut] = (FACETh->VfIn)/NFREFMAX;

	get_Indsf(FACETh,&sfIn,&sfOut);
	VOLUMEc->FACET[sfOut] = FACETh;
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

static void coarse_update(const struct S_VOLUME *VOLUME)
{
// Use this routine to update connectivity as well instead of redoing the set up. => Change function name (ToBeDeleted)
	unsigned int i, iMax, f, Nf, VType, vh, vhMin, vhMax, IndVc[NSUBFMAX],IndVf[NSUBFMAX], recombine_FACETs;

	struct S_ELEMENT *ELEMENT;
	struct S_VOLUME  *VOLUMEc, *VOLUMEp, **VOLUMEc_list;

	ELEMENT = get_ELEMENT_type(VOLUMEp->type);
	Nf = ELEMENT->Nf;

	get_vh_range(VOLUME,&vhMin,&vhMax);
	iMax = vhMax-vhMin+1;
	VOLUMEc_list = malloc(iMax * sizeof *VOLUMEc_list); // free

	VOLUMEc = VOLUME->child0;
	for (i = 0; i < iMax; i++) {
		VOLUMEc_list[i] = VOLUMEc;
		VOLUMEc = VOLUMEc->next;
	}


	VType = VOLUME->type;
	for (f = 0; f < Nf; f++) {
		switch (VType) {
		case TRI:
			// Supported: Isotropic refinement
			sfMax = 2;
			switch (f) {
			default: // f = 0
				IndVc[0] = 1;          IndVc[1] = 2;
				IndVf[0] = 0*NFREFMAX; IndVf[1] = 0*NFREFMAX;
				break;
			case 1:
				IndVc[0] = 1;          IndVc[1] = 2;
				IndVf[0] = 1*NFREFMAX; IndVf[1] = 1*NFREFMAX;
				break;
			case 2:
				IndVc[0] = 0;          IndVc[1] = 1;
				IndVf[0] = 2*NFREFMAX; IndVf[1] = 2*NFREFMAX;
				break;
			}
			break;
		default:
			printf("Error: Unsupported VType in recombine_FACETs.\n"), exit(1);
			break;
		}

		recombine_FACETs = 0;
// Only need to check sf = 0 here as all sub FACETs should have the same level gap if coarsening was allowed.
// (ToBeDeleted).
		for (sf = 0; sf < sfMax; sf++) {
			FACET = VOLUMEc_list[IndVc[sf]]->FACET[IndVf[sf]];
			if (FACET->VIn->level != FACET->VOut->level) {
				recombine_FACETs = 1;
				break;
			}
		}

		if (recombined_FACETs) {
// Already initialized?			VOLUMEp->FACET[f*NFREFMAX] = FACET->parent;
			VOLUMEp->NsubF[f] = 1;

			for (sf = 0; sf < sfMax; sf++) {
				VOLUMEc_list[IndVc[sf]]->FACET[IndVf[sf]]->update = 1;


		} else {
			VOLUMEp->NsubF[f] = sfMax;

		}
	}







	free(VOLUMEc_list);
	return recombine_FACETs;
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
	             Indf, IndVInh, adapt_type, vhMin, vhMax, level, same_level,
	             *NsubF;

	struct S_ELEMENT *ELEMENT;
	struct S_VOLUME  *VOLUME, *VOLUMEc, *VIn, *VOut;
	struct S_FACET   *FACET, *FACETh, *FACETnext;

	// silence
	IndVInh = Vfh = 0;
	FACETh = NULL;

	switch (Adapt) {
	default: // ADAPT_HP
		break;
	case ADAPT_H:
		// HREFINE
		for (l = 0; l < LevelsMax; l++) {
			for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
				if (VOLUME->update && VOLUME->adapt_type == HREFINE && l == VOLUME->level) {
					NsubF = VOLUME->NsubF;

//printf("%d\n",VOLUME->indexg);
					// External FACETs
					ELEMENT = get_ELEMENT_type(VOLUME->type);
					Nf = ELEMENT->Nf;

					for (f = 0; f < Nf; f++) {
						Indf = f*NSUBFMAX;
						sfMax = NsubF[f];
if (VOLUME->indexg == 3) {
//	printf("sfMax: %d %d\n",f,sfMax);
}
						if (sfMax == 1) { // Create new FACETs between two VOLUMEs which were previously on the same level.
							FACET = VOLUME->FACET[Indf];
//printf("342 %d\n",FACET->indexg);
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
//printf("362 %d\n",FACETh->indexg);
								FACETh->level  = (VOLUME->level)+1;
								FACETh->BC     = FACET->BC;

								// Find out if VOLUME == VIn or VOut
if (FACET->level == 1) {
//	printf("%d %d %d %d\n",FACET->indexg,VIn->indexg,VOut->indexg,VOLUME->indexg);
//	printf("%d\n",is_VOLUME_VIn(VOLUME->indexg,VIn->indexg));
//	exit(1);
}
								if (is_VOLUME_VIn(VOLUME->indexg,VIn->indexg)) {
// If condition can be outside of fh loop
									get_FACET_IndVIn(FACET->VfIn,fh,VIn->type,&IndVInh,&Vfh);

									VOLUMEc = VOLUME->child0;
									for (vh = 0; vh < IndVInh; vh++)
										VOLUMEc = VOLUMEc->next;

									FACETh->VIn  = VOLUMEc;
									FACETh->VfIn = Vfh*NFREFMAX;

									FACETh->VOut  = VOut;
									FACETh->VfOut = get_FACET_VfOut(fh,FACET->IndOrdOutIn,VIn->neigh_f[f*NFREFMAX],VIn->type);
									// Valid for TET/PYR? (ToBeDeleted)
									FACETh->IndOrdInOut = FACET->IndOrdInOut;
									FACETh->IndOrdOutIn = FACET->IndOrdOutIn;

									// Update VOLUME related parameters
// Add to else below.
// Potentially add VOLUME->(neigh,neigh_f,FACET)
//									VOut->NsubF[(FACETh->VfOut)/NFREFMAX] = fhMax;
if (VOLUME->indexg == 1 && f == 1) {
//	printf("Vh6: %d %d %d\n",VOut->indexg,FACETh->VfOut,VIn->neigh_f[f*NFREFMAX]);
}
								} else {
									get_FACET_IndVIn(FACET->VfOut,fh,VOut->type,&IndVInh,&Vfh);

									VOLUMEc = VOLUME->child0;
									for (vh = 0; vh < IndVInh; vh++)
										VOLUMEc = VOLUMEc->next;

									FACETh->VIn  = VOLUMEc;
									FACETh->VfIn = Vfh*NFREFMAX;

									FACETh->VOut  = VIn;
									FACETh->VfOut = get_FACET_VfOut(fh,FACET->IndOrdInOut,VOut->neigh_f[f*NFREFMAX],VOut->type);
//									FACETh->VfOut = get_FACET_VfOut(fh,FACET->IndOrdOutIn,VOLUMEc->neigh_f[f*NFREFMAX],VOLUME->type);

									FACETh->IndOrdInOut = FACET->IndOrdOutIn;
									FACETh->IndOrdOutIn = FACET->IndOrdInOut;

								}
								FACETh->VOut->NsubF[(FACETh->VfOut)/NFREFMAX] = fhMax;
								VOLUMEc->neigh_f[FACETh->VfIn] = (FACETh->VfOut)/NFREFMAX;

if (FACETh->VOut->indexg == 3) {
//	printf("FVOut3: %d %d %d\n",fh,FACETh->VfOut,(FACETh->VfOut)/NFREFMAX);
}

//printf("%d %d %d %d %d %d\n",FACET->indexg,VIn->indexg,VOut->indexg,FACETh->VfIn,FACETh->VfOut,FACETh->BC);
if (FACETh->indexg == 99) {
//	exit(1);
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
//								printf("Connect: %d %d %d ",FACETh->indexg,VOLUME->indexg,f);
//								printf("%d\n",FACETh->VfOut);

								// Only connectivity needs updating in FACETh->VOut
								set_FACET_Out_External(FACETh,VOLUME);
							}
//printf("Exiting\n");
//exit(1);
						}
					}

					// Internal FACETs (Always created)
					vh = 0;
					for (VOLUMEc = VOLUME->child0; VOLUMEc != NULL; VOLUMEc = VOLUMEc->next) {
						ELEMENT = get_ELEMENT_type(VOLUMEc->type);
						Nf = ELEMENT->Nf;

						for (f = 0; f < Nf; f++) {
							if (!VOLUMEc->FACET[f*NSUBFMAX]) {
								while (FACETh->next)
									FACETh = FACETh->next;
								FACETh->next = New_FACET();
								FACETh = FACETh->next;
								FACETh->update = 1;
								FACETh->parent = NULL;
								FACETh->indexg = NGF++;
								FACETh->level  = (VOLUME->level)+1;
								FACETh->BC     = 0; // Internal (Note: May have a curved edge in 3D)

								VOLUMEc->NsubF[f] = 1;

								FACETh->VIn = VOLUMEc;
								FACETh->VfIn = f*NFREFMAX;

								set_FACET_Out(vh,FACETh,VOLUME);
/*
//if (FACETh->indexg == 48) {
if (VOLUMEc->indexg == 20) {
	printf("*******************\n");
	printf("%d\n",f);
	printf("%d %d\n",FACETh->VIn->indexg,FACETh->VfIn);
	printf("%d %d\n",FACETh->VOut->indexg,FACETh->VfOut);
	printf("*******************\n");
}
*/
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
//printf("FUp: %d ",FACET->indexg);
//if (FACETnext)
//	printf("%d %d ",FACETnext->indexg,FACETnext->update);
//printf("\n");
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

		// HCOARSE
		for (l = LevelsMax; l > 0; l--) {
		for (VOLUMEc = DB.VOLUME; VOLUMEc != NULL; VOLUMEc = VOLUMEc->next) {
		if (VOLUMEc->update && VOLUMEc->adapt_type == HCOARSE && l == VOLUMEc->level) {
		if (VOLUMEc == VOLUMEc->parent->child0) {
			VOLUMEp = VOLUME->parent;

			if (VOLUMEp->update) {
				coarse_update(VOLUMEp);
			}
		}}}}
// Don't forget to update FACET list (ToBeDeleted)
// Don't forget to free children of coarsed VOLUMEs and FACETs (ToBeDeleted)
		// Fix list head if necessary
		FACET = DB.FACET;

		if (FACET->update) {
			adapt_type = FACET->adapt_type;
			if (adapt_type == HCOARSE) {
				DB.FACET = FACET->parent;
				for (FACETh = FACET; FACETh->parent == DB.FACET; FACETh = FACETh->next)
					FACETh->update = 0;
				DB.FACET->next = FACETh;
			}
		}

		// Fix remainder of list
		for (FACET = DB.FACET; FACET != NULL; FACET = FACET->next) {
			FACETnext = FACET->next;
//printf("FUp: %d ",FACET->indexg);
//if (FACETnext)
//	printf("%d %d ",FACETnext->indexg,FACETnext->update);
//printf("\n");
			if (FACETnext && FACETnext->update) {
				adapt_type = FACETnext->adapt_type;
				if (adapt_type == HCOARSE) {
					FACET->next = FACETnext->parent;
					for (FACETh = FACETnext; FACETh->next->parent == FACET->next; FACETh = FACETh->next)
						FACETh->update = 0;
					FACET->next->next = FACETh->next;
				}
			}
		}

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
