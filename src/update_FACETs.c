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
 *		Clean up this function: function names, variable names, comments, potentially combine the COARSE and FINE FACET
 *		list updating (ToBeDeleted)
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

static void set_FACET_Out(const unsigned int vh, struct S_FACET *FACETc, struct S_VOLUME *VOLUME)
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

	FACETc->VOut  = VOLUMEc;
	FACETc->VfOut = f*NFREFMAX;

	FACETc->IndOrdInOut = IndOrdInOut;
	FACETc->IndOrdOutIn = IndOrdOutIn;

	get_Indsf(FACETc,&sfIn,&sfOut);
//printf("sf: %d %d\n",sfIn,sfOut);
	FACETc->VIn->FACET[sfIn] = FACETc;
	FACETc->VOut->FACET[sfOut] = FACETc;
	FACETc->type = get_FACET_type(FACETc);

	FACETc->VIn->neigh_f[FACETc->VfIn]   = (FACETc->VfOut)/NFREFMAX;
	FACETc->VOut->neigh_f[FACETc->VfOut] = (FACETc->VfIn)/NFREFMAX;
}

static void set_FACET_Out_External(struct S_FACET *FACETc, struct S_VOLUME *VOLUME)
{
	unsigned int i, VType, IndVhOut, f, sfIn, sfOut;

	struct S_VOLUME *VOLUMEc;

	VType = VOLUME->type;
	switch (VType) {
	case TRI:
		// Isotropic refinement only.
		switch (FACETc->VfOut) {
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
printf("%d %d\n",FACETc->VIn->indexg,FACETc->VOut->indexg);
			printf("Error: Unsupported VfOut = %d in set_FACET_Out_External.\n",FACETc->VfOut), exit(1);
			break;
		}
		f = (FACETc->VfOut)/NFREFMAX;
		break;
	default:
		printf("Error: Unsupported VType in set_FACET_Out_External.\n"), exit(1);
		break;
	}

	VOLUMEc = VOLUME->child0;
	for (i = 0; i < IndVhOut; i++)
		VOLUMEc = VOLUMEc->next;

	FACETc->VOut  = VOLUMEc;
	FACETc->VfOut = f*NFREFMAX;
//printf("\t\t\t\tset_FACET_Ex: %d %d %d\n",VOLUME->indexg,VOLUMEc->indexg-8,FACETc->VfOut);
	FACETc->VOut->neigh_f[FACETc->VfOut] = (FACETc->VfIn)/NFREFMAX;

	get_Indsf(FACETc,&sfIn,&sfOut);
	VOLUMEc->FACET[sfOut] = FACETc;
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

static void coarse_update(struct S_VOLUME *VOLUME)
{
	unsigned int i, iMax, f, Nf, VType, vhMin, vhMax, sf, sfMax, sfMax_i, fIn, fOut,
	             IndVc[NVISUBFMAX],Indsf[NVISUBFMAX], recombine_FACETs,
	             dummy_ui;

	struct S_ELEMENT *ELEMENT;
	struct S_VOLUME  *VOLUMEc, **VOLUMEc_list, *VIn, *VOut;
	struct S_FACET   *FACET, *FACETp;

	ELEMENT = get_ELEMENT_type(VOLUME->type);
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
				Indsf[0] = 0*NSUBFMAX; Indsf[1] = 0*NSUBFMAX;
				break;
			case 1:
				IndVc[0] = 0;          IndVc[1] = 2;
				Indsf[0] = 1*NSUBFMAX; Indsf[1] = 1*NSUBFMAX;
				break;
			case 2:
				IndVc[0] = 0;          IndVc[1] = 1;
				Indsf[0] = 2*NSUBFMAX; Indsf[1] = 2*NSUBFMAX;
				break;
			}
			break;
		default:
			printf("Error: Unsupported VType in coarse_update.\n"), exit(1);
			break;
		}

		recombine_FACETs = 0;
// Only need to check sf = 0 here as all sub FACETs should have the same level gap if coarsening was allowed.
// (ToBeDeleted).
		for (sf = 0; sf < sfMax; sf++) {
			FACET = VOLUMEc_list[IndVc[sf]]->FACET[Indsf[sf]];
/*
printf("%d %d %p\n",IndVc[sf],Indsf[sf],FACET);
for (int i = 0; i < NFMAX*NSUBFMAX; i++) {
	if (VOLUMEc_list[IndVc[sf]]->FACET[i])
		printf("%d\n",i);
}
*/
			VIn  = FACET->VIn;
			VOut = FACET->VOut;

//printf("Fcoa: %d %d %d %d\n",VOLUME->indexg,f,VIn->level,VOut->level);
			if (VIn == VOut || VIn->level != VOut->level) {
//			if (VIn->level != VOut->level) { Not sure why this is working => check (ToBeDeleted)
				recombine_FACETs = 1;
				break;
			}
		}
//printf("Flag: %d %d %d\n",VOLUME->indexg,f,recombine_FACETs);

// Need to update FACET->VIn and FACET->VOut;
		if (recombine_FACETs) {
			VOLUME->NsubF[f] = 1;

			for (sf = 0; sf < sfMax; sf++) {
				VOLUMEc = VOLUMEc_list[IndVc[sf]];
				FACET = VOLUMEc->FACET[Indsf[sf]];
//printf("Recombup: %d %d %d\n",f,sf,FACET->indexg);
				FACET->update = 1;
				FACET->adapt_type = HCOARSE;
			}

			VOLUME->FACET[f*NSUBFMAX] = FACET->parent;
			FACETp = FACET->parent;
			if (is_VOLUME_VIn(VOLUME->indexg,FACETp->VIn->indexg)) {
				fOut = (FACETp->VfOut)/NFREFMAX;
				FACETp->VOut->NsubF[fOut] = 1;
				FACETp->VOut->FACET[fOut*NSUBFMAX] = FACETp;
			} else {
				fIn = (FACETp->VfIn)/NFREFMAX;
				FACETp->VIn->NsubF[fIn] = 1;
				FACETp->VIn->FACET[fIn*NSUBFMAX] = FACETp;
			}
		} else {
			for (sf = 0; sf < sfMax; sf++) {
				VOLUMEc = VOLUMEc_list[IndVc[sf]];
				FACET = VOLUMEc->FACET[Indsf[sf]];
//				FACET->adapt_type = HCOARSE;

				if (is_VOLUME_VIn(VOLUMEc->indexg,FACET->VIn->indexg)) {
					FACET->VIn  = FACET->VOut;
					FACET->VfIn = FACET->VfOut;

					dummy_ui           = FACET->IndOrdInOut;
					FACET->IndOrdInOut = FACET->IndOrdOutIn;
					FACET->IndOrdOutIn = dummy_ui;

					FACET->VOut  = VOLUME;
					FACET->VfOut = f*NFREFMAX+sf+1; // CHANGE THIS TO BE GENERAL! (ToBeDeleted)

					VOLUME->FACET[f*NSUBFMAX+sf] = FACET;

					// Recompute XYZ_fS, n_fS, and detJF_fS
// Need not be recomputed, just reordered (and potentially negated) (ToBeDeleted)
					free(FACET->XYZ_fS);
					setup_FACET_XYZ(FACET);

					free(FACET->n_fS);
					free(FACET->detJF_fS);
					setup_normals(FACET);
				} else {
					FACET->VOut  = VOLUME;
					FACET->VfOut = f*NFREFMAX+sf+1;

					VOLUME->FACET[f*NSUBFMAX+sf] = FACET;
				}
			}
			VOLUME->NsubF[f] = sfMax;

		}
	}

	// Mark internal FACETs for updating
	// Note: sfMax_i = ((#ELEMENTc)*(Nfc)-(Nfp)*(Nsubfp))/2; (c)hild, (p)arent
	switch (VType) {
	case TRI:
		// Supported: Isotropic refinement
		sfMax_i = 3;

		IndVc[0] = 3; IndVc[1] = 3; IndVc[2] = 3; 
		Indsf[0] = 0; Indsf[1] = 1; Indsf[2] = 2;
		for (i = 0; i < sfMax_i; i++)
			Indsf[i] *= NSUBFMAX;
		break;
	default:
		printf("Error: Unsupported VType in coarse_update (internal).\n"), exit(1);
		break;
	}

	for (sf = 0; sf < sfMax_i; sf++) {
		VOLUMEc = VOLUMEc_list[IndVc[sf]];
		FACET = VOLUMEc->FACET[Indsf[sf]];
		FACET->update = 1;
		FACET->adapt_type = HDELETE;
	}

	free(VOLUMEc_list);
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
	struct S_VOLUME  *VOLUME, *VOLUMEc, *VOLUMEp, *VIn, *VOut;
	struct S_FACET   *FACET, *FACETc, *FACETnext, *FACETtmp;

	// silence
	IndVInh = Vfh = 0;
	FACETc = NULL;

	switch (Adapt) {
	default: // ADAPT_HP
		break;
	case ADAPT_H:
		// HREFINE
		for (l = 0; l < LevelsMax; l++) {
			for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
				if (VOLUME->update && VOLUME->adapt_type == HREFINE && l == VOLUME->level) {
					NsubF = VOLUME->NsubF;

//printf("In HREF: %d %d\n",VOLUME->indexg,VOLUME->adapt_type);
					// External FACETs
					ELEMENT = get_ELEMENT_type(VOLUME->type);
					Nf = ELEMENT->Nf;

					for (f = 0; f < Nf; f++) {
//printf("VOLind,f,Nf,NsubF: %d %d %d %d\n",VOLUME->indexg,f,Nf,NsubF[f]);
						Indf = f*NSUBFMAX;
						sfMax = NsubF[f];
//printf("  f, sfMax: %d %d\n",f,sfMax);
						if (sfMax == 1) { // Create new FACETs between two VOLUMEs which were previously on the same level.
							FACET = VOLUME->FACET[Indf];
//printf("342 %d\n",FACET->indexg);
							FACET->update = 1;
							FACET->adapt_type = HREFINE;
/*
printf("FACETpoint: %p",FACET);
printf("FACETpoint: %p",FACET->VIn);
printf("FACETpoint: %p\n",FACET->VOut);
*/
							VIn   = FACET->VIn;
							VOut  = FACET->VOut;

							fhMax = get_fhMax(VOLUME->type,VOLUME->hrefine_type);
							for (fh = 0; fh < fhMax; fh++) {
								if (fh) {
									FACETc->next = New_FACET();
									FACETc = FACETc->next;
								} else {
									FACETc = New_FACET();
									FACET->child0 = FACETc;
								}
								FACETc->update = 1;
								FACETc->parent = FACET;

								FACETc->indexg = NGF++;
//printf("362 %d\n",FACETc->indexg);
								FACETc->level  = (VOLUME->level)+1;
								FACETc->BC     = FACET->BC;

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

									FACETc->VIn  = VOLUMEc;
									FACETc->VfIn = Vfh*NFREFMAX;

									FACETc->VOut  = VOut;
									FACETc->VfOut = get_FACET_VfOut(fh,FACET->IndOrdOutIn,VIn->neigh_f[f*NFREFMAX],VIn->type);
								// Valid for TET/PYR? (ToBeDeleted)
									FACETc->IndOrdInOut = FACET->IndOrdInOut;
									FACETc->IndOrdOutIn = FACET->IndOrdOutIn;

									// Update VOLUME related parameters
// Add to else below.
// Potentially add VOLUME->(neigh,neigh_f,FACET)
//									VOut->NsubF[(FACETc->VfOut)/NFREFMAX] = fhMax;
if (VOLUME->indexg == 1 && f == 1) {
//	printf("Vh6: %d %d %d\n",VOut->indexg,FACETc->VfOut,VIn->neigh_f[f*NFREFMAX]);
}
								} else {
									get_FACET_IndVIn(FACET->VfOut,fh,VOut->type,&IndVInh,&Vfh);

									VOLUMEc = VOLUME->child0;
									for (vh = 0; vh < IndVInh; vh++)
										VOLUMEc = VOLUMEc->next;

									FACETc->VIn  = VOLUMEc;
									FACETc->VfIn = Vfh*NFREFMAX;

									FACETc->VOut  = VIn;
									FACETc->VfOut = get_FACET_VfOut(fh,FACET->IndOrdInOut,VOut->neigh_f[f*NFREFMAX],VOut->type);
//									FACETc->VfOut = get_FACET_VfOut(fh,FACET->IndOrdOutIn,VOLUMEc->neigh_f[f*NFREFMAX],VOLUME->type);

									FACETc->IndOrdInOut = FACET->IndOrdOutIn;
									FACETc->IndOrdOutIn = FACET->IndOrdInOut;

								}
								FACETc->VIn->NsubF[(FACETc->VfIn)/NFREFMAX] = 1;
								FACETc->VOut->NsubF[(FACETc->VfOut)/NFREFMAX] = fhMax;
								VOLUMEc->neigh_f[FACETc->VfIn] = (FACETc->VfOut)/NFREFMAX;

if (FACETc->VOut->indexg == 3) {
//	printf("FVOut3: %d %d %d\n",fh,FACETc->VfOut,(FACETc->VfOut)/NFREFMAX);
}

//printf("%d %d %d %d %d %d\n",FACET->indexg,VIn->indexg,VOut->indexg,FACETc->VfIn,FACETc->VfOut,FACETc->BC);
if (FACETc->indexg == 99) {
//	exit(1);
}

								get_Indsf(FACETc,&sfIn,&sfOut);
//printf("sf: %d %d %d %d\n",FACETc->VfIn,FACETc->VfOut,sfIn,sfOut);
								FACETc->VIn->FACET[sfIn] = FACETc;
								FACETc->VOut->FACET[sfOut] = FACETc;
								FACETc->type = get_FACET_type(FACETc);
							}
						} else { // Connect to existing FACETs
							// VOLUME = VOut
							for (sf = 0; sf < sfMax; sf++) {
								FACETc = VOLUME->FACET[Indf+sf];
//								printf("Connect: %d %d %d ",FACETc->indexg,VOLUME->indexg,FACETc->level);
//								printf("%d\n",FACETc->VfOut);

								// Only connectivity needs updating in FACETc->VOut
//printf("VOL: %d\n",VOLUME->indexg);
								set_FACET_Out_External(FACETc,VOLUME);
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
								while (FACETc->next)
									FACETc = FACETc->next;
								FACETc->next = New_FACET();
								FACETc = FACETc->next;
								FACETc->update = 1;
								FACETc->indexg = NGF++;
								FACETc->level  = (VOLUME->level)+1;
								FACETc->BC     = 0; // Internal (Note: May have a curved edge in 3D)

								VOLUMEc->NsubF[f] = 1;

								FACETc->VIn = VOLUMEc;
								FACETc->VfIn = f*NFREFMAX;

								set_FACET_Out(vh,FACETc,VOLUME);
/*
//if (FACETc->indexg == 48) {
if (VOLUMEc->indexg == 20) {
	printf("*******************\n");
	printf("%d\n",f);
	printf("%d %d\n",FACETc->VIn->indexg,FACETc->VfIn);
	printf("%d %d\n",FACETc->VOut->indexg,FACETc->VfOut);
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
/*
for (FACET = DB.FACET; FACET; FACET = FACET->next) {
	printf("FTestREF: %d %d %d %d\n",FACET->indexg,FACET->level,FACET->update,FACET->adapt_type);
	if (FACET->indexg == 240)
		printf("\n");
}
printf("\n\n\n");
*/
// No need to look at HCOARSE here (ToBeDeleted).
// This is done inelegantly because it was required to keep the pointer to the parent FACET.

		// Fix list head if necessary
		FACET = DB.FACET;

		if (FACET->update) {
			adapt_type = FACET->adapt_type;
			if (adapt_type == HREFINE) {
				DB.FACET = FACET->child0;
				for (FACETc = DB.FACET; FACETc->next; FACETc = FACETc->next)
					;
				FACETc->next = FACET->next;
//			} else if (adapt_type == HCOARSE) {
//				DB.FACET = FACET->parent;
//				for (FACETc = FACET; FACETc->parent == DB.FACET; FACETc = FACETc->next)
//					FACETc->update = 0;
//				DB.FACET->next = FACETc;
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
					for (FACETc = FACET->next; FACETc->next; FACETc = FACETc->next)
						;
					FACETc->next = FACETnext->next;
//				} else if (adapt_type == HCOARSE) {
//					FACET->next = FACETnext->parent;
//					for (FACETc = FACETnext; FACETc->next->parent == FACET->next; FACETc = FACETc->next)
//						FACETc->update = 0;
//					FACET->next->next = FACETc->next;
				}
			}
		}

		for (FACET = DB.FACET; FACET != NULL; FACET = FACET->next) {
			if (FACET->update) {
				FACET->update = 0;
				if (FACET->parent && FACET->parent->update)
					FACET->parent->update = 0;
				VIn  = FACET->VIn;
				VOut = FACET->VOut;
if (max(VIn->level,VOut->level)-min(VIn->level,VOut->level) > 1) {
	printf("%d %d %d\n",FACET->indexg,VIn->indexg,VOut->indexg);
	printf("Error: More than 1-irregular VOLUMEs (update_FACET).\n"), exit(1);
}
				FACET->P      = max(VIn->P,VOut->P);
				FACET->curved = max(VIn->curved,VOut->curved);

				// Compute XYZ_fS, n_fS, and detJF_fS
				setup_FACET_XYZ(FACET);
				setup_normals(FACET);
			}
		}
//exit(1);

		// HCOARSE
		for (l = LevelsMax; l > 0; l--) {
		for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
		if (VOLUME->update && VOLUME->adapt_type == HCOARSE && l == VOLUME->level) {
		if (VOLUME == VOLUME->parent->child0) {
			VOLUMEp = VOLUME->parent;

			if (VOLUMEp->update) {
				coarse_update(VOLUMEp);
			}
		}}}}

		for (FACET = DB.FACET; FACET != NULL; FACET = FACET->next) {
//			printf("Fpar: %d %d %d %d %p\n",FACET->indexg,FACET->level,FACET->update,FACET->adapt_type,FACET->parent);
		}

		// Fix list head if necessary
/*
for (FACET = DB.FACET; FACET; FACET = FACET->next) {
	printf("FTestCOA: %d %d %d\n",FACET->indexg,FACET->update,FACET->adapt_type);
}
printf("\n\n\n");
output_to_paraview("Geomadapt");
*/
		FACET = DB.FACET;

		if (FACET->update) {
			adapt_type = FACET->adapt_type;
			if (adapt_type == HCOARSE) {
				// DB.FACET->parent can never be NULL.
				DB.FACET = FACET->parent;
//				for (FACETc = FACET; FACETc->next && FACETc->next->parent == DB.FACET; FACETc = FACETc->next)
				for (FACETc = FACET; FACETc->next->parent == DB.FACET; FACETc = FACETc->next)
					;
				while (FACETc->next && FACETc->next->adapt_type == HDELETE) {
					FACETtmp = FACETc->next->next;
					memory_destructor_F(FACETc->next);
					FACETc->next = FACETtmp;
				}
				DB.FACET->next = FACETc->next;
				FACETc->next = NULL;
			} else if (adapt_type == HDELETE) {
				printf("Error: Should not be entering HDELETE while updating FACET head.\n"), exit(1);
			}
		}

		// Fix remainder of list


		for (FACET = DB.FACET; FACET; FACET = FACET->next) {
//printf("%d %d\n",FACET->indexg,FACET->level);
			FACETnext = FACET->next;
			if (FACETnext && FACETnext->update) {
				adapt_type = FACETnext->adapt_type;

				if (adapt_type == HDELETE) {
					while (FACETnext) {
						if (FACETnext->adapt_type == HDELETE) {
							FACETtmp = FACETnext->next;
							memory_destructor_F(FACETnext);
							FACETnext = FACETtmp;
						} else {
							break;
						}
					}
					if (FACETnext)
						adapt_type = FACETnext->adapt_type;
				}

				if (adapt_type == HCOARSE) {
					FACET->next = FACETnext->parent;
					for (FACETc = FACETnext; FACETc->next && FACETc->next->parent == FACET->next; FACETc = FACETc->next)
						;
					while (FACETc->next && FACETc->next->adapt_type == HDELETE) {
//printf("%d ",FACETc->next->adapt_type);
						FACETtmp = FACETc->next->next;
						memory_destructor_F(FACETc->next);
						FACETc->next = FACETtmp;
//printf("%d %d\n",FACETc->next->adapt_type,FACETtmp->next->adapt_type);
					}

					FACET->next->next = FACETc->next;
					FACETc->next = NULL;
				} else {
					FACET->next = FACETnext;
				}
			}
		}

		for (FACET = DB.FACET; FACET; FACET = FACET->next) {
//			printf("End h: %d\n",FACET->indexg);
		}
//exit(1);
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
