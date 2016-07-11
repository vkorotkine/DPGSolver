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
			*Vfh = 1;
			break;
		case 2:
			switch (fh) {
			default: // fh == 0
				*IndVInh = 0; break;
			case 1:
				*IndVInh = 1; break;
			*Vfh = 2;
			break;
		}
		break;
	default:
		printf("Error: Unsupported VType in get_FACET_IndVIn.\n");
		break;
	}
}

static unsigned int get_VfOut(const unsigned int neigh_f, const unsigned int VType)
{
	switch (VType) {
	case TRI:
		if (IndOrdInOut == 0)
			return neigh_f*NFREFMAX+fh;
		else
			return neigh_f*NFREFMAX+((fh+1)%2);

		break;
	default:
		printf("Error: Unsupported VType in get_VfOut.\n"), exit(1);
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

static void is_VOLUME_VIn(const unsigned int indexgVOLUME, const unsigned int indexgVIn)
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
	unsigned int Adapt = DB.Adapt,
	             NGF   = DB.NGF;

	// Standard datatypes
	unsigned int dummy_ui;

	struct S_VOLUME    *VOLUME, *VIn, *VOut, *dummyPtr_V;
	struct S_FACET     *FACET, *FACETh;

	switch (Adapt) {
	default: // ADAPT_HP
		break;
	case ADAPT_H:
// Different for HREFINE and HCOARSE
		for (l = 0; l < LevelsMax; l++) {
			for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
				if (VOLUME->update && VOLUME->adapt_type == HREFINE && l = VOLUME->level) {
					ELEMENT = get_ELEMENT_type(VOLUME->type);
					Nf = ELEMENT->Nf;

					// External FACETs
					for (f = 0; f < Nf; f++) {
						fInd = f*NSUBFMAX;
						sfMax = NsubF[f];
						if (sfMax == 1) { // Create new FACETs
							FACET = VOLUME->FACET[fInd];
							VIn   = FACET->VIn;
							VOut  = FACET->VOut;

							fhMax = get_fhMax(VOLUME->hrefine_type);
							for (fh = 0; fh < fhMax; fh++) {
								if (fh) {
									FACETh->next = New_FACET();
									FACETh = FACETh->next;
								} else {
									FACETh = New_FACET();
									FACET->child0 = FACETh;
								}
								FACETh->update = 1;

								FACETh->indexg = NGF++;
								FACETh->level  = (FACET->level)+1;
								FACETh->BC     = FACET->BC;

								// Find out if VOLUME == VIn or VOut
								if (is_VOLUME_VIn(VOLUME->indexg,VIn->indexg)) {
// If condition can be outside of fh loop
									get_FACET_IndVIn(FACET->VfIn,fh,VIn->type,&IndVInh,&Vfh);

									VOLUMEh = VOLUME->child0;
									for (vh = 0; vh < IndVInh; vh++)
										VOLUMEh = VOLUMEh->next;

									FACETh->VIn  = VOLUMEh;
									FACETh->VfIn = f*NFREFMAX;

									FACETh->VOut  = VOut;
									FACETh->VfOut = get_FACET_VfOut(VIn->neigh_f[f*NFREFMAX],VIn->type);

									// Valid for TET/PYR? (ToBeDeleted)
									FACETh->IndOrdInOut = FACET->IndOrdInOut;
									FACETh->IndOrdOutIn = FACET->IndOrdOutIn;



								} else {
									// DO NOT SWAP VIn and VOut

								}
							}
						} else { // Connect to existing FACETs
							for (sf = 0; sh < sfMax; sf++) {
								FACETh = VOLUME->FACET[fInd+sf];
							}

						}
					}


					// Internal FACETs (Always created)

				}
			}
		}
		DB.NGF = NGF;

		for (FACET = DB.FACET; FACET != NULL; FACET = FACET->next) {
			if (FACET->update) {
				FACET->update = 0;
				FACET->P = max(FACET->VIn->P,FACET->VOut->P);
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
