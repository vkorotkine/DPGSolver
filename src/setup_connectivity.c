// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "setup_connectivity.h"

#include <stdlib.h>
#include <stdio.h>

#include "petscsys.h"
 
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"

#include "array_sort.h"
#include "array_print.h"

/*
 *	Purpose:
 *		Build global connectivity arrays. The main algorithm is based on tiConnect3D from Hesthaven's nodal DG code.
 *
 *	Comments:
 *		FToVe is sorted such that each row holds the list of vertices for the corresponding facet in ascending order.
 *		This ordering allows for comparison between facets from different volumes such that they can be matched up in
 *		VToV.
 *		If this function is found to be slow after profiling, it should be parallelized. (ToBeDeleted)
 *
 *	Notation:
 *		Fs   : (F)acet (s)tart index
 *		Vs   : (V)olume (s)tart index
 *
 *		VToVe : (V)olume to (Ve)rtex correspondence
 *		VType : (V)olume type numbering
 *
 *		NBF  : (N)umber of (B)oundary (F)acets
 *
 *		IndicesGF : (G)lobal (F)acet indices
 *		IndicesBF : (B)oundary (F)acet indices
 *
 *		FToVe  : (F)acet to (Ve)rtices correspondence.
 *		         This potentially includes non-existent facets. In this case, all facet vertices are set equal to NVe,
 *		         which is not a possible value.
 *		FNve   : (N)umber of local (ve)rtices on each (F)acet.
 *		BFToVe : (B)oundary (F)acet to (Ve)rtices correspondence.
 *		BTags  : (B)oundary Tags
 *
 *		NGF    : (N)umber of (G)lobal (F)acets
 *		         This number is modified as the connectivity/periodicity is established, and non-existent facets are
 *		         identified and eliminated.
 *		NVC    : (N)umber of (V)OLUMEs which are (C)urved
 *		NGFC   : (N)umber of (G)lobal (F)acets which are (C)urved
 *		VToV   : (V)olume to (V)olume connectivity
 *		VToF   : (V)olume to local (F)acet connectivity
 *		VToGF  : (V)olume to (G)lobal (F)acet correspondence
 *		VToBC  : (V)olume to (B)oundary (C)ondition flags
 *		GFToVe : (G)lobal (F)acet to (Ve)rtex correspondence
 *		VC     : List of (V)OLUMEs which are (C)urved
 *		GFC    : List of (G)lobal (F)acets which are (C)urved
 *
 *	References:
 *		Hesthaven (Nodal DG Code): https://github.com/tcew/nodal-dg
 *
 */

void setup_connectivity(void)
{
	// Initialize DB Parameters
	unsigned int d        = DB.d,
	             NVe      = DB.NVe,
	             NfMax    = DB.NfMax,
	             NfveMax  = DB.NfveMax,
	             *NE      = DB.NE,
	             *EToVe   = DB.EToVe,
	             *EType   = DB.EType,
	             *ETags   = DB.ETags;

	unsigned int PrintTesting = 0;

	// Standard datatypes
	unsigned int i, j, count, iMax,
	             IndB, IndM, iIn, iOut, NRows, NCols,
	             Fs, Vs, NV, v, f, gf, ve, vc, gfc, IndFixed[3],
	             *VeFcon, Nf, *Nfve,
	             NBF, NGF, *IndicesGF, *IndicesBF,
	             *VToVe, *VType, *FToVe, *FNve, *VToV, *VToF, *VToBC,
	             *BFToVe, *BTags, *VToGF, *GFToVeOver, *GFToVe, *NveGF,
	             *FNveSwap, *VToVSwap, *VToFSwap,
	             *IndicesMatchingOver, fNve[2], Match, *IndicesMatching, *matchIn, *matchOut,
	             vNeigh, fNeigh, Nve,
	             *VCOver, *GFCOver, *VC, *GFC;

	struct S_ELEMENT *ELEMENT;

	NV = NE[d];

	Fs = 0; for (i = 0; i < d-1; i++) Fs += NE[i];
	Vs = 0; for (i = 0; i < d; i++)   Vs += NE[i];

	VToVe = EToVe + Vs*8;
	VType = EType + Vs;

	NGF = NV*NfMax;

	FToVe = malloc(NGF*NfveMax * sizeof *FToVe); // free
	FNve  = malloc(NGF         * sizeof *FNve); // free

	// Initialize all entries of FToVe to NVe (Needed below for sorting)
	for (i = 0; i < NGF; i++)         FNve[i]  = 0;
	for (i = 0; i < NGF*NfveMax; i++) FToVe[i] = NVe;

	for (v = 0; v < NV; v++) {
		ELEMENT = DB.ELEMENT; while (ELEMENT->type != VType[v]) ELEMENT = ELEMENT->next;
		VeFcon = ELEMENT->VeFcon;
		Nf      = ELEMENT->Nf;
		Nfve    = ELEMENT->Nfve;

		IndFixed[0] = v*NfMax*NfveMax;
		for (f = 0; f < Nf; f++) {
			FNve[v*NfMax+f] = Nfve[f];

			IndFixed[1] = IndFixed[0]+f*NfveMax;

			for (ve = 0; ve < Nfve[f]; ve++)
				FToVe[IndFixed[1]+ve] = VToVe[v*8+VeFcon[f*4+ve]];

			PetscSortInt(Nfve[f],(int *)&FToVe[IndFixed[1]]);
		}
	}

	IndicesGF = malloc(NGF * sizeof *IndicesGF); // free
	for (i = 0; i < NGF; i++)
		IndicesGF[i] = i;

	array_sort_ui(NGF,NfveMax,FToVe,IndicesGF,'R','T');

	FNveSwap = malloc(NGF * sizeof *FNveSwap); // free
	VToVSwap = malloc(NGF * sizeof *VToVSwap); // free
	VToFSwap = malloc(NGF * sizeof *VToFSwap); // free

	for (v = 0; v < NV; v++) {
	for (f = 0; f < NfMax; f++) {
		FNveSwap[v*NfMax+f] = FNve[v*NfMax+f];
		VToVSwap[v*NfMax+f] = v;
		VToFSwap[v*NfMax+f] = f;
	}}

	VToV = malloc(NGF * sizeof *VToV); // keep
	VToF = malloc(NGF * sizeof *VToF); // keep

	// Sort based on Global Facet indexing obtained above
	for (v = count = 0; v < NV; v++) {
	for (f = 0; f < NfMax; f++) {
		FNve[v*NfMax+f] = FNveSwap[IndicesGF[count]];
		VToV[v*NfMax+f] = VToVSwap[IndicesGF[count]];
		VToF[v*NfMax+f] = VToFSwap[IndicesGF[count]];
		count++;
	}}

	free(FNveSwap);
	free(VToVSwap);
	free(VToFSwap);

	IndicesMatchingOver = malloc(NGF * sizeof *IndicesMatchingOver); // free

	for (i = 1, IndM = 0; i < NGF; i++) {
		fNve[0] = FNve[i-1];
		fNve[1] = FNve[i];

		if (FToVe[(i-1)*NfveMax] != NVe && (fNve[0] == fNve[1])) {
			for (j = 0, Match = 1; j < fNve[0]; j++) {
				if (FToVe[(i-1)*NfveMax+j] != FToVe[i*NfveMax+j]) {
					Match = 0;
					break;
				}
			}

			if (Match) {
				IndicesMatchingOver[IndM] = i-1;
				IndM++;
			}
		}
	}

	IndicesMatching = malloc(IndM * sizeof *IndicesMatching); // free
	for (i = 0; i < IndM; i++)
		IndicesMatching[i] = IndicesMatchingOver[i];
	free(IndicesMatchingOver);

	// Make links reflexive
	matchIn  = malloc(3*2*IndM * sizeof *matchIn); // free
	matchOut = malloc(3*2*IndM * sizeof *matchOut); // free

	for (i = 0, iIn = 0, iOut = 1; i < 2; i++) {
		if (i == 1) {
			iIn = 1;
			iOut = 0;
		}
		for (j = 0; j < IndM; j++) {
			matchIn[i*IndM*3+j*3+0]  = IndicesGF[IndicesMatching[j]+iIn];
			matchIn[i*IndM*3+j*3+1]  = VToV[IndicesMatching[j]+iIn];
			matchIn[i*IndM*3+j*3+2]  = VToF[IndicesMatching[j]+iIn];

			matchOut[i*IndM*3+j*3+0] = IndicesGF[IndicesMatching[j]+iOut];
			matchOut[i*IndM*3+j*3+1] = VToV[IndicesMatching[j]+iOut];
			matchOut[i*IndM*3+j*3+2] = VToF[IndicesMatching[j]+iOut];
		}
	}
	free(IndicesMatching);

	// Compute VToV and VToF
	for (v = 0; v < NV; v++) {
		ELEMENT = DB.ELEMENT; while (ELEMENT->type != VType[v]) ELEMENT = ELEMENT->next;
		Nf   = ELEMENT->Nf;

		for (f = 0; f < NfMax; f++) {
			if (f < Nf) {
				VToV[v*NfMax+f] = v;
				VToF[v*NfMax+f] = f;
			} else {
				VToV[v*NfMax+f] = NE[d];
				VToF[v*NfMax+f] = 6;
			}
		}
	}

	for (i = 0, iMax = 2*IndM; i < iMax; i++) {
		VToV[matchIn[i*3+0]] = matchOut[i*3+1];
		VToF[matchIn[i*3+0]] = matchOut[i*3+2];
	}

	free(matchIn);
	free(matchOut);

	// Setup Boundary Condition Information and find Curved VOLUMEs/FACETs
	NBF = NE[d-1];

	VToBC = malloc(NGF * sizeof *VToBC); // keep
	for (i = 0; i < NGF; i++)
		VToBC[i] = 0;
	BFToVe = malloc(NBF*NfveMax * sizeof *BFToVe); // free

	for (i = 0; i < NBF; i++) {
		for (j = 0; j < NfveMax; j++)
			BFToVe[i*NfveMax+j] = EToVe[(Fs+i)*8+j];
		PetscSortInt(NfveMax,(int *)&BFToVe[i*NfveMax]);
	}

	IndicesBF = malloc(NBF * sizeof *IndicesBF); // free
	for (i = 0; i < NBF; i++)
		IndicesBF[i] = i;

	array_sort_ui(NBF,NfveMax,BFToVe,IndicesBF,'R','T');

	BTags = malloc(NE[d-1]*(NfveMax+1) * sizeof *BTags); // free

	NRows = NE[d-1];
	NCols = NfveMax+1;
	for (i = 0; i < NRows; i++) {
	for (j = 0; j < NCols; j++) {
		if (j < NfveMax) BTags[i*NCols+j] = BFToVe[i*NfveMax+j];
		else             BTags[i*NCols+j] = ETags[(Fs+IndicesBF[i])*2+0];
	}}
	free(IndicesBF);

	for (gf = 0, IndB = 0; gf < NGF; gf++) {
		fNve[0] = FNve[gf];

		for (j = 0, Match = 1; j < fNve[0]; j++) {
			if ( BFToVe[IndB*NfveMax+j] == NVe ||
			    (BFToVe[IndB*NfveMax+j] != FToVe[gf*NfveMax+j])) {
					Match = 0;
					break;
			}
		}

		if (Match) {
			VToBC[VToV[IndicesGF[gf]]*NfMax+VToF[IndicesGF[gf]]] = BTags[IndB*(NfveMax+1)+NfveMax];
			IndB++;
		}
		if (IndB == NBF)
			break;
	}
	free(IndicesGF);

	// Setup VToGF and GFToVe (ToBeModified: Depending on what is needed for the element curving functions)
	VToGF      = malloc(NGF         * sizeof *VToGF); // keep
	GFToVeOver = malloc(NGF*NfveMax * sizeof *GFToVe); // free
	NveGF      = malloc(NGF         * sizeof *GFToVe); // free

	for (i = 0; i < NGF; i++)
		VToGF[i] = NGF;

	for (v = 0, gf = 0; v < NV; v++) {
	for (f = 0; f < NfMax; f++) {
		vNeigh = VToV[v*NfMax+f];
		fNeigh = VToF[v*NfMax+f];

		if (vNeigh == NV)
			continue;

		if (VToGF[v*NfMax+f] == NGF) {
			VToGF[v*NfMax+f] = gf;
			VToGF[vNeigh*NfMax+fNeigh] = gf;

			ELEMENT = DB.ELEMENT; while(ELEMENT->type != VType[v]) ELEMENT = ELEMENT->next;
			Nfve    = ELEMENT->Nfve;
			VeFcon = ELEMENT->VeFcon;

			Nve = Nfve[f];
			NveGF[gf] = Nve;

			for (ve = 0; ve < Nve; ve++)
				GFToVeOver[gf*NfveMax+ve] = VToVe[v*8+VeFcon[f*4+ve]];

			gf++;
		}
	}}
	NGF = gf;

	GFToVe = malloc(NGF*NfveMax * sizeof *GFToVe); // keep
	for (gf = 0; gf < NGF; gf++) {
		for (ve = 0; ve < NveGF[gf]; ve++)
			GFToVe[gf*NfveMax+ve] = GFToVeOver[gf*NfveMax+ve];
		for (ve = NveGF[gf]; ve < NfveMax; ve++)
			GFToVe[gf*NfveMax+ve] = NVe;
	}
	free(GFToVeOver);
	free(NveGF);

	// Store list of Volumes and Global Facets which are curved
	VCOver  = malloc(NV       * sizeof *VCOver); // free
	GFCOver = malloc(NV*NfMax * sizeof *GFCOver); // free
	for (v = 0, vc = 0, gfc = 0; v < NV; v++) {
	for (f = 0; f < NfMax; f++) {
		if (VToBC[v*NfMax+f] >= 20000) {
			if (vc == 0 || VCOver[vc-1] != v) {
				VCOver[vc] = v;
				vc++;
			}
			GFCOver[gfc] = VToGF[v*NfMax+f];
			gfc++;
		}
	}}

	VC  = malloc(vc  * sizeof *VC); // keep
	GFC = malloc(gfc * sizeof *GFC); // keep
	for (v = 0; v < vc; v++)     VC[v]   = VCOver[v];
	for (gf = 0; gf < gfc; gf++) GFC[gf] = GFCOver[gf];
	free(VCOver);
	free(GFCOver);

	// Assign DB Parameters
	DB.NV     = NV;
	DB.NGF    = NGF;
	DB.NVC    = vc;
	DB.NGFC   = gfc;
	DB.VToV   = VToV;
	DB.VToF   = VToF;
	DB.VToGF  = VToGF;
	DB.VToBC  = VToBC;
	DB.GFToVe = GFToVe;
	DB.VC     = VC;
	DB.GFC    = GFC;

	if (PrintTesting && !DB.MPIrank) {
		printf("VToVe:\n");          array_print_ui(DB.NV,8,VToVe,'R');
		//printf("VType:\n");          array_print_ui(DB.NV,1,VType,'R');
		printf("FNve (Sorted):\n");  array_print_ui(DB.NGF,1,FNve,'R');
		printf("FToVe (Sorted):\n"); array_print_ui(DB.NGF,DB.NfveMax,FToVe,'R');
		printf("VToV:\n");           array_print_ui(DB.NV,DB.NfMax,VToV,'R');
		printf("VToF:\n");           array_print_ui(DB.NV,DB.NfMax,VToF,'R');
		printf("BFToVe:\n");         array_print_ui(DB.NE[d-1],DB.NfveMax,BFToVe,'R');
		//printf("BTags:\n");          array_print_ui(DB.NE[d-1],DB.NfveMax+1,BTags,'R');
		printf("VToBC:\n");          array_print_ui(DB.NV,DB.NfMax,VToBC,'R');
		printf("VToGF:\n");          array_print_ui(DB.NV,DB.NfMax,VToGF,'R');
		printf("GFToVe:\n");         array_print_ui(DB.NGF,DB.NfveMax,GFToVe,'R');
		printf("VC:\n");             array_print_ui(1,DB.NVC,VC,'R');
		printf("GFC:\n");            array_print_ui(1,DB.NGFC,GFC,'R');
	}

	free(FToVe);
	free(FNve);
	free(BTags);
	free(BFToVe);
}
