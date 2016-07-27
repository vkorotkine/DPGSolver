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

static void get_FACET_IndVIn(const unsigned int Vf, const unsigned int fh, const unsigned int VType,
                             unsigned int *IndVInh, unsigned int *Vfh)
{
	/*
	 *	Comments:
	 *		Vfh != f for PYR refinement. (ToBeDeleted)
	 */

	// Standard datatypes
	unsigned int f;

	f = Vf/NFREFMAX;

	switch (VType) {
	case TRI:
		switch (f) {
		default: // FACE 0
			switch (fh) {
			case 0: *IndVInh = 1; break;
			case 1: *IndVInh = 2; break;
			default: printf("Error: Unsupported (%d, %d, %d) in get_FACET_IndVIn.\n",VType,f,fh), exit(1); break;
			}
			*Vfh = 0;
			break;
		case 1:
			switch (fh) {
			case 0: *IndVInh = 0; break;
			case 1: *IndVInh = 2; break;
			default: printf("Error: Unsupported (%d, %d, %d) in get_FACET_IndVIn.\n",VType,f,fh), exit(1); break;
			}
			*Vfh = 1;
			break;
		case 2:
			switch (fh) {
			case 0: *IndVInh = 0; break;
			case 1: *IndVInh = 1; break;
			default: printf("Error: Unsupported (%d, %d, %d) in get_FACET_IndVIn.\n",VType,f,fh), exit(1); break;
			}
			*Vfh = 2;
			break;
		}
		break;
	case QUAD:
		switch (f) {
		default: // FACE 0
			switch (fh) {
			case 0: *IndVInh = 0; break;
			case 1: *IndVInh = 2; break;
			default: printf("Error: Unsupported (%d, %d, %d) in get_FACET_IndVIn.\n",VType,f,fh), exit(1); break;
			}
			*Vfh = 0;
			break;
		case 1:
			switch (fh) {
			case 0: *IndVInh = 1; break;
			case 1: *IndVInh = 3; break;
			default: printf("Error: Unsupported (%d, %d, %d) in get_FACET_IndVIn.\n",VType,f,fh), exit(1); break;
			}
			*Vfh = 1;
			break;
		case 2:
			switch (fh) {
			case 0: *IndVInh = 0; break;
			case 1: *IndVInh = 1; break;
			default: printf("Error: Unsupported (%d, %d, %d) in get_FACET_IndVIn.\n",VType,f,fh), exit(1); break;
			}
			*Vfh = 2;
			break;
		case 3:
			switch (fh) {
			case 0: *IndVInh = 2; break;
			case 1: *IndVInh = 3; break;
			default: printf("Error: Unsupported (%d, %d, %d) in get_FACET_IndVIn.\n",VType,f,fh), exit(1); break;
			}
			*Vfh = 3;
			break;
		}
		break;
	case TET:
		switch (f) {
		default: // FACE 0
			switch (fh) {
			case 0: *IndVInh = 1; break;
			case 1: *IndVInh = 2; break;
			case 2: *IndVInh = 3; break;
			case 3: *IndVInh = 4; break;
			default: printf("Error: Unsupported (%d, %d, %d) in get_FACET_IndVIn.\n",VType,f,fh), exit(1); break;
			}
			*Vfh = 0;
			break;
		case 1:
			switch (fh) {
			case 0: *IndVInh = 0; break;
			case 1: *IndVInh = 2; break;
			case 2: *IndVInh = 3; break;
			case 3: *IndVInh = 5; break;
			default: printf("Error: Unsupported (%d, %d, %d) in get_FACET_IndVIn.\n",VType,f,fh), exit(1); break;
			}
			*Vfh = 1;
			break;
		case 2:
			switch (fh) {
			case 0: *IndVInh = 0; break;
			case 1: *IndVInh = 1; break;
			case 2: *IndVInh = 3; break;
			case 3: *IndVInh = 6; break;
			default: printf("Error: Unsupported (%d, %d, %d) in get_FACET_IndVIn.\n",VType,f,fh), exit(1); break;
			}
			*Vfh = 2;
			break;
		case 3:
			switch (fh) {
			case 0: *IndVInh = 0; break;
			case 1: *IndVInh = 1; break;
			case 2: *IndVInh = 2; break;
			case 3: *IndVInh = 7; break;
			default: printf("Error: Unsupported (%d, %d, %d) in get_FACET_IndVIn.\n",VType,f,fh), exit(1); break;
			}
			*Vfh = 3;
			break;
		}
		break;
	case HEX:
		switch (f) {
		default: // FACE 0
			switch (fh) {
			case 0: *IndVInh = 0; break;
			case 1: *IndVInh = 2; break;
			case 2: *IndVInh = 4; break;
			case 3: *IndVInh = 6; break;
			default: printf("Error: Unsupported (%d, %d, %d) in get_FACET_IndVIn.\n",VType,f,fh), exit(1); break;
			}
			*Vfh = 0;
			break;
		case 1:
			switch (fh) {
			case 0: *IndVInh = 1; break;
			case 1: *IndVInh = 3; break;
			case 2: *IndVInh = 5; break;
			case 3: *IndVInh = 7; break;
			default: printf("Error: Unsupported (%d, %d, %d) in get_FACET_IndVIn.\n",VType,f,fh), exit(1); break;
			}
			*Vfh = 1;
			break;
		case 2:
			switch (fh) {
			case 0: *IndVInh = 0; break;
			case 1: *IndVInh = 1; break;
			case 2: *IndVInh = 4; break;
			case 3: *IndVInh = 5; break;
			default: printf("Error: Unsupported (%d, %d, %d) in get_FACET_IndVIn.\n",VType,f,fh), exit(1); break;
			}
			*Vfh = 2;
			break;
		case 3:
			switch (fh) {
			case 0: *IndVInh = 2; break;
			case 1: *IndVInh = 3; break;
			case 2: *IndVInh = 6; break;
			case 3: *IndVInh = 7; break;
			default: printf("Error: Unsupported (%d, %d, %d) in get_FACET_IndVIn.\n",VType,f,fh), exit(1); break;
			}
			*Vfh = 3;
			break;
		case 4:
			switch (fh) {
			case 0: *IndVInh = 0; break;
			case 1: *IndVInh = 1; break;
			case 2: *IndVInh = 2; break;
			case 3: *IndVInh = 3; break;
			default: printf("Error: Unsupported (%d, %d, %d) in get_FACET_IndVIn.\n",VType,f,fh), exit(1); break;
			}
			*Vfh = 4;
			break;
		case 5:
			switch (fh) {
			case 0: *IndVInh = 4; break;
			case 1: *IndVInh = 5; break;
			case 2: *IndVInh = 6; break;
			case 3: *IndVInh = 7; break;
			default: printf("Error: Unsupported (%d, %d, %d) in get_FACET_IndVIn.\n",VType,f,fh), exit(1); break;
			}
			*Vfh = 5;
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
	unsigned int Vfl;

	// silence
	Vfl = 0;

	switch (VType) {
	case TRI:
	case QUAD:
		// Isotropic refinement only
		switch (IndOrd) {
		default: // case 0
			Vfl = fh+1;
			break;
		case 1:
			if      (fh == 0) Vfl = 2;
			else if (fh == 1) Vfl = 1;
			else
				printf("Error: Unsupported fh for VType %d, IndOrd %d in get_FACET_VfOut.\n",VType,IndOrd), exit(1);
			break;
		}
		break;
	case TET:
		// Isotropic refinement only
		switch (IndOrd) {
		default: // case 0
			Vfl = fh+1;
			break;
		case 1:
			if      (fh == 0) Vfl = 2;
			else if (fh == 1) Vfl = 3;
			else if (fh == 2) Vfl = 1;
			else if (fh == 3) Vfl = fh+1;
			break;
		case 2:
			if      (fh == 0) Vfl = 3;
			else if (fh == 1) Vfl = 1;
			else if (fh == 2) Vfl = 2;
			else if (fh == 3) Vfl = fh+1;
			break;
		case 3:
			if      (fh == 0) Vfl = 1;
			else if (fh == 1) Vfl = 3;
			else if (fh == 2) Vfl = 2;
			else if (fh == 3) Vfl = fh+1;
			break;
		case 4:
			if      (fh == 0) Vfl = 3;
			else if (fh == 1) Vfl = 2;
			else if (fh == 2) Vfl = 1;
			else if (fh == 3) Vfl = fh+1;
			break;
		case 5:
			if      (fh == 0) Vfl = 2;
			else if (fh == 1) Vfl = 1;
			else if (fh == 2) Vfl = 3;
			else if (fh == 3) Vfl = fh+1;
			break;
		}
		break;
	case HEX:
		// Isotropic refinement only
		switch (IndOrd) {
		default: // case 0
			Vfl = fh+1;
			break;
		case 1:
			if      (fh == 0) Vfl = 1;
			else if (fh == 1) Vfl = 0;
			else if (fh == 2) Vfl = 3;
			else if (fh == 3) Vfl = 2;
			break;
		case 2:
			if      (fh == 0) Vfl = 2;
			else if (fh == 1) Vfl = 3;
			else if (fh == 2) Vfl = 0;
			else if (fh == 3) Vfl = 1;
			break;
		case 3:
			if      (fh == 0) Vfl = 3;
			else if (fh == 1) Vfl = 2;
			else if (fh == 2) Vfl = 1;
			else if (fh == 3) Vfl = 0;
			break;
		case 4:
			if      (fh == 0) Vfl = 0;
			else if (fh == 1) Vfl = 2;
			else if (fh == 2) Vfl = 1;
			else if (fh == 3) Vfl = 3;
			break;
		case 5:
			if      (fh == 0) Vfl = 2;
			else if (fh == 1) Vfl = 0;
			else if (fh == 2) Vfl = 3;
			else if (fh == 3) Vfl = 1;
			break;
		case 6:
			if      (fh == 0) Vfl = 1;
			else if (fh == 1) Vfl = 3;
			else if (fh == 2) Vfl = 0;
			else if (fh == 3) Vfl = 2;
			break;
		case 7:
			if      (fh == 0) Vfl = 3;
			else if (fh == 1) Vfl = 1;
			else if (fh == 2) Vfl = 2;
			else if (fh == 3) Vfl = 0;
			break;
		}
		break;
	default:
		printf("Error: Unsupported VType in get_VfOut.\n"), exit(1);
		break;
	}
	return neigh_f*NFREFMAX+Vfl;
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
	case QUAD:
		// Supported href_type: 0 (Isotropic)
		return 2;
		break;
	case TET:
		// Supported href_type: 0 (Isotropic)
		return 4;
		break;
	case HEX:
		// Supported href_type: 0 (Isotropic)
		return 4;
		break;
	default:
		printf("Error: Unsupported VType in get_fhMax.\n"), exit(1);
		return 0;
		break;
	}
}

static void get_Indsf(struct S_FACET *FACET, unsigned int *sfIn, unsigned int *sfOut);

static void set_FACET_Out(const unsigned int vh, const unsigned int fIn, struct S_FACET *FACETc, struct S_VOLUME *VOLUME)
{
	unsigned int i, VType, IndVhOut, f, IndOrdInOut, IndOrdOutIn, sfIn, sfOut;

	struct S_VOLUME *VOLUMEc;

	// silence
	IndVhOut = IndOrdInOut = IndOrdOutIn = 0;

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
	case QUAD:
		// Isotropic refinement only.
		if (vh == 0) {
			if      (fIn == 1) { IndVhOut = 1; f = 0; }
			else if (fIn == 3) { IndVhOut = 2; f = 2; }
			else               printf("Error: Unsupported (%d %d %d) in set_FACET_Out.\n",VType,vh,fIn), exit(1);
		} else if (vh == 1) {
			IndVhOut = 3; f = 2;
		} else if (vh == 2) {
			IndVhOut = 3; f = 0;
		} else { // Should already have found all FACETs
			printf("Error: Should not be entering set_FACET_Out for vh %d for VType %d.\n",vh,VType), exit(1);
		}
		IndOrdInOut = 0; // Same
		IndOrdOutIn = 0; // Same
		break;
	case TET:
		// Isotropic refinement only.
		switch (vh) {
		case 0: IndVhOut = 6; f = 0; IndOrdInOut = 4; IndOrdOutIn = 4; break; // fIn = 0
		case 1: IndVhOut = 7; f = 1; IndOrdInOut = 5; IndOrdOutIn = 5; break; // fIn = 1
		case 2: IndVhOut = 4; f = 2; IndOrdInOut = 3; IndOrdOutIn = 3; break; // fIn = 2
		case 3: IndVhOut = 5; f = 3; IndOrdInOut = 4; IndOrdOutIn = 4; break; // fIn = 3
		case 4:
			if      (fIn == 1) { IndVhOut = 7; IndOrdInOut = 2; IndOrdOutIn = 1; } // fIn = 1
			else if (fIn == 3) { IndVhOut = 5; IndOrdInOut = 1; IndOrdOutIn = 2; } // fIn = 3
			f = 2;
			break;
		case 5: IndVhOut = 6; f = 3; IndOrdInOut = 0; IndOrdOutIn = 0; break; // fIn = 0
		case 6: IndVhOut = 7; f = 0; IndOrdInOut = 1; IndOrdOutIn = 2; break; // fIn = 1
		default: // Should already have found all FACETs
			printf("Error: Should not be entering set_FACET_Out for vh %d for VType %d.\n",vh,VType), exit(1);
			break;
		}
		break;
	case HEX:
		// Isotropic refinement only.
		switch (vh) {
		case 0:
			if      (fIn == 1) { IndVhOut = 1; f = 0; }
			else if (fIn == 3) { IndVhOut = 2; f = 2; }
			else if (fIn == 5) { IndVhOut = 4; f = 4; }
			break;
		case 1:
			if      (fIn == 3) { IndVhOut = 3; f = 2; }
			else if (fIn == 5) { IndVhOut = 5; f = 4; }
			break;
		case 2:
			if      (fIn == 1) { IndVhOut = 3; f = 0; }
			else if (fIn == 5) { IndVhOut = 6; f = 4; }
			break;
		case 3: IndVhOut = 7; f = 4; break; // fIn = 5
		case 4:
			if      (fIn == 1) { IndVhOut = 5; f = 0; } // fIn = 1
			else if (fIn == 3) { IndVhOut = 6; f = 2; } // fIn = 3
			break;
		case 5: IndVhOut = 7; f = 2; break; // fIn = 3
		case 6: IndVhOut = 7; f = 0; break; // fIn = 1
		default: // Should already have found all FACETs
			printf("Error: Should not be entering set_FACET_Out for vh %d for VType %d.\n",vh,VType), exit(1);
			break;
		}
		IndOrdInOut = 0; // Same
		IndOrdOutIn = 0; // Same
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
	unsigned int i, VType, IndVhOut, VfOut, f, sfIn, sfOut;

	struct S_VOLUME *VOLUMEc;

	VType = VOLUME->type;
	VfOut = FACETc->VfOut;
	switch (VType) {
	case TRI:
		// Isotropic refinement only.
		switch (VfOut) {
		case NFREFMAX+1:
		case 2*NFREFMAX+1:
			IndVhOut = 0; break;
		case 1:
		case 2*NFREFMAX+2:
			IndVhOut = 1; break;
		case 2:
		case NFREFMAX+2:
			IndVhOut = 2; break;
		default:
//printf("%d %d\n",FACETc->VIn->indexg,FACETc->VOut->indexg);
			printf("Error: Unsupported VfOut = %d in set_FACET_Out_External.\n",VfOut), exit(1);
			break;
		}
		f = VfOut/NFREFMAX;
		break;
	case QUAD:
		// Isotropic refinement only.
		switch (VfOut) {
		case 1:
		case 2*NFREFMAX+1:
			IndVhOut = 0; break;
		case NFREFMAX+1:
		case 2*NFREFMAX+2:
			IndVhOut = 1; break;
		case 2:
		case 3*NFREFMAX+1:
			IndVhOut = 2; break;
		case NFREFMAX+2:
		case 3*NFREFMAX+2:
			IndVhOut = 3; break;
		default:
			printf("Error: Unsupported VfOut = %d in set_FACET_Out_External.\n",VfOut), exit(1);
			break;
		}
		f = VfOut/NFREFMAX;
		break;
	case TET:
		// Isotropic refinement only.
		switch (VfOut) {
		case NFREFMAX+1:
		case 2*NFREFMAX+1:
		case 3*NFREFMAX+1:
			IndVhOut = 0; break;
		case 1:
		case 2*NFREFMAX+2:
		case 3*NFREFMAX+2:
			IndVhOut = 1; break;
		case 2:
		case NFREFMAX+2:
		case 3*NFREFMAX+3:
			IndVhOut = 2; break;
		case 3:
		case NFREFMAX+3:
		case 2*NFREFMAX+3:
			IndVhOut = 3; break;
		case 4:
			IndVhOut = 4; break;
		case NFREFMAX+4:
			IndVhOut = 5; break;
		case 2*NFREFMAX+4:
			IndVhOut = 6; break;
		case 3*NFREFMAX+4:
			IndVhOut = 7; break;
		default:
			printf("Error: Unsupported VfOut = %d in set_FACET_Out_External.\n",VfOut), exit(1);
			break;
		}
		f = VfOut/NFREFMAX;
		break;
	case HEX:
		// Isotropic refinement only.
		switch (VfOut) {
		case 1:
		case 2*NFREFMAX+1:
		case 4*NFREFMAX+1:
			IndVhOut = 0; break;
		case NFREFMAX+1:
		case 2*NFREFMAX+2:
		case 4*NFREFMAX+2:
			IndVhOut = 1; break;
		case 2:
		case 3*NFREFMAX+1:
		case 4*NFREFMAX+3:
			IndVhOut = 2; break;
		case NFREFMAX+2:
		case 3*NFREFMAX+2:
		case 4*NFREFMAX+4:
			IndVhOut = 3; break;
		case 3:
		case 2*NFREFMAX+3:
		case 5*NFREFMAX+1:
			IndVhOut = 4; break;
		case NFREFMAX+3:
		case 2*NFREFMAX+4:
		case 5*NFREFMAX+2:
			IndVhOut = 5; break;
		case 4:
		case 3*NFREFMAX+3:
		case 5*NFREFMAX+3:
			IndVhOut = 6; break;
		case NFREFMAX+4:
		case 3*NFREFMAX+4:
		case 5*NFREFMAX+4:
			IndVhOut = 7; break;
		default:
			printf("Error: Unsupported VfOut = %d in set_FACET_Out_External.\n",VfOut), exit(1);
			break;
		}
		f = VfOut/NFREFMAX;
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
	unsigned int i, VType, Vf, f, Vfl, sf[2];

	for (i = 0; i < 2; i++) {
		if (i == 0) {
			Vf    = FACET->VfIn;
			VType = FACET->VIn->type;
		} else {
			Vf    = FACET->VfOut;
			VType = FACET->VOut->type;
		}
		f     = Vf/NFREFMAX;
		Vfl   = Vf%NFREFMAX;

		switch (VType) {
		case TRI:
		case QUAD:
			// Isotropic refinement only.
			switch (Vfl) {
			case 0:
				sf[i] = 0;
				break;
			case 1:
			case 2:
				sf[i] = Vfl-1;
				break;
			default:
				printf("Error: Unsupported Vfl in case %d of get_Indsf (%d).\n",VType,i), exit(1);
				break;
			}
			break;
		case TET:
			// Isotropic refinement only.
			switch (Vfl) {
			case 0:
				sf[i] = 0;
				break;
			case 1:
			case 2:
			case 3:
			case 4:
				sf[i] = Vfl-1;
				break;
			default:
				printf("Error: Unsupported Vfl in case %d of get_Indsf (%d).\n",VType,i), exit(1);
				break;
			}
			break;
		case HEX:
			// Isotropic refinement only.
			switch (Vfl) {
			case 0:
				sf[i] = 0;
				break;
			case 1:
			case 2:
			case 3:
			case 4:
				sf[i] = Vfl-1;
				break;
			default:
				printf("Error: Unsupported Vfl in case %d of get_Indsf (%d).\n",VType,i), exit(1);
				break;
			}
			break;
		default:
			printf("Error: Unsupported VType in get_Indsf (%d).\n",i), exit(1);
			break;
		}
		sf[i] += f*NSUBFMAX;
	}
	*sfIn  = sf[0];
	*sfOut = sf[1];
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
		case QUAD:
			// Supported: Isotropic refinement
			sfMax = 2;
			switch (f) {
			default: // f = 0
				IndVc[0] = 0;          IndVc[1] = 2;
				Indsf[0] = 0*NSUBFMAX; Indsf[1] = 0*NSUBFMAX;
				break;
			case 1:
				IndVc[0] = 1;          IndVc[1] = 3;
				Indsf[0] = 1*NSUBFMAX; Indsf[1] = 1*NSUBFMAX;
				break;
			case 2:
				IndVc[0] = 0;          IndVc[1] = 1;
				Indsf[0] = 2*NSUBFMAX; Indsf[1] = 2*NSUBFMAX;
				break;
			case 3:
				IndVc[0] = 2;          IndVc[1] = 3;
				Indsf[0] = 3*NSUBFMAX; Indsf[1] = 3*NSUBFMAX;
				break;
			}
			break;
		case TET:
			// Supported: Isotropic refinement
			sfMax = 4;
			switch (f) {
			default: IndVc[0] = 1; IndVc[1] = 2; IndVc[2] = 3; IndVc[3] = 4; break;
			case 1:  IndVc[0] = 0; IndVc[1] = 2; IndVc[2] = 3; IndVc[3] = 5; break;
			case 2:  IndVc[0] = 0; IndVc[1] = 1; IndVc[2] = 3; IndVc[3] = 6; break;
			case 3:  IndVc[0] = 0; IndVc[1] = 1; IndVc[2] = 2; IndVc[3] = 7; break;
			}
			for (sf = 0; sf < sfMax; sf++)
				Indsf[sf] = f*NSUBFMAX;
			break;
		case HEX:
			// Supported: Isotropic refinement
			sfMax = 4;
			switch (f) {
			default: IndVc[0] = 0; IndVc[1] = 2; IndVc[2] = 4; IndVc[3] = 6; break;
			case 1:  IndVc[0] = 1; IndVc[1] = 3; IndVc[2] = 5; IndVc[3] = 7; break;
			case 2:  IndVc[0] = 0; IndVc[1] = 1; IndVc[2] = 4; IndVc[3] = 5; break;
			case 3:  IndVc[0] = 2; IndVc[1] = 3; IndVc[2] = 6; IndVc[3] = 7; break;
			case 4:  IndVc[0] = 0; IndVc[1] = 1; IndVc[2] = 2; IndVc[3] = 3; break;
			case 5:  IndVc[0] = 4; IndVc[1] = 5; IndVc[2] = 6; IndVc[3] = 7; break;
			}
			for (sf = 0; sf < sfMax; sf++)
				Indsf[sf] = f*NSUBFMAX;
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
	case QUAD:
		sfMax_i = 4;

		IndVc[0] = 0; IndVc[1] = 0; IndVc[2] = 1; IndVc[3] = 2;
		Indsf[0] = 1; Indsf[1] = 3; Indsf[2] = 3; Indsf[3] = 1;
		for (i = 0; i < sfMax_i; i++)
			Indsf[i] *= NSUBFMAX;
		break;
	case TET:
		sfMax_i = 8;

		IndVc[0] = 0; IndVc[1] = 1; IndVc[2] = 2; IndVc[3] = 3;
		IndVc[4] = 4; IndVc[5] = 4; IndVc[6] = 5; IndVc[7] = 6;

		Indsf[0] = 0; Indsf[1] = 1; Indsf[2] = 2; Indsf[3] = 3;
		Indsf[4] = 1; Indsf[5] = 3; Indsf[6] = 0; Indsf[7] = 1;
		for (i = 0; i < sfMax_i; i++)
			Indsf[i] *= NSUBFMAX;
		break;
	case HEX:
		sfMax_i = 12;

		IndVc[0] = 0; IndVc[1] = 0; IndVc[2]  = 0; IndVc[3]  = 1;
		IndVc[4] = 1; IndVc[5] = 2; IndVc[6]  = 2; IndVc[7]  = 3;
		IndVc[8] = 4; IndVc[9] = 4; IndVc[10] = 5; IndVc[11] = 6;

		Indsf[0] = 1; Indsf[1] = 3; Indsf[2]  = 5; Indsf[3]  = 3;
		Indsf[4] = 5; Indsf[5] = 1; Indsf[6]  = 5; Indsf[7]  = 5;
		Indsf[8] = 1; Indsf[9] = 3; Indsf[10] = 3; Indsf[11] = 1;
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
	             Indf, IndVInh, adapt_type, BC, internal_BC,
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
							VIn  = FACET->VIn;
							VOut = FACET->VOut;
							BC   = FACET->BC;

							internal_BC = (BC == 0 || (BC % BC_STEP_SC > 50));

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
								FACETc->BC     = BC;

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

									if (internal_BC) {
										FACETc->VOut  = VOut;
										FACETc->VfOut = get_FACET_VfOut(fh,FACET->IndOrdOutIn,
										                                VIn->neigh_f[f*NFREFMAX],VIn->type);
									} else {
										FACETc->VOut  = FACETc->VIn;
										FACETc->VfOut = FACETc->VfIn;
									}
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

									FACETc->IndOrdInOut = FACET->IndOrdOutIn;
									FACETc->IndOrdOutIn = FACET->IndOrdInOut;

								}
								FACETc->VIn->NsubF[(FACETc->VfIn)/NFREFMAX] = 1;
								if (internal_BC)
									FACETc->VOut->NsubF[(FACETc->VfOut)/NFREFMAX] = fhMax;
								VOLUMEc->neigh_f[FACETc->VfIn] = (FACETc->VfOut)/NFREFMAX;

if (FACETc->VOut->indexg == 3) {
//	printf("FVOut3: %d %d %d\n",fh,FACETc->VfOut,(FACETc->VfOut)/NFREFMAX);
}

//printf("%d %d %d %d %d %d\n",FACET->indexg,VIn->indexg,VOut->indexg,FACETc->VfIn,FACETc->VfOut,FACETc->BC);

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
					for (VOLUMEc = VOLUME->child0; VOLUMEc; VOLUMEc = VOLUMEc->next) {
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
//printf("upF sF: %d %d\n",vh,f);
								set_FACET_Out(vh,f,FACETc,VOLUME);
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

		for (FACET = DB.FACET; FACET; FACET = FACET->next) {
			if (FACET->update) {
				FACET->update = 0;
				if (FACET->parent && FACET->parent->update)
					FACET->parent->update = 0;
				VIn  = FACET->VIn;
				VOut = FACET->VOut;
if (max(VIn->level,VOut->level)-min(VIn->level,VOut->level) > 1) {
	printf("%d %d %d %d %d\n",FACET->indexg,VIn->indexg,VOut->indexg,VIn->level,VOut->level);
	if (VIn->parent)
		printf("VIn:  %d %d %d\n",VIn->parent->indexg,VIn->parent->Vadapt,VIn->parent->adapt_type);
	if (VOut->parent)
		printf("VOut: %d %d %d\n",VOut->parent->indexg,VOut->parent->Vadapt,VOut->parent->adapt_type);
printf("%d\n",FACET->BC);
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
