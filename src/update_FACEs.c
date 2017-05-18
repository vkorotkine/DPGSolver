// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "update_FACEs.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "S_FACE.h"

#include "element_functions.h"
#include "adaptation.h"
#include "setup_geometry.h"
#include "setup_normals.h"
#include "memory_constructors.h"
#include "memory_destructors.h"

/*
 *	Purpose:
 *		Update FACE information/operators in ELEMENTs which have undergone hp refinement.
 *
 *	Comments:
 *		H-adaptation FACE updating is done from the coarsest to finest for HREFINE and in the opposite sense for
 *		HCOARSE in order to ensure that neighbouring elements are no more than 1-irregular.
 *		Clean up this function: function names, variable names, comments, potentially combine the COARSE and FINE FACE
 *		list updating (ToBeDeleted)
 *		coarse_update and get_FACE_IndVL use the same IndVLh => combine these during cleanup (ToBeDeleted)
 *		Likely include FACE renumbering as well to avoid overflow for long-time simulations. (ToBeDeleted)
 *
 *	Notation:
 *
 *	References:
 */

static void get_FACE_IndVL(const unsigned int Vf, const unsigned int fh, const unsigned int VType,
                             unsigned int *IndVLh, unsigned int *Vfh)
{
	/*
	 *	Comments:
	 *		Vfh != f for TET6 and PYR refinement.
	 */

	// Initialize DB Parameters
	unsigned int TETrefineType = DB.TETrefineType;

	// Standard datatypes
	unsigned int f;

	f = Vf/NFREFMAX;

	switch (VType) {
	case TRI:
		switch (f) {
		default: // FACE 0
			switch (fh) {
			case 0: *IndVLh = 1; break;
			case 1: *IndVLh = 2; break;
			default: printf("Error: Unsupported (%d, %d, %d).\n",VType,f,fh), EXIT_MSG; break;
			}
			*Vfh = 0;
			break;
		case 1:
			switch (fh) {
			case 0: *IndVLh = 0; break;
			case 1: *IndVLh = 2; break;
			default: printf("Error: Unsupported (%d, %d, %d).\n",VType,f,fh), EXIT_MSG; break;
			}
			*Vfh = 1;
			break;
		case 2:
			switch (fh) {
			case 0: *IndVLh = 0; break;
			case 1: *IndVLh = 1; break;
			default: printf("Error: Unsupported (%d, %d, %d).\n",VType,f,fh), EXIT_MSG; break;
			}
			*Vfh = 2;
			break;
		}
		break;
	case QUAD:
		switch (f) {
		default: // FACE 0
			switch (fh) {
			case 0: *IndVLh = 0; break;
			case 1: *IndVLh = 2; break;
			default: printf("Error: Unsupported (%d, %d, %d).\n",VType,f,fh), EXIT_MSG; break;
			}
			*Vfh = 0;
			break;
		case 1:
			switch (fh) {
			case 0: *IndVLh = 1; break;
			case 1: *IndVLh = 3; break;
			default: printf("Error: Unsupported (%d, %d, %d).\n",VType,f,fh), EXIT_MSG; break;
			}
			*Vfh = 1;
			break;
		case 2:
			switch (fh) {
			case 0: *IndVLh = 0; break;
			case 1: *IndVLh = 1; break;
			default: printf("Error: Unsupported (%d, %d, %d).\n",VType,f,fh), EXIT_MSG; break;
			}
			*Vfh = 2;
			break;
		case 3:
			switch (fh) {
			case 0: *IndVLh = 2; break;
			case 1: *IndVLh = 3; break;
			default: printf("Error: Unsupported (%d, %d, %d).\n",VType,f,fh), EXIT_MSG; break;
			}
			*Vfh = 3;
			break;
		}
		break;
	case TET:
		if (TETrefineType == TET8 || TETrefineType == TET12) {
			switch (f) {
			default: // FACE 0
				switch (fh) {
				case 0: *IndVLh = 2; break;
				case 1: *IndVLh = 3; break;
				case 2: *IndVLh = 1; break;
				case 3: *IndVLh = 4; break;
				default: printf("Error: Unsupported (%d, %d, %d).\n",VType,f,fh), EXIT_MSG; break;
				}
				*Vfh = 0;
				break;
			case 1:
				switch (fh) {
				case 0: *IndVLh = 2; break;
				case 1: *IndVLh = 3; break;
				case 2: *IndVLh = 0; break;
				case 3: *IndVLh = 5; break;
				default: printf("Error: Unsupported (%d, %d, %d).\n",VType,f,fh), EXIT_MSG; break;
				}
				*Vfh = 1;
				break;
			case 2:
				switch (fh) {
				case 0: *IndVLh = 0; break;
				case 1: *IndVLh = 1; break;
				case 2: *IndVLh = 3; break;
				case 3: *IndVLh = 6; break;
				default: printf("Error: Unsupported (%d, %d, %d).\n",VType,f,fh), EXIT_MSG; break;
				}
				*Vfh = 2;
				break;
			case 3:
				switch (fh) {
				case 0: *IndVLh = 0; break;
				case 1: *IndVLh = 1; break;
				case 2: *IndVLh = 2; break;
				case 3: *IndVLh = 7; break;
				default: printf("Error: Unsupported (%d, %d, %d).\n",VType,f,fh), EXIT_MSG; break;
				}
				*Vfh = 3;
				break;
			}
		} else if (TETrefineType == TET6) {
			switch (f) {
			default: // FACE 0
				switch (fh) {
				case 0: *IndVLh = 2; *Vfh = 0; break;
				case 1: *IndVLh = 3; *Vfh = 0; break;
				case 2: *IndVLh = 1; *Vfh = 0; break;
				case 3: *IndVLh = 5; *Vfh = 1; break;
				default: printf("Error: Unsupported (%d, %d, %d).\n",VType,f,fh), EXIT_MSG; break;
				}
				break;
			case 1:
				switch (fh) {
				case 0: *IndVLh = 2; *Vfh = 1; break;
				case 1: *IndVLh = 3; *Vfh = 1; break;
				case 2: *IndVLh = 0; *Vfh = 1; break;
				case 3: *IndVLh = 5; *Vfh = 0; break;
				default: printf("Error: Unsupported (%d, %d, %d).\n",VType,f,fh), EXIT_MSG; break;
				}
				break;
			case 2:
				switch (fh) {
				case 0: *IndVLh = 0; *Vfh = 2; break;
				case 1: *IndVLh = 1; *Vfh = 2; break;
				case 2: *IndVLh = 3; *Vfh = 2; break;
				case 3: *IndVLh = 4; *Vfh = 2; break;
				default: printf("Error: Unsupported (%d, %d, %d).\n",VType,f,fh), EXIT_MSG; break;
				}
				break;
			case 3:
				switch (fh) {
				case 0: *IndVLh = 0; *Vfh = 3; break;
				case 1: *IndVLh = 1; *Vfh = 3; break;
				case 2: *IndVLh = 2; *Vfh = 3; break;
				case 3: *IndVLh = 4; *Vfh = 3; break;
				default: printf("Error: Unsupported (%d, %d, %d).\n",VType,f,fh), EXIT_MSG; break;
				}
				break;
			}
		} else {
			printf("Error: Unsupported.\n"), EXIT_MSG;
		}
		break;
	case HEX:
		switch (f) {
		default: // FACE 0
			switch (fh) {
			case 0: *IndVLh = 0; break;
			case 1: *IndVLh = 2; break;
			case 2: *IndVLh = 4; break;
			case 3: *IndVLh = 6; break;
			default: printf("Error: Unsupported (%d, %d, %d).\n",VType,f,fh), EXIT_MSG; break;
			}
			*Vfh = 0;
			break;
		case 1:
			switch (fh) {
			case 0: *IndVLh = 1; break;
			case 1: *IndVLh = 3; break;
			case 2: *IndVLh = 5; break;
			case 3: *IndVLh = 7; break;
			default: printf("Error: Unsupported (%d, %d, %d).\n",VType,f,fh), EXIT_MSG; break;
			}
			*Vfh = 1;
			break;
		case 2:
			switch (fh) {
			case 0: *IndVLh = 0; break;
			case 1: *IndVLh = 1; break;
			case 2: *IndVLh = 4; break;
			case 3: *IndVLh = 5; break;
			default: printf("Error: Unsupported (%d, %d, %d).\n",VType,f,fh), EXIT_MSG; break;
			}
			*Vfh = 2;
			break;
		case 3:
			switch (fh) {
			case 0: *IndVLh = 2; break;
			case 1: *IndVLh = 3; break;
			case 2: *IndVLh = 6; break;
			case 3: *IndVLh = 7; break;
			default: printf("Error: Unsupported (%d, %d, %d).\n",VType,f,fh), EXIT_MSG; break;
			}
			*Vfh = 3;
			break;
		case 4:
			switch (fh) {
			case 0: *IndVLh = 0; break;
			case 1: *IndVLh = 1; break;
			case 2: *IndVLh = 2; break;
			case 3: *IndVLh = 3; break;
			default: printf("Error: Unsupported (%d, %d, %d).\n",VType,f,fh), EXIT_MSG; break;
			}
			*Vfh = 4;
			break;
		case 5:
			switch (fh) {
			case 0: *IndVLh = 4; break;
			case 1: *IndVLh = 5; break;
			case 2: *IndVLh = 6; break;
			case 3: *IndVLh = 7; break;
			default: printf("Error: Unsupported (%d, %d, %d).\n",VType,f,fh), EXIT_MSG; break;
			}
			*Vfh = 5;
			break;
		}
		break;
	case WEDGE:
		switch (f) {
		default: // FACE 0
			switch (fh) {
			case 0: *IndVLh = 1; break;
			case 1: *IndVLh = 2; break;
			case 2: *IndVLh = 5; break;
			case 3: *IndVLh = 6; break;
			default: printf("Error: Unsupported (%d, %d, %d).\n",VType,f,fh), EXIT_MSG; break;
			}
			*Vfh = 0;
			break;
		case 1:
			switch (fh) {
			case 0: *IndVLh = 0; break;
			case 1: *IndVLh = 2; break;
			case 2: *IndVLh = 4; break;
			case 3: *IndVLh = 6; break;
			default: printf("Error: Unsupported (%d, %d, %d).\n",VType,f,fh), EXIT_MSG; break;
			}
			*Vfh = 1;
			break;
		case 2:
			switch (fh) {
			case 0: *IndVLh = 0; break;
			case 1: *IndVLh = 1; break;
			case 2: *IndVLh = 4; break;
			case 3: *IndVLh = 5; break;
			default: printf("Error: Unsupported (%d, %d, %d).\n",VType,f,fh), EXIT_MSG; break;
			}
			*Vfh = 2;
			break;
		case 3:
			switch (fh) {
			case 0: *IndVLh = 0; break;
			case 1: *IndVLh = 1; break;
			case 2: *IndVLh = 2; break;
			case 3: *IndVLh = 3; break;
			default: printf("Error: Unsupported (%d, %d, %d).\n",VType,f,fh), EXIT_MSG; break;
			}
			*Vfh = 3;
			break;
		case 4:
			switch (fh) {
			case 0: *IndVLh = 4; break;
			case 1: *IndVLh = 5; break;
			case 2: *IndVLh = 6; break;
			case 3: *IndVLh = 7; break;
			default: printf("Error: Unsupported (%d, %d, %d).\n",VType,f,fh), EXIT_MSG; break;
			}
			*Vfh = 4;
			break;
		}
		break;
	case PYR:
		switch (f) {
		default: // FACE 0
			switch (fh) {
			case 0: *IndVLh = 0; *Vfh = 0; break;
			case 1: *IndVLh = 2; *Vfh = 0; break;
			case 2: *IndVLh = 9; *Vfh = 0; break;
			case 3: *IndVLh = 4; *Vfh = 1; break;
			default: printf("Error: Unsupported (%d, %d, %d).\n",VType,f,fh), EXIT_MSG; break;
			}
			break;
		case 1:
			switch (fh) {
			case 0: *IndVLh = 1; *Vfh = 1; break;
			case 1: *IndVLh = 3; *Vfh = 1; break;
			case 2: *IndVLh = 9; *Vfh = 1; break;
			case 3: *IndVLh = 5; *Vfh = 0; break;
			default: printf("Error: Unsupported (%d, %d, %d).\n",VType,f,fh), EXIT_MSG; break;
			}
			break;
		case 2:
			switch (fh) {
			case 0: *IndVLh = 0; break;
			case 1: *IndVLh = 1; break;
			case 2: *IndVLh = 9; break;
			case 3: *IndVLh = 6; break;
			default: printf("Error: Unsupported (%d, %d, %d).\n",VType,f,fh), EXIT_MSG; break;
			}
			*Vfh = 2;
			break;
		case 3:
			switch (fh) {
			case 0: *IndVLh = 2; break;
			case 1: *IndVLh = 3; break;
			case 2: *IndVLh = 9; break;
			case 3: *IndVLh = 7; break;
			default: printf("Error: Unsupported (%d, %d, %d).\n",VType,f,fh), EXIT_MSG; break;
			}
			*Vfh = 3;
			break;
		case 4:
			switch (fh) {
			case 0: *IndVLh = 0; break;
			case 1: *IndVLh = 1; break;
			case 2: *IndVLh = 2; break;
			case 3: *IndVLh = 3; break;
			default: printf("Error: Unsupported (%d, %d, %d).\n",VType,f,fh), EXIT_MSG; break;
			}
			*Vfh = 4;
			break;
		}
		break;
	default:
		printf("Error: Unsupported VType.\n"), EXIT_MSG;
		break;
	}
}

static unsigned int get_FACE_VfR(const unsigned int fh, const unsigned int IndOrd, const unsigned int neigh_f,
                                    const unsigned int FType)
{
	/*
	 *	Purpose:
	 *		Returns VfR based on IndOrd of neighbouring parent elements.
	 *
	 *	Comments:
	 *		In the current implementation, the IndOrd of the refined external faces is the same as that of the parent
	 *		faces. If TET6 refinement algorithms with internal PYRs in arbitrary orientations are desired, this would
	 *		no longer be true and IndOrd would also need to be returned.
	 */

	unsigned int Vfl;

	// silence
	Vfl = 0;

	switch (FType) {
	case LINE:
		// Isotropic refinement only
		switch (IndOrd) {
		default: // case 0
			Vfl = fh+1;
			break;
		case 1:
			if      (fh == 0) Vfl = 2;
			else if (fh == 1) Vfl = 1;
			else
				printf("Error: Unsupported fh for FType %d, IndOrd %d.\n",FType,IndOrd), EXIT_MSG;
			break;
		}
		break;
	case TRI:
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
	case QUAD:
		// Isotropic refinement only
		switch (IndOrd) {
		default: // case 0
			Vfl = fh+1;
			break;
		case 1:
			if      (fh == 0) Vfl = 2;
			else if (fh == 1) Vfl = 1;
			else if (fh == 2) Vfl = 4;
			else if (fh == 3) Vfl = 3;
			break;
		case 2:
			if      (fh == 0) Vfl = 3;
			else if (fh == 1) Vfl = 4;
			else if (fh == 2) Vfl = 1;
			else if (fh == 3) Vfl = 2;
			break;
		case 3:
			if      (fh == 0) Vfl = 4;
			else if (fh == 1) Vfl = 3;
			else if (fh == 2) Vfl = 2;
			else if (fh == 3) Vfl = 1;
			break;
		case 4:
			if      (fh == 0) Vfl = 1;
			else if (fh == 1) Vfl = 3;
			else if (fh == 2) Vfl = 2;
			else if (fh == 3) Vfl = 4;
			break;
		case 5:
			if      (fh == 0) Vfl = 3;
			else if (fh == 1) Vfl = 1;
			else if (fh == 2) Vfl = 4;
			else if (fh == 3) Vfl = 2;
			break;
		case 6:
			if      (fh == 0) Vfl = 2;
			else if (fh == 1) Vfl = 4;
			else if (fh == 2) Vfl = 1;
			else if (fh == 3) Vfl = 3;
			break;
		case 7:
			if      (fh == 0) Vfl = 4;
			else if (fh == 1) Vfl = 2;
			else if (fh == 2) Vfl = 3;
			else if (fh == 3) Vfl = 1;
			break;
		}
		break;
	default:
		printf("Error: Unsupported FType.\n"), EXIT_MSG;
		break;
	}
	return neigh_f*NFREFMAX+Vfl;
}

static unsigned int get_FACE_type(struct S_FACE *FACE)
{
	// Initialize DB Parameters
	const unsigned int d = DB.d;

	// Standard datatypesc
	unsigned int VType, VfL, fIn;

	struct S_VOLUME *VL = FACE->VL;

	VType = VL->type;
	VfL  = FACE->VfL;
	fIn   = VfL/NFREFMAX;

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
	if (href_type)
		printf("Error: Unsupported href_type (%d).\n",href_type), EXIT_MSG;

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
	case WEDGE:
		// Supported href_type: 0 (Isotropic)
		return 4;
		break;
	case PYR:
		// Supported href_type: 0 (Isotropic)
		return 4;
		break;
	default:
		printf("Error: Unsupported VType.\n"), EXIT_MSG;
		return 0;
		break;
	}
}

static void get_Indsf(struct S_FACE *FACE, unsigned int *sfIn, unsigned int *sfOut);

static void set_FACE_Out(const unsigned int vh, const unsigned int fIn, struct S_FACE *FACEc, struct S_VOLUME *VOLUME)
{
	// Initialize DB Parameters
	unsigned int TETrefineType = DB.TETrefineType;

	// Standard datatypes
	unsigned int i, VType, IndVhOut, f, IndOrdLR, IndOrdRL, sfIn, sfOut;

	struct S_VOLUME *VOLUMEc;

	// silence
	IndVhOut = IndOrdLR = IndOrdRL = f = 0;

	VType = VOLUME->type;
	switch (VType) {
	case TRI:
		// Isotropic refinement only.
		IndVhOut = 3;
		f = vh;
		IndOrdLR = 1; // Reversed
		IndOrdRL = 1; // Reversed
		if (vh == 3)
			printf("Error: Should not be entering for vh %d for VType %d.\n",vh,VType), EXIT_MSG;
			// Should already have found all FACEs
		break;
	case QUAD:
		// Isotropic refinement only.
		if (vh == 0) {
			if      (fIn == 1) { IndVhOut = 1; f = 0; }
			else if (fIn == 3) { IndVhOut = 2; f = 2; }
			else               printf("Error: Unsupported (%d %d %d).\n",VType,vh,fIn), EXIT_MSG;
		} else if (vh == 1) {
			IndVhOut = 3; f = 2;
		} else if (vh == 2) {
			IndVhOut = 3; f = 0;
		} else { // Should already have found all FACEs
			printf("Error: Should not be entering for vh %d for VType %d.\n",vh,VType), EXIT_MSG;
		}
		IndOrdLR = 0; // Same
		IndOrdRL = 0; // Same
		break;
	case TET:
		// Isotropic refinement only.
		if (TETrefineType == TET8) {
			switch (vh) {
			case 0: IndVhOut = 5; f = 0; IndOrdLR = 5; IndOrdRL = 5; break; // fIn = 0
			case 1: IndVhOut = 4; f = 1; IndOrdLR = 5; IndOrdRL = 5; break; // fIn = 1
			case 2: IndVhOut = 7; f = 2; IndOrdLR = 5; IndOrdRL = 5; break; // fIn = 2
			case 3: IndVhOut = 6; f = 3; IndOrdLR = 5; IndOrdRL = 5; break; // fIn = 3
			case 4:
				if      (fIn == 2) { IndVhOut = 7; IndOrdLR = 0; IndOrdRL = 0; }
				else if (fIn == 3) { IndVhOut = 6; IndOrdLR = 5; IndOrdRL = 5; }
				f = 1;
				break;
			case 5:
				if      (fIn == 2) { IndVhOut = 7; IndOrdLR = 5; IndOrdRL = 5; }
				else if (fIn == 3) { IndVhOut = 6; IndOrdLR = 0; IndOrdRL = 0; }
				f = 0;
				break;
			default: // Should already have found all FACEs
				printf("Error: Should not be entering for vh %d for VType %d.\n",vh,VType), EXIT_MSG;
				break;
			}
		} else if (TETrefineType == TET12) {
			switch (vh) {
			case 0: IndVhOut = 8;  f = 0; IndOrdLR = 5; IndOrdRL = 5; break; // fIn = 0
			case 1: IndVhOut = 9;  f = 1; IndOrdLR = 5; IndOrdRL = 5; break; // fIn = 1
			case 2: IndVhOut = 10; f = 2; IndOrdLR = 5; IndOrdRL = 5; break; // fIn = 2
			case 3: IndVhOut = 11; f = 3; IndOrdLR = 5; IndOrdRL = 5; break; // fIn = 3
			case 4:
				if      (fIn == 1) { IndVhOut = 9;  f = 0; IndOrdLR = 0; IndOrdRL = 0; }
				else if (fIn == 2) { IndVhOut = 10; f = 1; IndOrdLR = 0; IndOrdRL = 0; }
				else if (fIn == 3) { IndVhOut = 11; f = 1; IndOrdLR = 5; IndOrdRL = 5; }
				break;
			case 5:
				if      (fIn == 0) { IndVhOut = 8;  f = 1; IndOrdLR = 0; IndOrdRL = 0; }
				else if (fIn == 2) { IndVhOut = 10; f = 0; IndOrdLR = 5; IndOrdRL = 5; }
				else if (fIn == 3) { IndVhOut = 11; f = 0; IndOrdLR = 0; IndOrdRL = 0; }
				break;
			case 6:
				if      (fIn == 0) { IndVhOut = 8;  f = 3; IndOrdLR = 0; IndOrdRL = 0; }
				else if (fIn == 1) { IndVhOut = 9;  f = 3; IndOrdLR = 5; IndOrdRL = 5; }
				else if (fIn == 3) { IndVhOut = 11; f = 2; IndOrdLR = 0; IndOrdRL = 0; }
				break;
			case 7:
				if      (fIn == 0) { IndVhOut = 8;  f = 2; IndOrdLR = 5; IndOrdRL = 5; }
				else if (fIn == 1) { IndVhOut = 9;  f = 2; IndOrdLR = 0; IndOrdRL = 0; }
				else if (fIn == 2) { IndVhOut = 10; f = 3; IndOrdLR = 0; IndOrdRL = 0; }
				break;
			default: // Should already have found all FACEs
				printf("Error: Should not be entering for vh %d for VType %d.\n",vh,VType), EXIT_MSG;
				break;
			}
		} else if (TETrefineType == TET6) {
			switch (vh) {
			case 0: IndVhOut = 4; f = 1; IndOrdLR = 5; IndOrdRL = 5; break; // fIn = 0
			case 1: IndVhOut = 4; f = 0; IndOrdLR = 5; IndOrdRL = 5; break; // fIn = 1
			case 2: IndVhOut = 5; f = 3; IndOrdLR = 0; IndOrdRL = 0; break; // fIn = 2
			case 3: IndVhOut = 5; f = 2; IndOrdLR = 0; IndOrdRL = 0; break; // fIn = 3
			case 4: IndVhOut = 5; f = 4; IndOrdLR = 1; IndOrdRL = 1; break; // fIn = 4
			default: // Should already have found all FACEs
				printf("Error: Should not be entering for vh %d for VType %d.\n",vh,VType), EXIT_MSG;
				break;
			}
		} else {
			printf("Error: Unsupported.\n"), EXIT_MSG;
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
		default: // Should already have found all FACEs
			printf("Error: Should not be entering for vh %d for VType %d.\n",vh,VType), EXIT_MSG;
			break;
		}
		IndOrdLR = 0; // Same
		IndOrdRL = 0; // Same
		break;
	case WEDGE:
		// Isotropic refinement only.
		switch (vh) {
		case 0:
			if      (fIn == 0) { IndVhOut = 3; f = 0; IndOrdLR = 1; IndOrdRL = 1; }
			else if (fIn == 4) { IndVhOut = 4; f = 3; IndOrdLR = 0; IndOrdRL = 0; }
			break;
		case 1:
			if      (fIn == 1) { IndVhOut = 3; f = 1; IndOrdLR = 1; IndOrdRL = 1; }
			else if (fIn == 4) { IndVhOut = 5; f = 3; IndOrdLR = 0; IndOrdRL = 0; }
			break;
		case 2:
			if      (fIn == 2) { IndVhOut = 3; f = 2; IndOrdLR = 1; IndOrdRL = 1; }
			else if (fIn == 4) { IndVhOut = 6; f = 3; IndOrdLR = 0; IndOrdRL = 0; }
			break;
		case 3: IndVhOut = 7; f = 3; IndOrdLR = 0; IndOrdRL = 0; break; // fIn = 4
		case 4: IndVhOut = 7; f = 0; IndOrdLR = 1; IndOrdRL = 1; break; // fIn = 0
		case 5: IndVhOut = 7; f = 1; IndOrdLR = 1; IndOrdRL = 1; break; // fIn = 1
		case 6: IndVhOut = 7; f = 2; IndOrdLR = 1; IndOrdRL = 1; break; // fIn = 2
		default: // Should already have found all FACEs
			printf("Error: Should not be entering for vh %d for VType %d.\n",vh,VType), EXIT_MSG;
			break;
		}
		break;
	case PYR:
		// Isotropic refinement only.
		switch (vh) {
		case 0:
			if      (fIn == 1) { IndVhOut = 6; f = 0; IndOrdLR = 5; IndOrdRL = 5; }
			else if (fIn == 3) { IndVhOut = 4; f = 2; IndOrdLR = 0; IndOrdRL = 0; }
			break;
		case 1:
			if      (fIn == 0) { IndVhOut = 6; f = 1; IndOrdLR = 5; IndOrdRL = 5; }
			else if (fIn == 3) { IndVhOut = 5; f = 2; IndOrdLR = 0; IndOrdRL = 0; }
			break;
		case 2:
			if      (fIn == 1) { IndVhOut = 7; f = 0; IndOrdLR = 5; IndOrdRL = 5; }
			else if (fIn == 2) { IndVhOut = 4; f = 3; IndOrdLR = 0; IndOrdRL = 0; }
			break;
		case 3:
			if      (fIn == 0) { IndVhOut = 7; f = 1; IndOrdLR = 5; IndOrdRL = 5; }
			else if (fIn == 2) { IndVhOut = 5; f = 3; IndOrdLR = 0; IndOrdRL = 0; }
			break;
		case 4: IndVhOut = 8; f = 1; IndOrdLR = 5; IndOrdRL = 5; break; // fIn = 0
		case 5: IndVhOut = 8; f = 0; IndOrdLR = 5; IndOrdRL = 5; break; // fIn = 1
		case 6: IndVhOut = 8; f = 2; IndOrdLR = 0; IndOrdRL = 0; break; // fIn = 3
		case 7: IndVhOut = 8; f = 3; IndOrdLR = 0; IndOrdRL = 0; break; // fIn = 2
		case 8: IndVhOut = 9; f = 4; IndOrdLR = 1; IndOrdRL = 1; break; // fIn = 4
		default: // Should already have found all FACEs
			printf("Error: Should not be entering for vh %d for VType %d.\n",vh,VType), EXIT_MSG;
			break;
		}
		break;
	default:
		printf("Error: Unsupported VType.\n"), EXIT_MSG;
		break;
	}

	VOLUMEc = VOLUME->child0;
	for (i = 0; i < IndVhOut; i++)
		VOLUMEc = VOLUMEc->next;

	FACEc->VR  = VOLUMEc;
	FACEc->VfR = f*NFREFMAX;

	FACEc->IndOrdLR = IndOrdLR;
	FACEc->IndOrdRL = IndOrdRL;

	get_Indsf(FACEc,&sfIn,&sfOut);

	FACEc->VL->FACE[sfIn] = FACEc;
	FACEc->VR->FACE[sfOut] = FACEc;
	FACEc->type = get_FACE_type(FACEc);

	FACEc->VL->neigh_f[FACEc->VfL]   = (FACEc->VfR)/NFREFMAX;
	FACEc->VR->neigh_f[FACEc->VfR] = (FACEc->VfL)/NFREFMAX;
}

static void set_FACE_Out_External(struct S_FACE *FACEc, struct S_VOLUME *VOLUME)
{
	// Initialize DB Parameters
	unsigned int TETrefineType = DB.TETrefineType;

	// standard datatypes
	unsigned int i, VType, IndVhOut, VfR, f, sfIn, sfOut;

	struct S_VOLUME *VOLUMEc;

	VType = VOLUME->type;
	VfR = FACEc->VfR;
	switch (VType) {
	case TRI:
		// Isotropic refinement only.
		switch (VfR) {
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
			printf("Error: Unsupported VfR = %d.\n",VfR), EXIT_MSG;
			break;
		}
		f = VfR/NFREFMAX;
		break;
	case QUAD:
		// Isotropic refinement only.
		switch (VfR) {
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
			printf("Error: Unsupported VfR = %d.\n",VfR), EXIT_MSG;
			break;
		}
		f = VfR/NFREFMAX;
		break;
	case TET:
		// Isotropic refinement only.
		if (TETrefineType == TET8 || TETrefineType == TET12) {
			switch (VfR) {
			case NFREFMAX+3:
			case 2*NFREFMAX+1:
			case 3*NFREFMAX+1:
				IndVhOut = 0; break;
			case 3:
			case 2*NFREFMAX+2:
			case 3*NFREFMAX+2:
				IndVhOut = 1; break;
			case 1:
			case NFREFMAX+1:
			case 3*NFREFMAX+3:
				IndVhOut = 2; break;
			case 2:
			case NFREFMAX+2:
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
				printf("Error: Unsupported VfR = %d.\n",VfR), EXIT_MSG;
				break;
			}
			f = VfR/NFREFMAX;
		} else if (TETrefineType == TET6) {
			switch (VfR) {
			case NFREFMAX+3:
			case 2*NFREFMAX+1:
			case 3*NFREFMAX+1:
				IndVhOut = 0; f = VfR/NFREFMAX; break;
			case 3:
			case 2*NFREFMAX+2:
			case 3*NFREFMAX+2:
				IndVhOut = 1; f = VfR/NFREFMAX; break;
			case 1:
			case NFREFMAX+1:
			case 3*NFREFMAX+3:
				IndVhOut = 2; f = VfR/NFREFMAX; break;
			case 2:
			case NFREFMAX+2:
			case 2*NFREFMAX+3:
				IndVhOut = 3; f = VfR/NFREFMAX; break;
			case 4:
				IndVhOut = 5; f = 1; break;
			case NFREFMAX+4:
				IndVhOut = 5; f = 0; break;
			case 2*NFREFMAX+4:
				IndVhOut = 4; f = 2; break;
			case 3*NFREFMAX+4:
				IndVhOut = 4; f = 3; break;
			default:
				printf("Error: Unsupported VfR = %d.\n",VfR), EXIT_MSG;
				break;
			}
		} else {
			printf("Error: Unsupported.\n"), EXIT_MSG;
		}
		break;
	case HEX:
		// Isotropic refinement only.
		switch (VfR) {
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
			printf("Error: Unsupported VfR = %d.\n",VfR), EXIT_MSG;
			break;
		}
		f = VfR/NFREFMAX;
		break;
	case WEDGE:
		// Isotropic refinement only.
		switch (VfR) {
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
		case 3*NFREFMAX+4:
			IndVhOut = 3; break;
		case NFREFMAX+3:
		case 2*NFREFMAX+3:
		case 4*NFREFMAX+1:
			IndVhOut = 4; break;
		case 3:
		case 2*NFREFMAX+4:
		case 4*NFREFMAX+2:
			IndVhOut = 5; break;
		case 4:
		case NFREFMAX+4:
		case 4*NFREFMAX+3:
			IndVhOut = 6; break;
		case 4*NFREFMAX+4:
			IndVhOut = 7; break;
		default:
			printf("Error: Unsupported VfR = %d.\n",VfR), EXIT_MSG;
			break;
		}
		f = VfR/NFREFMAX;
		break;
	case PYR:
		// Isotropic refinement only.
		switch (VfR) {
		case 1:
		case 2*NFREFMAX+1:
		case 4*NFREFMAX+1:
			IndVhOut = 0; f = VfR/NFREFMAX; break;
		case NFREFMAX+1:
		case 2*NFREFMAX+2:
		case 4*NFREFMAX+2:
			IndVhOut = 1; f = VfR/NFREFMAX; break;
		case 2:
		case 3*NFREFMAX+1:
		case 4*NFREFMAX+3:
			IndVhOut = 2; f = VfR/NFREFMAX; break;
		case NFREFMAX+2:
		case 3*NFREFMAX+2:
		case 4*NFREFMAX+4:
			IndVhOut = 3; f = VfR/NFREFMAX; break;
		case 4:
			IndVhOut = 4; f = 1; break;
		case NFREFMAX+4:
			IndVhOut = 5; f = 0; break;
		case 2*NFREFMAX+4:
			IndVhOut = 6; f = VfR/NFREFMAX; break;
		case 3*NFREFMAX+4:
			IndVhOut = 7; f = VfR/NFREFMAX; break;
		case 3:
		case NFREFMAX+3:
		case 2*NFREFMAX+3:
		case 3*NFREFMAX+3:
			IndVhOut = 9; f = VfR/NFREFMAX; break;
		default:
			printf("Error: Unsupported VfR = %d.\n",VfR), EXIT_MSG;
			break;
		}
		break;
	default:
		printf("Error: Unsupported VType.\n"), EXIT_MSG;
		break;
	}

	VOLUMEc = VOLUME->child0;
	for (i = 0; i < IndVhOut; i++)
		VOLUMEc = VOLUMEc->next;

	FACEc->VR  = VOLUMEc;
	FACEc->VfR = f*NFREFMAX;

	FACEc->VL->neigh_f[FACEc->VfL]   = (FACEc->VfR)/NFREFMAX;
	FACEc->VR->neigh_f[FACEc->VfR] = (FACEc->VfL)/NFREFMAX;

	get_Indsf(FACEc,&sfIn,&sfOut);
	VOLUMEc->FACE[sfOut] = FACEc;
}

static void get_Indsf(struct S_FACE *FACE, unsigned int *sfIn, unsigned int *sfOut)
{
	unsigned int i, VType, Vf, f, Vfl, sf[2];

	for (i = 0; i < 2; i++) {
		if (i == 0) {
			Vf    = FACE->VfL;
			VType = FACE->VL->type;
		} else {
			Vf    = FACE->VfR;
			VType = FACE->VR->type;
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
				printf("Error: Unsupported Vfl in case %d of get_Indsf (%d).\n",VType,i), EXIT_MSG;
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
				printf("Error: Unsupported Vfl in case %d of get_Indsf (%d).\n",VType,i), EXIT_MSG;
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
				printf("Error: Unsupported Vfl in case %d of get_Indsf (%d).\n",VType,i), EXIT_MSG;
				break;
			}
			break;
		case WEDGE:
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
				printf("Error: Unsupported Vfl in case %d of get_Indsf (%d).\n",VType,i), EXIT_MSG;
				break;
			}
			break;
		case PYR:
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
				printf("Error: Unsupported Vfl in case %d of get_Indsf (%d).\n",VType,i), EXIT_MSG;
				break;
			}
			break;
		default:
			printf("Error: Unsupported VType in get_Indsf (%d).\n",i), EXIT_MSG;
			break;
		}
		sf[i] += f*NSUBFMAX;
	}
	*sfIn  = sf[0];
	*sfOut = sf[1];
}

static void update_memory_FACE(struct S_FACE *const FACE);

static unsigned int is_VOLUME_VL(const unsigned int indexgVOLUME, const unsigned int indexgVL)
{
	if (indexgVOLUME == indexgVL)
		return 1;
	return 0;
}

static void coarse_update(struct S_VOLUME *VOLUME)
{
	// Initialize DB Parameters
	unsigned int TETrefineType = DB.TETrefineType;

	// Standard datatypes
	unsigned int i, iMax, f, Nf, VType, vhMin, vhMax, sf, sfMax, sfMax_i, fIn, fOut,
	             IndVc[NVISUBFMAX],Indsf[NVISUBFMAX], recombine_FACEs,
	             dummy_ui;

	struct S_ELEMENT *ELEMENT;
	struct S_VOLUME  *VOLUMEc, **VOLUMEc_list, *VL, *VR;
	struct S_FACE   *FACE, *FACEp;

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
			// NOTE: The computation of FACE->VfR must be modified if anisotropic refinement is enabled.
			//       (ToBeDeleted)
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
			for (sf = 0; sf < sfMax; sf++)
				Indsf[sf] = f*NSUBFMAX;
			if (TETrefineType == TET8 || TETrefineType == TET12) {
				switch (f) {
				default: IndVc[0] = 2; IndVc[1] = 3; IndVc[2] = 1; IndVc[3] = 4; break;
				case 1:  IndVc[0] = 2; IndVc[1] = 3; IndVc[2] = 0; IndVc[3] = 5; break;
				case 2:  IndVc[0] = 0; IndVc[1] = 1; IndVc[2] = 3; IndVc[3] = 6; break;
				case 3:  IndVc[0] = 0; IndVc[1] = 1; IndVc[2] = 2; IndVc[3] = 7; break;
				}
			} else if (TETrefineType == TET6) {
				switch (f) {
				default: IndVc[0] = 2; IndVc[1] = 3; IndVc[2] = 1; IndVc[3] = 5; Indsf[sfMax-1] = 1*NSUBFMAX; break;
				case 1:  IndVc[0] = 2; IndVc[1] = 3; IndVc[2] = 0; IndVc[3] = 5; Indsf[sfMax-1] = 0*NSUBFMAX; break;
				case 2:  IndVc[0] = 0; IndVc[1] = 1; IndVc[2] = 3; IndVc[3] = 4; break;
				case 3:  IndVc[0] = 0; IndVc[1] = 1; IndVc[2] = 2; IndVc[3] = 4; break;
				}
			} else {
				printf("Error: Unsupported.\n"), EXIT_MSG;
			}
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
		case WEDGE:
			// Supported: Isotropic refinement
			sfMax = 4;
			switch (f) {
			default: IndVc[0] = 1; IndVc[1] = 2; IndVc[2] = 5; IndVc[3] = 6; break;
			case 1:  IndVc[0] = 0; IndVc[1] = 2; IndVc[2] = 4; IndVc[3] = 6; break;
			case 2:  IndVc[0] = 0; IndVc[1] = 1; IndVc[2] = 4; IndVc[3] = 5; break;
			case 3:  IndVc[0] = 0; IndVc[1] = 1; IndVc[2] = 2; IndVc[3] = 3; break;
			case 4:  IndVc[0] = 4; IndVc[1] = 5; IndVc[2] = 6; IndVc[3] = 7; break;
			}
			for (sf = 0; sf < sfMax; sf++)
				Indsf[sf] = f*NSUBFMAX;
			break;
		case PYR:
			// Supported: Isotropic refinement
			sfMax = 4;

			for (sf = 0; sf < sfMax; sf++)
				Indsf[sf] = f*NSUBFMAX;

			switch (f) {
			default: IndVc[0] = 0; IndVc[1] = 2; IndVc[2] = 9; IndVc[3] = 4; Indsf[sfMax-1] = 1*NSUBFMAX; break;
			case 1:  IndVc[0] = 1; IndVc[1] = 3; IndVc[2] = 9; IndVc[3] = 5; Indsf[sfMax-1] = 0*NSUBFMAX; break;
			case 2:  IndVc[0] = 0; IndVc[1] = 1; IndVc[2] = 9; IndVc[3] = 6; break;
			case 3:  IndVc[0] = 2; IndVc[1] = 3; IndVc[2] = 9; IndVc[3] = 7; break;
			case 4:  IndVc[0] = 0; IndVc[1] = 1; IndVc[2] = 2; IndVc[3] = 3; break;
			}
			break;
		default:
			printf("Error: Unsupported VType.\n"), EXIT_MSG;
			break;
		}

		recombine_FACEs = 0;
// Only need to check sf = 0 here as all sub FACEs should have the same level gap if coarsening was allowed.
// (ToBeDeleted).
		for (sf = 0; sf < sfMax; sf++) {
			FACE = VOLUMEc_list[IndVc[sf]]->FACE[Indsf[sf]];

			VL  = FACE->VL;
			VR = FACE->VR;

			if (VL == VR || VL->level != VR->level) {
				recombine_FACEs = 1;
				break;
			}
		}

		// Need to update FACE->VL and FACE->VR;
		if (recombine_FACEs) {
			VOLUME->NsubF[f] = 1;

			for (sf = 0; sf < sfMax; sf++) {
				VOLUMEc = VOLUMEc_list[IndVc[sf]];
				FACE = VOLUMEc->FACE[Indsf[sf]];

				FACE->update = 1;
				FACE->adapt_type = HCOARSE;
			}

			VOLUME->FACE[f*NSUBFMAX] = FACE->parent;
			FACEp = FACE->parent;
			if (is_VOLUME_VL(VOLUME->indexg,FACEp->VL->indexg)) {
				fOut = (FACEp->VfR)/NFREFMAX;
				FACEp->VR->NsubF[fOut] = 1;
				FACEp->VR->FACE[fOut*NSUBFMAX] = FACEp;
			} else {
				fIn = (FACEp->VfL)/NFREFMAX;
				FACEp->VL->NsubF[fIn] = 1;
				FACEp->VL->FACE[fIn*NSUBFMAX] = FACEp;
			}
//			update_memory_FACE(FACEp);
		} else {
			for (sf = 0; sf < sfMax; sf++) {
				VOLUMEc = VOLUMEc_list[IndVc[sf]];
				FACE = VOLUMEc->FACE[Indsf[sf]];

				if (is_VOLUME_VL(VOLUMEc->indexg,FACE->VL->indexg)) {
					FACE->VL  = FACE->VR;
					FACE->VfL = FACE->VfR;

					dummy_ui           = FACE->IndOrdLR;
					FACE->IndOrdLR = FACE->IndOrdRL;
					FACE->IndOrdRL = dummy_ui;

					FACE->VR  = VOLUME;
					FACE->VfR = f*NFREFMAX+sf+1; // CHANGE THIS TO BE GENERAL! (ToBeDeleted)

					VOLUME->FACE[f*NSUBFMAX+sf] = FACE;

					// Recompute XYZ_fS, n_fS, and detJF_fS
// Need not be recomputed, just reordered (and potentially negated) (ToBeDeleted)
					setup_FACE_XYZ(FACE);
					setup_normals(FACE);
				} else {
					FACE->VR  = VOLUME;
					FACE->VfR = f*NFREFMAX+sf+1;

					VOLUME->FACE[f*NSUBFMAX+sf] = FACE;
				}
			}
			VOLUME->NsubF[f] = sfMax;

		}
	}

	// Mark internal FACEs for updating
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
		if (TETrefineType == TET8) {
			sfMax_i = 8;

			IndVc[0] = 0; IndVc[1] = 1; IndVc[2] = 2; IndVc[3] = 3;
			IndVc[4] = 4; IndVc[5] = 4; IndVc[6] = 5; IndVc[7] = 5;

			Indsf[0] = 0; Indsf[1] = 1; Indsf[2] = 2; Indsf[3] = 3;
			Indsf[4] = 2; Indsf[5] = 3; Indsf[6] = 2; Indsf[7] = 3;
		} else if (TETrefineType == TET12) {
			sfMax_i = 16;

			IndVc[0]  = 0; IndVc[1]  = 1; IndVc[2]  = 2; IndVc[3]  = 3;
			IndVc[4]  = 4; IndVc[5]  = 4; IndVc[6]  = 4;
			IndVc[7]  = 5; IndVc[8]  = 5; IndVc[9]  = 5;
			IndVc[10] = 6; IndVc[11] = 6; IndVc[12] = 6;
			IndVc[13] = 7; IndVc[14] = 7; IndVc[15] = 7;

			Indsf[0]  = 0; Indsf[1]  = 1; Indsf[2]  = 2; Indsf[3]  = 3;
			Indsf[4]  = 1; Indsf[5]  = 2; Indsf[6]  = 3;
			Indsf[7]  = 0; Indsf[8]  = 2; Indsf[9]  = 3;
			Indsf[10] = 0; Indsf[11] = 1; Indsf[12] = 3;
			Indsf[13] = 0; Indsf[14] = 1; Indsf[15] = 2;
		} else if (TETrefineType == TET6) {
			sfMax_i = 5;

			IndVc[0] = 0; IndVc[1] = 1; IndVc[2] = 2; IndVc[3] = 3; IndVc[4] = 4;

			Indsf[0] = 0; Indsf[1] = 1; Indsf[2] = 2; Indsf[3] = 3; Indsf[4] = 4;
		} else {
			printf("Error: Unsupported.\n"), EXIT_MSG;
		}
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
	case WEDGE:
		sfMax_i = 10;

		IndVc[0] = 0; IndVc[1] = 0; IndVc[2] = 1; IndVc[3] = 1; IndVc[4] = 2;
		IndVc[5] = 2; IndVc[6] = 3; IndVc[7] = 4; IndVc[8] = 5; IndVc[9] = 6;

		Indsf[0] = 0; Indsf[1] = 4; Indsf[2] = 1; Indsf[3] = 4; Indsf[4] = 2;
		Indsf[5] = 4; Indsf[6] = 4; Indsf[7] = 0; Indsf[8] = 1; Indsf[9] = 2;
		for (i = 0; i < sfMax_i; i++)
			Indsf[i] *= NSUBFMAX;
		break;
	case PYR:
		sfMax_i = 13;

		IndVc[0]  = 0; IndVc[1]  = 0; IndVc[2]  = 1; IndVc[3] = 1; IndVc[4] = 2;
		IndVc[5]  = 2; IndVc[6]  = 3; IndVc[7]  = 3; IndVc[8] = 4; IndVc[9] = 5;
		IndVc[10] = 6; IndVc[11] = 7; IndVc[12] = 8;

		Indsf[0]  = 1; Indsf[1]  = 3; Indsf[2]  = 0; Indsf[3] = 3; Indsf[4] = 1;
		Indsf[5]  = 2; Indsf[6]  = 0; Indsf[7]  = 2; Indsf[8] = 0; Indsf[9] = 1;
		Indsf[10] = 3; Indsf[11] = 2; Indsf[12] = 4;
		for (i = 0; i < sfMax_i; i++)
			Indsf[i] *= NSUBFMAX;
		break;
	default:
		printf("Error: Unsupported VType.\n"), EXIT_MSG;
		break;
	}

	for (sf = 0; sf < sfMax_i; sf++) {
		VOLUMEc = VOLUMEc_list[IndVc[sf]];
		FACE = VOLUMEc->FACE[Indsf[sf]];
		FACE->update = 1;
		FACE->adapt_type = HDELETE;
	}

	free(VOLUMEc_list);
}

static void update_memory_FACE(struct S_FACE *const FACE)
{
	/*
	 *	Purpose:
	 *		Update amount of memory allocated to RHS/LHS arrays (used for non-vectorized functions).
	 */

	if (DB.Vectorized)
		EXIT_UNSUPPORTED;

	unsigned int const d    = DB.d,
	                   Nvar = DB.Nvar,
	                   Neq  = DB.Neq;

	// Left VOLUME
	unsigned int const NvnSL = FACE->VL->NvnS;
	if (NvnSL == 0)
		EXIT_UNSUPPORTED;

	// RHS/LHS
	if (FACE->RHSL != NULL)
		free(FACE->RHSL);
	FACE->RHSL = malloc(NvnSL*Nvar * sizeof *(FACE->RHSL)); // keep

	if (strstr(DB.SolverType,"Implicit")) {
		if (FACE->LHSLL != NULL)
			free(FACE->LHSLL);
		FACE->LHSLL = malloc(NvnSL*NvnSL*Nvar*Neq * sizeof *(FACE->LHSLL)); // keep
	}

	if (!(FACE->Boundary)) {
		unsigned int const NvnSR = FACE->VR->NvnS;
		if (FACE->RHSR != NULL)
			free(FACE->RHSR);
		FACE->RHSR = malloc(NvnSR*Nvar * sizeof *(FACE->RHSR)); // keep

		if (strstr(DB.SolverType,"Implicit")) {
			if (FACE->LHSLR != NULL)
				free(FACE->LHSLR);
			FACE->LHSLR = malloc(NvnSR*NvnSL*Nvar*Neq * sizeof *(FACE->LHSLR)); // keep
			if (FACE->LHSRL != NULL)
				free(FACE->LHSRL);
			FACE->LHSRL = malloc(NvnSL*NvnSR*Nvar*Neq * sizeof *(FACE->LHSRL)); // keep
			if (FACE->LHSRR != NULL)
				free(FACE->LHSRR);
			FACE->LHSRR = malloc(NvnSR*NvnSR*Nvar*Neq * sizeof *(FACE->LHSRR)); // keep
		}
	}

	// Other solver related arrays
	if (DB.Viscous) {
		for (size_t dim = 0; dim < d; dim++) {
			if (FACE->QhatL[dim] != NULL)
				free(FACE->QhatL[dim]);
			FACE->QhatL[dim] = malloc(NvnSL*Nvar * sizeof *(FACE->QhatL[dim])); // keep

			if (!FACE->Boundary) {
				unsigned int const NvnSR = FACE->VR->NvnS;
				if (FACE->QhatR[dim] != NULL)
					free(FACE->QhatR[dim]);
				FACE->QhatR[dim] = malloc(NvnSR*Nvar * sizeof *(FACE->QhatR[dim])); // keep
			}
		}

		if (strstr(DB.SolverType,"Implicit")) {
			for (size_t dim = 0; dim < d; dim++) {
				if (FACE->Boundary) {
					if (FACE->QhatL_WhatL[dim] != NULL)
						free(FACE->QhatL_WhatL[dim]);
					FACE->QhatL_WhatL[dim] = malloc(NvnSL*NvnSL*Nvar*Neq * sizeof *(FACE->QhatL_WhatL[dim])); // keep
				} else {
					unsigned int const NvnSR = FACE->VR->NvnS;
					if (FACE->QhatL_WhatL[dim] != NULL)
						free(FACE->QhatL_WhatL[dim]);
					FACE->QhatL_WhatL[dim] = malloc(NvnSL*NvnSL * sizeof *(FACE->QhatL_WhatL[dim])); // keep
					if (FACE->QhatL_WhatR[dim] != NULL)
						free(FACE->QhatL_WhatR[dim]);
					FACE->QhatL_WhatR[dim] = malloc(NvnSL*NvnSR * sizeof *(FACE->QhatL_WhatR[dim])); // keep
					if (FACE->QhatR_WhatL[dim] != NULL)
						free(FACE->QhatR_WhatL[dim]);
					FACE->QhatR_WhatL[dim] = malloc(NvnSR*NvnSL * sizeof *(FACE->QhatR_WhatL[dim])); // keep
					if (FACE->QhatR_WhatR[dim] != NULL)
						free(FACE->QhatR_WhatR[dim]);
					FACE->QhatR_WhatR[dim] = malloc(NvnSR*NvnSR * sizeof *(FACE->QhatR_WhatR[dim])); // keep
				}
			}
		}
	}
}

void update_memory_FACEs(void)
{
	for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next)
		update_memory_FACE(FACE);
}

static void free_memory_solver_FACE(struct S_FACE *const FACE)
{
	/*
	 *	Purpose:
	 *		Free memory associated with FACE solver arrays when a FACE is h-refined.
	 *
	 *	Comments:
	 *		This includes memory used to store the solution as well as memory used for RHS/LHS terms.
	 *		Memory addresses associated must be set to NULL after being freed such that it is not attempted to free them
	 *		again in update_memory_FACE in the case of the mesh being coarsened.
	 */

	if (!(FACE->adapt_type == HREFINE))
		EXIT_UNSUPPORTED;

	unsigned int const d = DB.d;

	// RHS/LHS
	free(FACE->RHSL);
	FACE->RHSL = NULL;

	if (strstr(DB.SolverType,"Implicit")) {
		free(FACE->LHSLL);
		FACE->LHSLL = NULL;
	}

	if (!FACE->Boundary) {
		free(FACE->RHSR);
		FACE->RHSR = NULL;

		if (strstr(DB.SolverType,"Implicit")) {
			free(FACE->LHSLR);
			FACE->LHSLR = NULL;
			free(FACE->LHSRL);
			FACE->LHSRL = NULL;
			free(FACE->LHSRR);
			FACE->LHSRR = NULL;
		}
	}

	if (DB.Viscous) {
		for (size_t dim = 0; dim < d; dim++) {
			free(FACE->QhatL[dim]);
			FACE->QhatL[dim] = NULL;

			if (!FACE->Boundary) {
				free(FACE->QhatR[dim]);
				FACE->QhatR[dim] = NULL;
			}
		}
		if (strstr(DB.SolverType,"Implicit")) {
			for (size_t dim = 0; dim < d; dim++) {
				free(FACE->QhatL_WhatL[dim]);
				FACE->QhatL_WhatL[dim] = NULL;

				if (!FACE->Boundary) {
					free(FACE->QhatL_WhatR[dim]);
					FACE->QhatL_WhatR[dim] = NULL;
					free(FACE->QhatR_WhatL[dim]);
					FACE->QhatR_WhatL[dim] = NULL;
					free(FACE->QhatR_WhatR[dim]);
					FACE->QhatR_WhatR[dim] = NULL;
				}
			}
		}
	}
}

void update_FACE_hp(void)
{
	/*
	 *	Comments:
	 *		Parameters to (potentially) modify:
	 *			P, typeInt, BC, indexg,
	 *			VL/VR, VfL/VfR, IndOrdLR/IndOrdRL,
	 *			n_fS, XYZ_fS, detJF_fS (ToBeModified (potentially to _fI))
	 */

	// Initialize DB Parameters
	unsigned int Adapt     = DB.Adapt,
	             NGF       = DB.NGF,
	             LevelsMax = DB.LevelsMax;

	// Standard datatypes
	unsigned int l, Nf, f, fh, Vfh, vh, sf, sfIn, sfOut, fhMax, sfMax, dummy_ui,
	             Indf, IndVLh, adapt_type, BC, internal_BC,
	             *NsubF;

	struct S_ELEMENT *ELEMENT;
	struct S_VOLUME  *VOLUME, *VOLUMEc, *VOLUMEp, *VL, *VR;
	struct S_FACE   *FACE, *FACEc, *FACEnext, *FACEtmp;

	// silence
	IndVLh = Vfh = 0;
	FACEc = NULL;

	switch (Adapt) {
	default: // ADAPT_HP
//      Pass adapt type to update_FACE_hp (ToBeModified)
		DB.Adapt = ADAPT_H; update_FACE_hp();
		DB.Adapt = ADAPT_P; update_FACE_hp();
		DB.Adapt = ADAPT_HP;
		break;
	case ADAPT_H:
		// HREFINE
		for (l = 0; l < LevelsMax; l++) {
			for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
				if (VOLUME->update && VOLUME->adapt_type == HREFINE && l == VOLUME->level) {
					NsubF = VOLUME->NsubF;

					// External FACEs
					ELEMENT = get_ELEMENT_type(VOLUME->type);
					Nf = ELEMENT->Nf;

					for (f = 0; f < Nf; f++) {
						Indf = f*NSUBFMAX;
						sfMax = NsubF[f];
						if (sfMax == 1) { // Create new FACEs between two VOLUMEs which were previously on the same level.
							FACE = VOLUME->FACE[Indf];
							FACE->update = 1;
							FACE->adapt_type = HREFINE;

							VL  = FACE->VL;
							VR = FACE->VR;
							BC   = FACE->BC;

							internal_BC = (BC == 0 || (BC % BC_STEP_SC > 50));

							fhMax = get_fhMax(VOLUME->type,VOLUME->hrefine_type);
							for (fh = 0; fh < fhMax; fh++) {
								if (fh) {
									FACEc->next = New_FACE();
									FACEc = FACEc->next;
								} else {
									FACEc = New_FACE();
									FACE->child0 = FACEc;
								}
								FACEc->update = 1;
								FACEc->parent = FACE;

								FACEc->indexg = NGF++;
								FACEc->level  = (VOLUME->level)+1;
								FACEc->BC       = BC;
								FACEc->Boundary = FACE->Boundary;

								// Ensuring that FACEc->Boundary is not being set to UINT_MAX
								if (FACE->Boundary != 0 && FACE->Boundary != 1)
									EXIT_UNSUPPORTED;

								// Find out if VOLUME == VL or VR
								// If condition can be outside of fh loop (ToBeDeleted)
								if (is_VOLUME_VL(VOLUME->indexg,VL->indexg)) {
									get_FACE_IndVL(FACE->VfL,fh,VL->type,&IndVLh,&Vfh);

									VOLUMEc = VOLUME->child0;
									for (vh = 0; vh < IndVLh; vh++)
										VOLUMEc = VOLUMEc->next;

									FACEc->VL  = VOLUMEc;
									FACEc->VfL = Vfh*NFREFMAX;
									FACEc->type = get_FACE_type(FACEc);

									if (internal_BC) {
										FACEc->VR  = VR;
										FACEc->VfR = get_FACE_VfR(fh,FACE->IndOrdRL,
										                                VL->neigh_f[f*NFREFMAX],FACEc->type);
									} else {
										FACEc->VR  = FACEc->VL;
										FACEc->VfR = FACEc->VfL;
									}

									FACEc->IndOrdLR = FACE->IndOrdLR;
									FACEc->IndOrdRL = FACE->IndOrdRL;
								} else {
									get_FACE_IndVL(FACE->VfR,fh,VR->type,&IndVLh,&Vfh);

									VOLUMEc = VOLUME->child0;
									for (vh = 0; vh < IndVLh; vh++)
										VOLUMEc = VOLUMEc->next;

									FACEc->VL  = VOLUMEc;
									FACEc->VfL = Vfh*NFREFMAX;
									FACEc->type = get_FACE_type(FACEc);

									FACEc->VR  = VL;
									FACEc->VfR = get_FACE_VfR(fh,FACE->IndOrdLR,VR->neigh_f[f*NFREFMAX],FACEc->type);

									FACEc->IndOrdLR = FACE->IndOrdRL;
									FACEc->IndOrdRL = FACE->IndOrdLR;

								}
								FACEc->VL->NsubF[(FACEc->VfL)/NFREFMAX] = 1;
								if (internal_BC)
									FACEc->VR->NsubF[(FACEc->VfR)/NFREFMAX] = fhMax;
								VOLUMEc->neigh_f[FACEc->VfL] = (FACEc->VfR)/NFREFMAX;

								get_Indsf(FACEc,&sfIn,&sfOut);

								FACEc->VL->FACE[sfIn] = FACEc;
								FACEc->VR->FACE[sfOut] = FACEc;
//								update_memory_FACE(FACEc);
							}
							free_memory_solver_FACE(FACE);
						} else { // Connect to existing FACEs
							// VOLUME = VR
							for (sf = 0; sf < sfMax; sf++) {
								FACEc = VOLUME->FACE[Indf+sf];

								// Only connectivity needs updating in FACEc->VR
								set_FACE_Out_External(FACEc,VOLUME);
							}
						}
					}

					// Internal FACEs (Always created)
					vh = 0;
					for (VOLUMEc = VOLUME->child0; VOLUMEc; VOLUMEc = VOLUMEc->next) {
						ELEMENT = get_ELEMENT_type(VOLUMEc->type);
						Nf = ELEMENT->Nf;

						for (f = 0; f < Nf; f++) {
							if (!VOLUMEc->FACE[f*NSUBFMAX]) {
								while (FACEc->next)
									FACEc = FACEc->next;
								FACEc->next = New_FACE();
								FACEc = FACEc->next;
								FACEc->update = 1;
								FACEc->indexg = NGF++;
								FACEc->level  = (VOLUME->level)+1;
								FACEc->BC     = 0; // Internal (Note: May have a curved edge in 3D)
								FACEc->Boundary = 0;

								VOLUMEc->NsubF[f] = 1;

								FACEc->VL = VOLUMEc;
								FACEc->VfL = f*NFREFMAX;
								set_FACE_Out(vh,f,FACEc,VOLUME);
//								update_memory_FACE(FACEc);
							}
						}
						vh++;
					}
				}
			}
		}
		DB.NGF = NGF;

		// Update FACE linked list (For HREFINE)
		// This is done inelegantly because it was required to keep the pointer to the parent FACE.

		// Fix list head if necessary
		FACE = DB.FACE;

		if (FACE->update) {
			adapt_type = FACE->adapt_type;
			if (adapt_type == HREFINE) {
				DB.FACE = FACE->child0;
				for (FACEc = DB.FACE; FACEc->next; FACEc = FACEc->next)
					;
				FACEc->next = FACE->next;
			}
		}

		// Fix remainder of list
		for (FACE = DB.FACE; FACE; FACE = FACE->next) {
			FACEnext = FACE->next;
			if (FACEnext && FACEnext->update) {
				adapt_type = FACEnext->adapt_type;
				if (adapt_type == HREFINE) {
					FACE->next = FACEnext->child0;
					for (FACEc = FACE->next; FACEc->next; FACEc = FACEc->next)
						;
					FACEc->next = FACEnext->next;
				}
			}
		}

		for (FACE = DB.FACE; FACE; FACE = FACE->next) {
			if (FACE->update) {
				FACE->update = 0;
				if (FACE->parent && FACE->parent->update)
					FACE->parent->update = 0;
				VL  = FACE->VL;
				VR = FACE->VR;

				FACE->P      = max(VL->P,VR->P);
				FACE->curved = max(VL->curved,VR->curved);
				if (!FACE->curved)
					FACE->typeInt = 's';
				else
					FACE->typeInt = 'c';

				// Compute XYZ_fS, n_fS, and detJF_fS (ToBeModified)
				setup_FACE_XYZ(FACE);
				setup_normals(FACE);
			}
		}

		// HCOARSE
		for (l = LevelsMax; l > 0; l--) {
		for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		if (VOLUME->update && VOLUME->adapt_type == HCOARSE && l == VOLUME->level) {
		if (VOLUME == VOLUME->parent->child0) {
			VOLUMEp = VOLUME->parent;

			if (VOLUMEp->update) {
				coarse_update(VOLUMEp);
			}
		}}}}

		// Fix list head if necessary
		FACE = DB.FACE;

		if (FACE->update) {
			adapt_type = FACE->adapt_type;
			if (adapt_type == HCOARSE) {
				// DB.FACE->parent can never be NULL.
				DB.FACE = FACE->parent;
				for (FACEc = FACE; FACEc->next->parent == DB.FACE; FACEc = FACEc->next)
					;
				while (FACEc->next && FACEc->next->adapt_type == HDELETE) {
					FACEtmp = FACEc->next->next;
					memory_destructor_F(FACEc->next);
					FACEc->next = FACEtmp;
				}
				DB.FACE->next = FACEc->next;
				FACEc->next = NULL;
			} else if (adapt_type == HDELETE) {
				printf("Error: Should not be entering HDELETE while updating FACE head.\n"), EXIT_MSG;
			}
		}

		// Fix remainder of list
		for (FACE = DB.FACE; FACE; FACE = FACE->next) {
			FACEnext = FACE->next;
			if (FACEnext && FACEnext->update) {
				adapt_type = FACEnext->adapt_type;

				if (adapt_type == HDELETE) {
					while (FACEnext) {
						if (FACEnext->adapt_type == HDELETE) {
							FACEtmp = FACEnext->next;
							memory_destructor_F(FACEnext);
							FACEnext = FACEtmp;
						} else {
							break;
						}
					}
					if (FACEnext)
						adapt_type = FACEnext->adapt_type;
				}

				if (adapt_type == HCOARSE) {
					FACE->next = FACEnext->parent;
					for (FACEc = FACEnext; FACEc->next && FACEc->next->parent == FACE->next; FACEc = FACEc->next)
						;
					while (FACEc->next && FACEc->next->adapt_type == HDELETE) {
						FACEtmp = FACEc->next->next;
						memory_destructor_F(FACEc->next);
						FACEc->next = FACEtmp;
					}

					FACE->next->next = FACEc->next;
					FACEc->next = NULL;
				} else {
					FACE->next = FACEnext;
				}
			}
		}

		// Slightly inefficient but did not put the effort to find why the individual calls commented above did not set
		// sufficient memory. ToBeModified.
		update_memory_FACEs();

		break;
	case ADAPT_P:
		/*	No modifications required for:
		 *		indexg, typeInt, BC
		 */

		for (FACE = DB.FACE; FACE; FACE = FACE->next) {
			VL  = FACE->VL;
			VR = FACE->VR;

			if (VL->update || VR->update) {
				FACE->P = max(VL->P,VR->P);

				// Ensure that VL is the highest order VOLUME
				if (VL->P < VR->P) {
					FACE->VL  = VR;
					FACE->VR = VL;

					dummy_ui     = FACE->VfL;
					FACE->VfL  = FACE->VfR;
					FACE->VfR = dummy_ui;

					dummy_ui           = FACE->IndOrdLR;
					FACE->IndOrdLR = FACE->IndOrdRL;
					FACE->IndOrdRL = dummy_ui;
				}

				// Recompute XYZ_fS, n_fS, and detJF_fS (ToBeModified)
				setup_FACE_XYZ(FACE);
				setup_normals(FACE);
			}
		}
		break;
	}
}
