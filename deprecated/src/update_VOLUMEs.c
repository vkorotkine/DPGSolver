// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

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
#include "S_FACE.h"

#include "adaptation.h"
#include "element_functions.h"
#include "matrix_functions.h"
#include "setup_ToBeCurved.h"
#include "setup_Curved.h"
#include "setup_geometry.h"
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
	unsigned int NvnGs, NvnGc, NvnS, *Nvve, NvnSP, NvnI, **VeMask;
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
	OPS->I_vGs_vGs   = ELEMENT->I_vGs_vGs[1][1];
	OPS->I_vGs_vGc   = ELEMENT->I_vGs_vGc[1][PNew][0];
	OPS->Ihat_vS_vS  = ELEMENT->Ihat_vS_vS[P][PNew];
	OPS->L2hat_vS_vS = ELEMENT->L2hat_vS_vS[P][PNew];
	OPS->VeMask      = ELEMENT->VeMask[1][2];
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

static void set_VOLUMEc_BC_Info(struct S_VOLUME *VOLUME, const unsigned int vh, unsigned int **BC)
{
	/*
	 *	Purpose:
	 *		Transfer BC information from parent to child VOLUME for FACEs (and EDGEs).
	 *
	 *	Comments:
	 *		For 3D elements need to update both FACE and EDGE BC information.
	 *
	 *	References:
	 *		Generate arrays here with python/documentation/h_refinement_info.py.
	 */

	unsigned int NF = 0, NE = 0, IndF[NFMAX], IndFP[NFMAX], IndE[NEMAX], IndEP[NEMAX], IndBC[NEMAX], F = 0, E = 1;

	switch (VOLUME->type) {
	case TRI:
		// FACE only
		if      (vh == 1) { NF = 2; IndFP[0] = 1; IndFP[1] = 2; }
		else if (vh == 2) { NF = 2; IndFP[0] = 0; IndFP[1] = 2; }
		else if (vh == 3) { NF = 2; IndFP[0] = 0; IndFP[1] = 1; }
		else if (vh == 4) { NF = 0; }
		else              { EXIT_UNSUPPORTED; }

		for (size_t f = 0; f < NF; f++)
			IndF[f] = IndFP[f];

		break;
	case QUAD:
		// FACE only
		if      (vh == 1) { NF = 2; IndFP[0] = 0; IndFP[1] = 2; }
		else if (vh == 2) { NF = 2; IndFP[0] = 1; IndFP[1] = 2; }
		else if (vh == 3) { NF = 2; IndFP[0] = 0; IndFP[1] = 3; }
		else if (vh == 4) { NF = 2; IndFP[0] = 1; IndFP[1] = 3; }
		else if (vh == 5) { NF = 3; IndFP[0] = 0; IndFP[1] = 2; IndFP[2] = 3; }
		else if (vh == 6) { NF = 3; IndFP[0] = 1; IndFP[1] = 2; IndFP[2] = 3; }
		else if (vh == 7) { NF = 3; IndFP[0] = 0; IndFP[1] = 1; IndFP[2] = 2; }
		else if (vh == 8) { NF = 3; IndFP[0] = 0; IndFP[1] = 1; IndFP[2] = 3; }
		else              { EXIT_UNSUPPORTED; }

		for (size_t f = 0; f < NF; f++)
			IndF[f] = IndFP[f];

		break;
	case TET:
		if (DB.TETrefineType == TET8) {
			// FACE
			if      (vh == 1) { NF = 3; IndFP[0] = 1; IndFP[1] = 2; IndFP[2] = 3; }
			else if (vh == 2) { NF = 3; IndFP[0] = 0; IndFP[1] = 2; IndFP[2] = 3; }
			else if (vh == 3) { NF = 3; IndFP[0] = 0; IndFP[1] = 1; IndFP[2] = 3; }
			else if (vh == 4) { NF = 3; IndFP[0] = 0; IndFP[1] = 1; IndFP[2] = 2; }
			else if (vh == 5) { NF = 1; IndFP[0] = 0; }
			else if (vh == 6) { NF = 1; IndFP[0] = 1; }
			else if (vh == 7) { NF = 1; IndFP[0] = 2; }
			else if (vh == 8) { NF = 1; IndFP[0] = 3; }
			else              { EXIT_UNSUPPORTED; }

			for (size_t f = 0; f < NF; f++)
				IndF[f] = IndFP[f];

			// EDGE
			if (vh == 1) {
				NE = 6;
				IndE[0]  = 0; IndE[1]  = 1; IndE[2]  = 2; IndE[3]  = 3; IndE[4]  = 4; IndE[5]  = 5;
				IndEP[0] = 3; IndEP[1] = 1; IndEP[2] = 2; IndEP[3] = 3; IndEP[4] = 2; IndEP[5] = 1;
				IndBC[0] = F; IndBC[1] = E; IndBC[2] = E; IndBC[3] = E; IndBC[4] = F; IndBC[5] = F;
			} else if (vh == 2) {
				NE = 6;
				IndE[0]  = 0; IndE[1]  = 1; IndE[2]  = 2; IndE[3]  = 3; IndE[4]  = 4; IndE[5]  = 5;
				IndEP[0] = 0; IndEP[1] = 3; IndEP[2] = 2; IndEP[3] = 2; IndEP[4] = 4; IndEP[5] = 0;
				IndBC[0] = E; IndBC[1] = F; IndBC[2] = E; IndBC[3] = F; IndBC[4] = E; IndBC[5] = F;
			} else if (vh == 3) {
				NE = 6;
				IndE[0]  = 0; IndE[1]  = 1; IndE[2]  = 2; IndE[3]  = 3; IndE[4]  = 4; IndE[5]  = 5;
				IndEP[0] = 0; IndEP[1] = 1; IndEP[2] = 3; IndEP[3] = 1; IndEP[4] = 0; IndEP[5] = 5;
				IndBC[0] = E; IndBC[1] = E; IndBC[2] = F; IndBC[3] = F; IndBC[4] = F; IndBC[5] = E;
			} else if (vh == 4) {
				NE = 6;
				IndE[0]  = 0; IndE[1]  = 1; IndE[2]  = 2; IndE[3]  = 3; IndE[4]  = 4; IndE[5]  = 5;
				IndEP[0] = 0; IndEP[1] = 1; IndEP[2] = 2; IndEP[3] = 3; IndEP[4] = 4; IndEP[5] = 5;
				IndBC[0] = F; IndBC[1] = F; IndBC[2] = F; IndBC[3] = E; IndBC[4] = E; IndBC[5] = E;
			} else if (vh == 5) {
				NE = 5;
				IndE[0]  = 0; IndE[1]  = 1; IndE[2]  = 3; IndE[3]  = 4; IndE[4]  = 5;
				IndEP[0] = 0; IndEP[1] = 2; IndEP[2] = 3; IndEP[3] = 0; IndEP[4] = 0;
				IndBC[0] = F; IndBC[1] = F; IndBC[2] = F; IndBC[3] = F; IndBC[4] = F;
			} else if (vh == 6) {
				NE = 5;
				IndE[0]  = 0; IndE[1]  = 1; IndE[2]  = 3; IndE[3]  = 4; IndE[4]  = 5;
				IndEP[0] = 2; IndEP[1] = 1; IndEP[2] = 1; IndEP[3] = 3; IndEP[4] = 1;
				IndBC[0] = F; IndBC[1] = F; IndBC[2] = F; IndBC[3] = F; IndBC[4] = F;
			} else if (vh == 7) {
				NE = 5;
				IndE[0]  = 0; IndE[1]  = 1; IndE[2]  = 2; IndE[3]  = 3; IndE[4]  = 4;
				IndEP[0] = 1; IndEP[1] = 0; IndEP[2] = 2; IndEP[3] = 2; IndEP[4] = 2;
				IndBC[0] = F; IndBC[1] = F; IndBC[2] = F; IndBC[3] = F; IndBC[4] = F;
			} else if (vh == 8) {
				NE = 5;
				IndE[0]  = 0; IndE[1]  = 1; IndE[2]  = 2; IndE[3]  = 3; IndE[4]  = 4;
				IndEP[0] = 3; IndEP[1] = 3; IndEP[2] = 3; IndEP[3] = 0; IndEP[4] = 1;
				IndBC[0] = F; IndBC[1] = F; IndBC[2] = F; IndBC[3] = F; IndBC[4] = F;
			} else {
				EXIT_UNSUPPORTED;
			}
		} else {
			EXIT_UNSUPPORTED;
		}
		break;
	case HEX:
		// FACE
		if      (vh == 1) { NF = 3; IndFP[0] = 0; IndFP[1] = 2; IndFP[2] = 4; }
		else if (vh == 2) { NF = 3; IndFP[0] = 1; IndFP[1] = 2; IndFP[2] = 4; }
		else if (vh == 3) { NF = 3; IndFP[0] = 0; IndFP[1] = 3; IndFP[2] = 4; }
		else if (vh == 4) { NF = 3; IndFP[0] = 1; IndFP[1] = 3; IndFP[2] = 4; }
		else if (vh == 5) { NF = 3; IndFP[0] = 0; IndFP[1] = 2; IndFP[2] = 5; }
		else if (vh == 6) { NF = 3; IndFP[0] = 1; IndFP[1] = 2; IndFP[2] = 5; }
		else if (vh == 7) { NF = 3; IndFP[0] = 0; IndFP[1] = 3; IndFP[2] = 5; }
		else if (vh == 8) { NF = 3; IndFP[0] = 1; IndFP[1] = 3; IndFP[2] = 5; }
		else              { EXIT_UNSUPPORTED; }

		for (size_t f = 0; f < NF; f++)
			IndF[f] = IndFP[f];

		// EDGE
		NE = 9;
		if (vh == 1) {
			IndE[0]  = 0; IndE[1]  = 1; IndE[2]  = 2; IndE[3]  = 4; IndE[4]  = 5; IndE[5]  = 6; IndE[6]  = 8; IndE[7]  = 9; IndE[8]  = 10;
			IndEP[0] = 0; IndEP[1] = 4; IndEP[2] = 2; IndEP[3] = 4; IndEP[4] = 4; IndEP[5] = 0; IndEP[6] = 8; IndEP[7] = 2; IndEP[8] = 0;
			IndBC[0] = E; IndBC[1] = F; IndBC[2] = F; IndBC[3] = E; IndBC[4] = F; IndBC[5] = F; IndBC[6] = E; IndBC[7] = F; IndBC[8] = F;
		} else if (vh == 2) {
			IndE[0]  = 0; IndE[1]  = 1; IndE[2]  = 2; IndE[3]  = 4; IndE[4]  = 5; IndE[5]  = 7; IndE[6]  = 8; IndE[7]  = 9; IndE[8]  = 11;
			IndEP[0] = 0; IndEP[1] = 4; IndEP[2] = 2; IndEP[3] = 4; IndEP[4] = 5; IndEP[5] = 1; IndEP[6] = 2; IndEP[7] = 9; IndEP[8] = 1;
			IndBC[0] = E; IndBC[1] = F; IndBC[2] = F; IndBC[3] = F; IndBC[4] = E; IndBC[5] = F; IndBC[6] = F; IndBC[7] = E; IndBC[8] = F;
		} else if (vh == 3) {
			IndE[0]  = 0; IndE[1]  = 1; IndE[2]  = 3; IndE[3]  = 4; IndE[4]  = 5; IndE[5]  = 6; IndE[6]  = 8; IndE[7]  = 10; IndE[8]  = 11;
			IndEP[0] = 4; IndEP[1] = 1; IndEP[2] = 3; IndEP[3] = 4; IndEP[4] = 4; IndEP[5] = 0; IndEP[6] = 0; IndEP[7] = 10; IndEP[8] = 3;
			IndBC[0] = F; IndBC[1] = E; IndBC[2] = F; IndBC[3] = E; IndBC[4] = F; IndBC[5] = F; IndBC[6] = F; IndBC[7] = E; IndBC[8] = F;
		} else if (vh == 4) {
			IndE[0]  = 0; IndE[1]  = 1; IndE[2]  = 3; IndE[3]  = 4; IndE[4]  = 5; IndE[5]  = 7; IndE[6]  = 9; IndE[7]  = 10; IndE[8]  = 11;
			IndEP[0] = 4; IndEP[1] = 1; IndEP[2] = 3; IndEP[3] = 4; IndEP[4] = 5; IndEP[5] = 1; IndEP[6] = 1; IndEP[7] = 3; IndEP[8] = 11;
			IndBC[0] = F; IndBC[1] = E; IndBC[2] = F; IndBC[3] = F; IndBC[4] = E; IndBC[5] = F; IndBC[6] = F; IndBC[7] = F; IndBC[8] = E;
		} else if (vh == 5) {
			IndE[0]  = 0; IndE[1]  = 2; IndE[2]  = 3; IndE[3]  = 4; IndE[4]  = 6; IndE[5]  = 7; IndE[6]  = 8; IndE[7]  = 9; IndE[8]  = 10;
			IndEP[0] = 2; IndEP[1] = 2; IndEP[2] = 5; IndEP[3] = 0; IndEP[4] = 6; IndEP[5] = 5; IndEP[6] = 8; IndEP[7] = 2; IndEP[8] = 0;
			IndBC[0] = F; IndBC[1] = E; IndBC[2] = F; IndBC[3] = F; IndBC[4] = E; IndBC[5] = F; IndBC[6] = E; IndBC[7] = F; IndBC[8] = F;
		} else if (vh == 6) {
			IndE[0]  = 0; IndE[1]  = 2; IndE[2]  = 3; IndE[3]  = 5; IndE[4]  = 6; IndE[5]  = 7; IndE[6]  = 8; IndE[7]  = 9; IndE[8]  = 11;
			IndEP[0] = 2; IndEP[1] = 2; IndEP[2] = 5; IndEP[3] = 1; IndEP[4] = 5; IndEP[5] = 7; IndEP[6] = 2; IndEP[7] = 9; IndEP[8] = 1;
			IndBC[0] = F; IndBC[1] = E; IndBC[2] = F; IndBC[3] = F; IndBC[4] = F; IndBC[5] = E; IndBC[6] = F; IndBC[7] = E; IndBC[8] = F;
		} else if (vh == 7) {
			IndE[0]  = 1; IndE[1]  = 2; IndE[2]  = 3; IndE[3]  = 4; IndE[4]  = 6; IndE[5]  = 7; IndE[6]  = 8; IndE[7]  = 10; IndE[8]  = 11;
			IndEP[0] = 3; IndEP[1] = 5; IndEP[2] = 3; IndEP[3] = 0; IndEP[4] = 6; IndEP[5] = 5; IndEP[6] = 0; IndEP[7] = 10; IndEP[8] = 3;
			IndBC[0] = F; IndBC[1] = F; IndBC[2] = E; IndBC[3] = F; IndBC[4] = E; IndBC[5] = F; IndBC[6] = F; IndBC[7] = E; IndBC[8] = F;
		} else if (vh == 8) {
			IndE[0]  = 1; IndE[1]  = 2; IndE[2]  = 3; IndE[3]  = 5; IndE[4]  = 6; IndE[5]  = 7; IndE[6]  = 9; IndE[7]  = 10; IndE[8]  = 11;
			IndEP[0] = 3; IndEP[1] = 5; IndEP[2] = 3; IndEP[3] = 1; IndEP[4] = 5; IndEP[5] = 7; IndEP[6] = 1; IndEP[7] = 3; IndEP[8] = 11;
			IndBC[0] = F; IndBC[1] = F; IndBC[2] = E; IndBC[3] = F; IndBC[4] = F; IndBC[5] = E; IndBC[6] = F; IndBC[7] = F; IndBC[8] = E;
		} else {
			EXIT_UNSUPPORTED;
		}
		break;
	case WEDGE:
		// FACE
		if      (vh == 1) { NF = 3; IndFP[0] = 1; IndFP[1] = 2; IndFP[2] = 3; }
		else if (vh == 2) { NF = 3; IndFP[0] = 0; IndFP[1] = 2; IndFP[2] = 3; }
		else if (vh == 3) { NF = 3; IndFP[0] = 0; IndFP[1] = 1; IndFP[2] = 3; }
		else if (vh == 4) { NF = 1; IndFP[0] = 3; }
		else if (vh == 5) { NF = 3; IndFP[0] = 1; IndFP[1] = 2; IndFP[2] = 4; }
		else if (vh == 6) { NF = 3; IndFP[0] = 0; IndFP[1] = 2; IndFP[2] = 4; }
		else if (vh == 7) { NF = 3; IndFP[0] = 0; IndFP[1] = 1; IndFP[2] = 4; }
		else if (vh == 8) { NF = 1; IndFP[0] = 4; }
		else              { EXIT_UNSUPPORTED; }

		for (size_t f = 0; f < NF; f++)
			IndF[f] = IndFP[f];

		// EDGE
		if (vh == 1) {
			NE = 8;
			IndE[0]  = 0; IndE[1]  = 1; IndE[2]  = 2; IndE[3]  = 4; IndE[4]  = 5; IndE[5]  = 6; IndE[6]  = 7; IndE[7]  = 8;
			IndEP[0] = 3; IndEP[1] = 1; IndEP[2] = 2; IndEP[3] = 1; IndEP[4] = 2; IndEP[5] = 6; IndEP[6] = 2; IndEP[7] = 1;
			IndBC[0] = F; IndBC[1] = E; IndBC[2] = E; IndBC[3] = F; IndBC[4] = F; IndBC[5] = E; IndBC[6] = F; IndBC[7] = F;
		} else if (vh == 2) {
			NE = 8;
			IndE[0]  = 0; IndE[1]  = 1; IndE[2]  = 2; IndE[3]  = 3; IndE[4]  = 5; IndE[5]  = 6; IndE[6]  = 7; IndE[7]  = 8;
			IndEP[0] = 0; IndEP[1] = 3; IndEP[2] = 2; IndEP[3] = 0; IndEP[4] = 2; IndEP[5] = 2; IndEP[6] = 7; IndEP[7] = 0;
			IndBC[0] = E; IndBC[1] = F; IndBC[2] = E; IndBC[3] = F; IndBC[4] = F; IndBC[5] = F; IndBC[6] = E; IndBC[7] = F;
		} else if (vh == 3) {
			NE = 8;
			IndE[0]  = 0; IndE[1]  = 1; IndE[2]  = 2; IndE[3]  = 3; IndE[4]  = 4; IndE[5]  = 6; IndE[6]  = 7; IndE[7]  = 8;
			IndEP[0] = 0; IndEP[1] = 1; IndEP[2] = 3; IndEP[3] = 0; IndEP[4] = 1; IndEP[5] = 1; IndEP[6] = 0; IndEP[7] = 8;
			IndBC[0] = E; IndBC[1] = E; IndBC[2] = F; IndBC[3] = F; IndBC[4] = F; IndBC[5] = F; IndBC[6] = F; IndBC[7] = E;
		} else if (vh == 4) {
			NE = 6;
			IndE[0]  = 0; IndE[1]  = 1; IndE[2]  = 2; IndE[3]  = 6; IndE[4]  = 7; IndE[5]  = 8;
			IndEP[0] = 3; IndEP[1] = 3; IndEP[2] = 3; IndEP[3] = 0; IndEP[4] = 1; IndEP[5] = 2;
			IndBC[0] = F; IndBC[1] = F; IndBC[2] = F; IndBC[3] = F; IndBC[4] = F; IndBC[5] = F;
		} else if (vh == 5) {
			NE = 8;
			IndE[0]  = 1; IndE[1]  = 2; IndE[2]  = 3; IndE[3]  = 4; IndE[4]  = 5; IndE[5]  = 6; IndE[6]  = 7; IndE[7]  = 8;
			IndEP[0] = 1; IndEP[1] = 2; IndEP[2] = 4; IndEP[3] = 4; IndEP[4] = 5; IndEP[5] = 6; IndEP[6] = 2; IndEP[7] = 1;
			IndBC[0] = F; IndBC[1] = F; IndBC[2] = F; IndBC[3] = E; IndBC[4] = E; IndBC[5] = E; IndBC[6] = F; IndBC[7] = F;
		} else if (vh == 6) {
			NE = 8;
			IndE[0]  = 0; IndE[1]  = 2; IndE[2]  = 3; IndE[3]  = 4; IndE[4]  = 5; IndE[5]  = 6; IndE[6]  = 7; IndE[7]  = 8;
			IndEP[0] = 0; IndEP[1] = 2; IndEP[2] = 3; IndEP[3] = 4; IndEP[4] = 5; IndEP[5] = 2; IndEP[6] = 7; IndEP[7] = 0;
			IndBC[0] = F; IndBC[1] = F; IndBC[2] = E; IndBC[3] = F; IndBC[4] = E; IndBC[5] = F; IndBC[6] = E; IndBC[7] = F;
		} else if (vh == 7) {
			NE = 8;
			IndE[0]  = 0; IndE[1]  = 1; IndE[2]  = 3; IndE[3]  = 4; IndE[4]  = 5; IndE[5]  = 6; IndE[6]  = 7; IndE[7]  = 8;
			IndEP[0] = 0; IndEP[1] = 1; IndEP[2] = 3; IndEP[3] = 4; IndEP[4] = 4; IndEP[5] = 1; IndEP[6] = 0; IndEP[7] = 8;
			IndBC[0] = F; IndBC[1] = F; IndBC[2] = E; IndBC[3] = E; IndBC[4] = F; IndBC[5] = F; IndBC[6] = F; IndBC[7] = E;
		} else if (vh == 8) {
			NE = 6;
			IndE[0]  = 3; IndE[1]  = 4; IndE[2]  = 5; IndE[3]  = 6; IndE[4]  = 7; IndE[5]  = 8;
			IndEP[0] = 4; IndEP[1] = 4; IndEP[2] = 4; IndEP[3] = 0; IndEP[4] = 1; IndEP[5] = 2;
			IndBC[0] = F; IndBC[1] = F; IndBC[2] = F; IndBC[3] = F; IndBC[4] = F; IndBC[5] = F;
		} else {
			EXIT_UNSUPPORTED;
		}
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}

	// This is assuming that indices of external FACEs of reference refined elements are the same as those of the parent
	for (size_t f = 0; f < NF; f++)
		VOLUME->BC[0][IndF[f]] = BC[0][IndFP[f]];

	// Note: NE = 0 for d = 2
	for (size_t e = 0; e < NE; e++)
		VOLUME->BC[1][IndE[e]] = BC[IndBC[e]][IndEP[e]];
}

static void update_memory_VOLUME(struct S_VOLUME *const VOLUME)
{
	/*
	 *	Purpose:
	 *		Update amount of memory allocated to RHS/LHS arrays (used for non-vectorized functions).
	 */

	if (DB.Vectorized)
		EXIT_UNSUPPORTED;

	unsigned int const d    = DB.d,
	                   Nvar = DB.Nvar,
	                   Neq  = DB.Neq,
	                   NvnS = VOLUME->NvnS;

	if (NvnS == 0)
		EXIT_UNSUPPORTED;

	// RHS/LHS
	if (VOLUME->RHS != NULL)
		free(VOLUME->RHS);
	VOLUME->RHS = malloc(NvnS*Nvar * sizeof *(VOLUME->RHS)); // keep

	if (strstr(DB.SolverType,"Implicit")) {
		if (VOLUME->LHS != NULL)
			free(VOLUME->LHS);
		VOLUME->LHS = malloc(NvnS*NvnS*Nvar*Neq * sizeof *(VOLUME->LHS)); // keep

		if (DB.Viscous) {
			for (size_t dim = 0; dim < d; dim++) {
				if (VOLUME->LHSQ[dim] != NULL)
					free(VOLUME->LHSQ[dim]);
				VOLUME->LHSQ[dim] = malloc(NvnS*NvnS*Nvar*Neq * sizeof *(VOLUME->LHSQ[dim])); // keep
			}
		}
	}

	// Other solver related arrays
	if (DB.Viscous) {
		for (size_t dim = 0; dim < d; dim++) {
			if (VOLUME->Qhat[dim] != NULL)
				free(VOLUME->Qhat[dim]);
			VOLUME->Qhat[dim]  = malloc(NvnS*Nvar * sizeof *(VOLUME->Qhat[dim])); // keep
		}
	}
}

void update_memory_VOLUMEs(void)
{
	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next)
		update_memory_VOLUME(VOLUME);
}

static void free_memory_solver_VOLUME(struct S_VOLUME *const VOLUME)
{
	/*
	 *	Purpose:
	 *		Free memory associated with VOLUME solver arrays when a VOLUME is h-refined.
	 *
	 *	Comments:
	 *		This includes memory used to store the solution as well as memory used for RHS/LHS terms.
	 *		Memory addresses associated must be set to NULL after being freed such that it is not attempted to free them
	 *		again in update_memory_VOLUME in the case of the mesh being coarsened.
	 */

	if (!(VOLUME->adapt_type == HREFINE))
		EXIT_UNSUPPORTED;

	unsigned int const d = DB.d;

	free(VOLUME->RHS);
	VOLUME->RHS = NULL;
	if (strstr(DB.SolverType,"Implicit")) {
		free(VOLUME->LHS);
		VOLUME->LHS = NULL;

		if (DB.Viscous) {
			for (size_t dim = 0; dim < d; dim++) {
				free(VOLUME->LHSQ[dim]);
				VOLUME->LHSQ[dim] = NULL;
			}
		}
	}

	if (DB.Viscous) {
		free(VOLUME->What);
		for (size_t dim = 0; dim < d; dim++) {
			free(VOLUME->Qhat[dim]);
			VOLUME->Qhat[dim] = NULL;
		}
	}
// Currently, explicit solver is called even for implicit runs to start off so this must be reset (ToBeModified)
	if (1||strstr(DB.SolverType,"Explicit"))
		free(VOLUME->RES);
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
	unsigned int i, j, ve, dim, iMax, P, PNew, f, level, adapt_type, vh, vhMin, vhMax, VType, Nf, Nve, NveP2,
	             IndEhref, NvnGs[2], NvnGc[2], NvnS[2], NvnSP, NCols, update, maxP, *VeInfo, cVeCount, **VeMask;
	double       *I_vGs_vGc[2], *XYZ_vV, *XYZ_vVP2, *XYZ_S,
	             **Ihat_vS_vS, **I_vGs_vGs, **L2hat_vS_vS, *What, *RES, *WhatP, *WhatH, *RESP, *RESH, *dummyPtr_d;

	struct S_OPERATORS *OPS, *OPSp;
	struct S_ELEMENT   *ELEMENT;
	struct S_VOLUME    *VOLUME, *VOLUMEc, *VOLUMEp;

	// silence
	I_vGs_vGs = NULL;
	VOLUMEc   = NULL;

	OPS  = malloc(sizeof *OPS);  // free
	OPSp = malloc(sizeof *OPSp); // free

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		if (!VOLUME->Vadapt)
			continue;

		P     = VOLUME->P;
		level = VOLUME->level;
		adapt_type = VOLUME->adapt_type;

		switch(adapt_type) {
		case PREFINE:
			if (P < PMax)
				PNew = P+1;
			else
				EXIT_UNSUPPORTED;
			VOLUME->PNew = PNew;
			break;
		case PCOARSE:
			if (P >= 1)
				PNew = P-1;
			else
				EXIT_UNSUPPORTED;
			VOLUME->PNew = PNew;
			break;
		case HREFINE:
			if (level == LevelsMax)
				EXIT_UNSUPPORTED;
			VOLUME->PNew = P;
			break;
		case HCOARSE:
			if (level == 0)
				EXIT_UNSUPPORTED;
			VOLUME->PNew = P;
			break;
		default:
			EXIT_UNSUPPORTED;
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

			NvnS[1] = OPS->NvnS;
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

				XYZ_vV = VOLUME->XYZ_vV;
				XYZ_S  = malloc(NvnGc[0]*NCols * sizeof *XYZ_S); // keep
				mm_CTN_d(NvnGc[0],NCols,NvnGs[0],I_vGs_vGc[0],XYZ_vV,XYZ_S);

				free(VOLUME->XYZ_S);
				VOLUME->XYZ_S = XYZ_S;
				VOLUME->NvnG  = NvnGc[0];

				free(VOLUME->XYZ);
				if (strstr(MeshType,"ToBeCurved")) {
					free(VOLUME->XYZ_vVc);
					setup_ToBeCurved(VOLUME);
				} else if (strstr(MeshType,"Curved")) {
					setup_Curved(VOLUME);
				}
			}

			free(VOLUME->detJV_vI);
			free(VOLUME->C_vI);
			setup_geom_factors(VOLUME);

			// Project solution coefficients (and RES if applicable)
			NvnS[0]    = OPS->NvnS;
			NvnSP      = OPS->NvnSP;

			VOLUME->NvnS = NvnSP;

			WhatP = malloc(NvnSP*Nvar * sizeof *WhatP); // keep
			if (adapt_type == PREFINE) {
				Ihat_vS_vS = OPS->Ihat_vS_vS;
				mm_CTN_d(NvnSP,Nvar,NvnS[0],Ihat_vS_vS[0],VOLUME->What,WhatP);
			} else {
				L2hat_vS_vS = OPS->L2hat_vS_vS;
				mm_CTN_d(NvnSP,Nvar,NvnS[0],L2hat_vS_vS[0],VOLUME->What,WhatP);
			}
			free(VOLUME->What);
			VOLUME->What = WhatP;

//			if (strstr(DB.SolverType,"Explicit")) {
// Currently, explicit solver is called even for implicit runs to start off so this must be reset (ToBeModified)
			if (1||strstr(DB.SolverType,"Explicit")) {
				switch (DB.ExplicitSolverType) {
				case EULER:
					free(VOLUME->RES);
					VOLUME->RES = NULL;
					break;
				case RK3_SSP:
					free(VOLUME->RES);
					VOLUME->RES = malloc(NvnSP*Nvar * sizeof *(VOLUME->RES)); // keep
					break;
				case RK4_LS:
					RESP = malloc(NvnSP*Nvar * sizeof *RESP); // keep
					if (adapt_type == PREFINE)
						mm_CTN_d(NvnSP,Nvar,NvnS[0],OPS->Ihat_vS_vS[0],VOLUME->RES,RESP);
					else
						mm_CTN_d(NvnSP,Nvar,NvnS[0],OPS->L2hat_vS_vS[0],VOLUME->RES,RESP);

					free(VOLUME->RES);
					VOLUME->RES = RESP;
					break;
				default:
					EXIT_UNSUPPORTED;
					break;
				}
			}

			// Update memory required by the non-vectorized solver
			update_memory_VOLUME(VOLUME);

			free(VOLUME->MInv);
			VOLUME->MInv = NULL;
			break;
		case HREFINE:
			VType = VOLUME->type;
			What  = VOLUME->What;
			RES   = VOLUME->RES;


			NvnGs[0]     = OPS->NvnGs;
			NvnGc[0]     = OPS->NvnGc;
			NvnS[0]      = OPS->NvnS;
			I_vGs_vGs    = OPS->I_vGs_vGs;
			I_vGs_vGc[0] = OPS->I_vGs_vGc;
			VeMask       = OPS->VeMask;

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

				// Update geometry
				IndEhref = get_IndEhref(VType,vh);

				VOLUMEc->XYZ_vV = malloc(NvnGs[IndEhref]*d * sizeof *XYZ_vV); // keep

				// May not need VeInfo for new VOLUMEs if avoiding usage for treating curved geometry (ToBeDeleted)
				Nve    = ELEMENT->Nve;
				VeInfo = VOLUMEc->VeInfo;

				if (AC) {
					mm_CTN_d(NvnGs[IndEhref],NCols,NvnGs[0],I_vGs_vGs[vh],VOLUME->XYZ_vV,VOLUMEc->XYZ_vV);
					VOLUMEc->curved = 1;

					// Set VOLUME BC Information
					set_VOLUMEc_BC_Info(VOLUMEc,vh,VOLUME->BC);
				} else if (VOLUME->curved) {
					// Determined VeInfo for VOLUMEc
					cVeCount = 0;
					for (ve = 0; ve < Nve; ve++) {
						VeInfo[ve+Nve*0] = 1;
						VeInfo[ve+Nve*1] = 1;
						for (j = 0; j < NvnGs[0]; j++) {
							if (fabs(I_vGs_vGs[vh][ve*NvnGs[0]+j]-1.0) < EPS) {
								// Already existing vertex
								for (i = 0; i < NVEINFO; i++)
									VeInfo[ve+Nve*i] = VOLUME->VeInfo[j+Nve*i];
								break;
							} else if (fabs(I_vGs_vGs[vh][ve*NvnGs[0]+j]) > EPS) {
								if (!VOLUME->VeInfo[j+Nve*0]) { // If not curved
									VeInfo[ve+Nve*0] = 0;
									VeInfo[ve+Nve*1] = 0;
									VeInfo[ve+Nve*2] = UINT_MAX;
									VeInfo[ve+Nve*3] = 0;
									break;
								} else {
									VeInfo[ve+Nve*2] = VOLUME->VeInfo[j+NvnGs[0]*2];
									// This is incorrect (currently not used) for vertices shared by two curved surfaces.
//									VeInfo[ve+Nve*3] = VOLUME->VeInfo[j+NvnGs[0]*3];
									VeInfo[ve+Nve*3] = UINT_MAX;
								}
							}
						}
						if (VeInfo[ve])
							cVeCount++;
					}

					if (d == DMAX && cVeCount == 2)
						VOLUMEc->curved = 2; // Curved EDGE
					else
						VOLUMEc->curved = 1; // Curved FACE

					// Ensure that vertices are place on the curved boundaries
					NveP2 = ELEMENT->NveP2;

					if (vh == vhMin)
						setup_Curved_vertices(VOLUME);

					if (DB.TETrefineType == TET12)
						printf("Error: VeMask not correct for this case.\n"), EXIT_MSG;

					XYZ_vV   = VOLUMEc->XYZ_vV;
					XYZ_vVP2 = VOLUME->XYZ_vVP2;

					for (ve = 0; ve < Nve; ve++) {
						for (dim = 0; dim < d; dim++)
							XYZ_vV[ve+Nve*dim] = XYZ_vVP2[VeMask[vh][ve]+NveP2*dim];
					}

					// Set VOLUME BC Information
					set_VOLUMEc_BC_Info(VOLUMEc,vh,VOLUME->BC);
				} else {
					mm_CTN_d(NvnGs[IndEhref],NCols,NvnGs[0],I_vGs_vGs[vh],VOLUME->XYZ_vV,VOLUMEc->XYZ_vV);
					for (ve = 0; ve < Nve; ve++) {
						VeInfo[ve+Nve*0] = 0;
						VeInfo[ve+Nve*1] = 0;
						VeInfo[ve+Nve*2] = UINT_MAX;
					}
					VOLUMEc->curved = 0;
				}

				XYZ_vV = VOLUMEc->XYZ_vV;
				if (!VOLUMEc->curved) {
					VOLUMEc->NvnG  = NvnGs[IndEhref];
					VOLUMEc->XYZ_S = malloc(NvnGs[IndEhref]*NCols * sizeof *XYZ_S); // keep
					XYZ_S = VOLUMEc->XYZ_S;
					for (unsigned int i = 0, iMax = NCols*NvnGs[IndEhref]; i < iMax; i++)
						XYZ_S[i] = XYZ_vV[i];
					setup_straight(VOLUMEc);
				} else {
					VOLUMEc->NvnG = NvnGc[IndEhref];

					VOLUMEc->XYZ_S = malloc(NvnGc[IndEhref]*NCols * sizeof *XYZ_S); // keep
					mm_CTN_d(NvnGc[IndEhref],NCols,NvnGs[IndEhref],I_vGs_vGc[IndEhref],XYZ_vV,VOLUMEc->XYZ_S);

					if (strstr(MeshType,"ToBeCurved"))
						setup_ToBeCurved(VOLUMEc);
					else if (strstr(MeshType,"Curved")) {
						setup_Curved(VOLUMEc);
					}
				}
				setup_geom_factors(VOLUMEc);

				// Project solution coefficients (and RES if applicable)
				Ihat_vS_vS = OPS->Ihat_vS_vS;
				VOLUMEc->NvnS = NvnS[IndEhref];

				WhatH = malloc(NvnS[IndEhref]*Nvar * sizeof *WhatH); // keep
				mm_CTN_d(NvnS[IndEhref],Nvar,NvnS[0],Ihat_vS_vS[vh],VOLUME->What,WhatH);
				VOLUMEc->What = WhatH;

//			if (strstr(DB.SolverType,"Explicit")) {
// Currently, explicit solver is called even for implicit runs to start off so this must be reset (ToBeModified)
				if (1||strstr(DB.SolverType,"Explicit")) {
					switch (DB.ExplicitSolverType) {
					case EULER:
						free(VOLUME->RES);
						VOLUME->RES = NULL;
						break;
					case RK3_SSP:
						free(VOLUMEc->RES);
						VOLUMEc->RES = malloc(NvnS[IndEhref]*Nvar * sizeof *(VOLUMEc->RES));  // keep
						break;
					case RK4_LS:
						RESH  = malloc(NvnS[IndEhref]*Nvar * sizeof *RESH);  // keep
						mm_CTN_d(NvnS[IndEhref],Nvar,NvnS[0],Ihat_vS_vS[vh],VOLUME->RES,RESH);
						free(VOLUMEc->RES);
						VOLUMEc->RES = RESH;
						break;
					default:
						EXIT_UNSUPPORTED;
						break;
					}
				}

				update_memory_VOLUME(VOLUMEc);
			}

			free_memory_solver_VOLUME(VOLUME);
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

					// If the parent VOLUME is curved and its order will be changed, its geometry must be updated.
					if (VOLUMEp->P != maxP && VOLUMEp->curved) {
						VOLUMEp->P    = maxP;
						VOLUMEp->PNew = maxP;

						init_ops(OPSp,VOLUMEp,0);

						NCols = d;

						XYZ_vV = VOLUMEp->XYZ_vV;
						XYZ_S  = malloc(OPSp->NvnGc*NCols * sizeof *XYZ_S); // keep
						mm_CTN_d(OPSp->NvnGc,NCols,OPSp->NvnGs,OPSp->I_vGs_vGc,XYZ_vV,XYZ_S);

						free(VOLUMEp->XYZ_S);
						VOLUMEp->XYZ_S = XYZ_S;
						VOLUMEp->NvnG  = OPSp->NvnGc;

						free(VOLUMEp->XYZ);
						if (strstr(MeshType,"ToBeCurved")) {
							free(VOLUMEp->XYZ_vVc);
							setup_ToBeCurved(VOLUMEp);
						} else if (strstr(MeshType,"Curved")) {
							setup_Curved(VOLUMEp);
						} else {
							EXIT_UNSUPPORTED;
						}

						free(VOLUMEp->detJV_vI);
						free(VOLUMEp->C_vI);
						setup_geom_factors(VOLUMEp);
					}

					VOLUMEp->P    = maxP;
					VOLUMEp->PNew = maxP;
					VOLUMEp->adapt_type = HCOARSE;

					// Project solution coefficients (and RES if applicable)
					NvnS[0] = OPS->NvnS;
					L2hat_vS_vS = OPS->L2hat_vS_vS;

					NCols = d;

					dummyPtr_d = malloc(NvnS[0]*Nvar * sizeof *dummyPtr_d); // free

					What = calloc(NvnS[0]*Nvar , sizeof *What); // keep

// Currently, explicit solver is called even for implicit runs to start off so this must be reset (ToBeModified)
					RES = NULL;
					if (1||strstr(DB.SolverType,"Explicit"))
						RES = calloc(NvnS[0]*Nvar , sizeof *RES); // keep

					VOLUMEc = VOLUME;
					for (vh = vhMin; vh <= vhMax; vh++) {
						IndEhref = get_IndEhref(VOLUMEp->type,vh);
						if (vh > vhMin)
							VOLUMEc = VOLUMEc->next;

						WhatH = VOLUMEc->What;
						mm_CTN_d(NvnS[0],Nvar,NvnS[IndEhref],L2hat_vS_vS[vh],WhatH,dummyPtr_d);
						for (i = 0, iMax = NvnS[0]*Nvar; i < iMax; i++)
							What[i] += dummyPtr_d[i];

						if (1||strstr(DB.SolverType,"Explicit")) {
							switch (DB.ExplicitSolverType) {
							case EULER:
								free(RES); RES = NULL;
								break;
							case RK3_SSP:
								; // Do nothing
								break;
							case RK4_LS:
								mm_CTN_d(NvnS[0],Nvar,NvnS[IndEhref],L2hat_vS_vS[vh],VOLUMEc->RES,dummyPtr_d);
								for (i = 0, iMax = NvnS[0]*Nvar; i < iMax; i++)
									RES[i] += dummyPtr_d[i];
								break;
							}
						}
					}
					free(dummyPtr_d);

					VOLUMEp->What = What;
					if (1||strstr(DB.SolverType,"Explicit"))
						VOLUMEp->RES = RES;
					update_memory_VOLUME(VOLUMEp);
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
	free(OPS);
	free(OPSp);
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
	/*
	 *	Purpose:
	 *		Compute the inverse of the mass matrix.
	 *
	 *	Comments:
	 *		The invserse is computed as it does not change unless the VOLUME undergoes hp adaptation and is thus
	 *		generally reused multiple times before being recomputed (if it is recomputed at all).
	 */
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

	if (DB.Collocated) {
		double *const MInv_diag = malloc(NvnS * sizeof *MInv_diag); // keep
		for (size_t i = 0; i < NvnS; i++)
			MInv_diag[i] = 1.0/wdetJV_vI[i];

		VOLUME->MInv_diag = MInv_diag;
		MInv = diag_d(MInv_diag,NvnS);
	} else {
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
	}

	free(VOLUME->MInv);
	VOLUME->MInv = MInv;

	free(OPS);
}

void update_VOLUME_Ops(void)
{
	// Initialize DB Parameters
	unsigned int Collocated  = DB.Collocated;

	// Standard datatypes
	struct S_VOLUME    *VOLUME;

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		if (VOLUME->update || VOLUME->MInv == NULL) {
			VOLUME->update = 0;
			if (strstr(DB.SolverType,"Explicit")) {
				if (!Collocated)
					compute_inverse_mass(VOLUME);
			} else {
				// Inverse mass matrix needed for viscous numerical flux.
				if (DB.Viscous)
					compute_inverse_mass(VOLUME);
			}
		}
	}
}

void update_VOLUME_finalize(void)
{
	unsigned int NV = 0;
	unsigned int VfL, VfR, fL, fR;

	struct S_VOLUME *VOLUME, *VL, *VR;
	struct S_FACE  *FACE;

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		VOLUME->indexg = NV++;
		VOLUME->Vadapt = 0;
		VOLUME->adapt_type = UINT_MAX; // Done for debugging only
		VOLUME->update = 0;
	}

	DB.NV = NV;
	DB.NVglobal = NV;

	for (FACE = DB.FACE; FACE; FACE = FACE->next) {
		VL  = FACE->VL;
		VfL = FACE->VfL;
		fL  = VfL/NFREFMAX;

		VR  = FACE->VR;
		VfR = FACE->VfR;
		fR  = VfR/NFREFMAX;

		FACE->Boundary = !((VL->indexg != VR->indexg) || (VL->indexg == VR->indexg && fL != fR));

		VL->neigh[VfL] = VR->indexg;
		VR->neigh[VfR] = VL->indexg;

		if (abs((int) VL->level - (int) VR->level) > 1.0) {
			printf("%d %d %d\n",VL->indexg,VR->indexg,VL->parent->indexg);
			printf("Error: Adjacent VOLUMEs are more than 1-irregular.\n"), EXIT_MSG;
		}
	}
}
