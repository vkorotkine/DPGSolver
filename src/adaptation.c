// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "adaptation.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"

#include "element_functions.h"
#include "array_sort.h"
#include "update_VOLUMEs.h"
#include "update_FACEs.h"
#include "setup_geometry.h"
#include "memory_free.h"

#include "array_print.h"

/*
 *	Purpose:
 *		Determine which VOLUMEs should be adapted based on specified error indicator and fixed fraction.
 *
 *	Comments:
 *		For the moment, adaptation is based only on the residual error. (ToBeDeleted)
 *		Don't forget to include MPI here to make sure that the fixed fraction applies to the global adaptation.
 *		(ToBeDeleted)
 *
 *	Notation:
 *
 *	References:
 */

struct S_VInfo {
	unsigned int *neigh, *type, *vh_range, *fh_range, *href_type, *p_levels, *h_levels, *h_siblings, *h_forbid_c,
	             *hp_refine, *hp_coarse, *hp_coarse_l, *hp_update;
} *VInfo;

void get_PS_range(size_t *const PSMin, size_t *const PSMax)
{
	/*
	 *	Purpose:
	 *		Return the range of orders used for the solution.
	 */

	// Initialize DB Parameters
	unsigned int PMax    = DB.PMax,
	             PGlobal = DB.PGlobal,
	             Adapt   = DB.Adapt;

	switch (Adapt) {
		default: // ADAPT_P or ADAPT_HP
			*PSMin = 0;
			*PSMax = PMax;
			break;
		case ADAPT_0:
		case ADAPT_H:
			*PSMin = PGlobal;
			*PSMax = PGlobal;
			break;
	}
}

void get_Pb_range(unsigned int const P, size_t *const PbMin, size_t *const PbMax)
{
	/*
	 *	Purpose:
	 *		Return the range of orders used for operators which interpolate between different orders.
	 *
	 *	Comments:
	 *		For Adapt == ADAPT_HP, the full range must be available such that FACE orders are acceptable if
	 *		h-coarsening is applied to a single neighbouring VOLUME having a range of orders.
	 */

	// Initialize DB Parameters
	unsigned int PMax  = DB.PMax,
	             Adapt = DB.Adapt;

	switch (Adapt) {
		default: // ADAPT_HP
			*PbMin = 0;
			*PbMax = PMax;
			break;
		case ADAPT_P:
			if      (P == 0)    *PbMin = P,   *PbMax = P+1;
			else if (P == PMax) *PbMin = P-1, *PbMax = PMax;
			else                *PbMin = P-1, *PbMax = P+1;
			break;
		case ADAPT_0:
		case ADAPT_H:
			*PbMin = P;
			*PbMax = P;
			break;
	}
}

void get_vh_range(const struct S_VOLUME *VOLUME, unsigned int *vhMin, unsigned int *vhMax)
{
	// Initialize DB Parameters
	unsigned int TETrefineType = DB.TETrefineType;

	// Standard datatypes
	unsigned int VType;

	VType = VOLUME->type;
//	href_type = VOLUME->hrefine_type;

	switch (VType) {
	case TRI:
	case QUAD:
		// Supported href_type: 0 (Isotropic)
		*vhMin = 1;
		*vhMax = 4;
		break;
	case TET:
		// Supported href_type: 0 (Isotropic - 1st internal TET orientation)
		*vhMin = 1;
		if (TETrefineType == TET8)
			*vhMax = 8;
		else if (TETrefineType == TET12)
			*vhMax = 12;
		else if (TETrefineType == TET6)
			*vhMax = 6;
		else
			printf("Error: Unsupported.\n"), EXIT_MSG;
		break;
	case HEX:
		// Supported href_type: 0 (Isotropic)
		*vhMin = 1;
		*vhMax = 8;
		break;
	case WEDGE:
		// Supported href_type: 0 (Isotropic)
		*vhMin = 1;
		*vhMax = 8;
		break;
	case PYR:
		// Supported href_type: 0 (Isotropic)
		*vhMin = 1;
		*vhMax = 10;
		break;
	default:
		printf("Error: Unsupported VType.\n"), EXIT_MSG;
		break;
	}
}

void get_fh_range(const struct S_VOLUME *VOLUME, const unsigned int f, unsigned int *fhMin, unsigned int *fhMax)
{
	unsigned int VType, NsubF;

	NsubF = VOLUME->NsubF[f];

	if (NsubF == 1) {
		*fhMin = 0;
		*fhMax = 0;
	} else {
		VType = VOLUME->type;
//		href_type = VOLUME->hrefine_type;

		switch (VType) {
		case TRI:
		case QUAD:
			// Supported href_type: 0 (Isotropic)
			*fhMin = 1;
			*fhMax = 2;
			break;
		case TET:
		case HEX:
		case WEDGE:
		case PYR:
			// Supported href_type: 0 (Isotropic)
			*fhMin = 1;
			*fhMax = 4;
			break;
		default:
			printf("Error: Unsupported VType.\n"), EXIT_MSG;
			break;
		}
	}
}

unsigned int get_VOLUMEc_type(const unsigned int VType, const unsigned int vh)
{
	// Initialize DB Parameters
	unsigned int TETrefineType = DB.TETrefineType;

	switch (VType) {
	case TET:
		if (TETrefineType == TET8 || TETrefineType == TET12) {
			return TET;
		} else if (TETrefineType == TET6) {
			if (vh < 5)
				return TET;
			else
				return PYR;
		} else {
			printf("Error: Unsupported.\n"), EXIT_MSG;
		}
		break;
	case PYR:
		if (vh < 5 || vh > 8)
			return PYR;
		else
			return TET;
		break;
	default:
		printf("Error: Unsupported VType.\n"), EXIT_MSG;
		break;
	}
}

unsigned int get_IndEhref(const unsigned int VType, const unsigned int vh)
{
	// Initialize DB Parameters
	unsigned int TETrefineType = DB.TETrefineType;

	switch (VType) {
	case POINT:
	case LINE:
	case TRI:
	case QUAD:
	case HEX:
	case WEDGE:
		return 0;
		break;
	case TET:
		if (TETrefineType == TET8 || TETrefineType == TET12) {
			return 0;
		} else if (TETrefineType == TET6) {
			if (vh < 5)
				return 0;
			else
				return 1;
		} else {
			printf("Error: Unsupported.\n"), EXIT_MSG;
		}
		break;
	case PYR:
		if (vh < 5 || vh > 8)
			return 0;
		else
			return 1;
		break;
	default:
		printf("Error: Unsupported VType.\n"), EXIT_MSG;
		break;
	}
}

static void check_levels_refine(const unsigned int indexg, const struct S_VInfo *VInfo, const char hp_type)
{
	// Standard datatypes
	unsigned int Nf, indexg_neigh, f, fh, fhMin, fhMax, Indf, Indfh,
	             *VNeigh, *VType, *hp_levels, *hp_refine_current, *fh_range;

	struct S_ELEMENT *ELEMENT;

	VNeigh = VInfo->neigh;
	VType  = VInfo->type;
	hp_refine_current = VInfo->hp_refine;

	switch (hp_type) {
	default: // p
		hp_levels = VInfo->p_levels;
		break;
	case 'h':
		hp_levels = VInfo->h_levels;
		break;
	}

	ELEMENT = get_ELEMENT_type(VType[indexg]);
	Nf = ELEMENT->Nf;

	hp_refine_current[indexg] = 1;
	hp_levels[indexg]++;

	Indf = indexg*NFREFMAX*NFMAX;
	switch (hp_type) {
	default: // p
		for (f = 0; f < Nf; f++) {
			indexg_neigh = VNeigh[Indf+f*NFREFMAX];
			if (!hp_refine_current[indexg_neigh] && ((int) hp_levels[indexg] - (int) hp_levels[indexg_neigh]) > 1)
				check_levels_refine(indexg_neigh,VInfo,'p');
		}
		break;
	case 'h':
		fh_range = VInfo->fh_range;
		for (f = 0; f < Nf; f++) {
			Indfh = Indf + f*NFREFMAX;
			fhMin = fh_range[(indexg*NFMAX+f)*2  ];
			fhMax = fh_range[(indexg*NFMAX+f)*2+1];
			for (fh = fhMin; fh <= fhMax; fh++) {
				indexg_neigh = VNeigh[Indfh+fh];
				if (!hp_refine_current[indexg_neigh] && ((int) hp_levels[indexg] - (int) hp_levels[indexg_neigh]) > 1)
					check_levels_refine(indexg_neigh,VInfo,'h');
			}
		}
		break;
	}
}

static void check_levels_coarse(const unsigned int indexg, const struct S_VInfo *VInfo, const char hp_type)
{
	// Standard datatypes
	unsigned int i, iMax, Nf, f, fh, vhMin, vhMax, fhMin, fhMax, Indf, Indfh, Indsib,
	             indexg_sib, indexg_neigh,
	             *VNeigh, *VType, *hp_levels, *hp_coarse_current, *fh_range, *vh_range, *h_siblings, *h_forbid_c;

	struct S_ELEMENT *ELEMENT;

	VNeigh = VInfo->neigh;
	VType  = VInfo->type;
	hp_coarse_current = VInfo->hp_coarse_l;

	switch (hp_type) {
	default: // p
		hp_levels = VInfo->p_levels;
		break;
	case 'h':
		hp_levels = VInfo->h_levels;
		break;
	}

	switch (hp_type) {
	default: // p
		ELEMENT = get_ELEMENT_type(VType[indexg]);
		Nf    = ELEMENT->Nf;

		hp_coarse_current[indexg] = 1;

		Indf = indexg*NFREFMAX*NFMAX;
		for (f = 0; f < Nf; f++) {
			indexg_neigh = VNeigh[Indf+f*NFREFMAX];
			if (!hp_coarse_current[indexg_neigh] && ((int) hp_levels[indexg_neigh] - (int) hp_levels[indexg]) > 0)
				check_levels_coarse(indexg_neigh,VInfo,'p');
		}
		break;
	case 'h':
		vh_range   = VInfo->vh_range;
		fh_range   = VInfo->fh_range;
		h_siblings = VInfo->h_siblings;
		h_forbid_c = VInfo->h_forbid_c;

		if (!h_forbid_c[indexg]) { // Check siblings and neighbours
			Indsib = indexg*NSIBMAX;
			vhMin = vh_range[indexg*2  ];
			vhMax = vh_range[indexg*2+1];
			for (i = 0, iMax = vhMax-vhMin; i <= iMax; i++) {
				indexg_sib = h_siblings[Indsib+i];

				ELEMENT = get_ELEMENT_type(VType[indexg_sib]);
				Nf    = ELEMENT->Nf;

				if (!hp_coarse_current[indexg_sib]) {
					hp_coarse_current[indexg_sib] = 1;

					Indf = indexg_sib*NFREFMAX*NFMAX;
					for (f = 0; f < Nf; f++) {
						Indfh = Indf + f*NFREFMAX;
						fhMin = fh_range[(indexg_sib*NFMAX+f)*2  ];
						fhMax = fh_range[(indexg_sib*NFMAX+f)*2+1];
						for (fh = fhMin; fh <= fhMax; fh++) {
							indexg_neigh = VNeigh[Indfh+fh];
							if (!hp_coarse_current[indexg_neigh] &&
								((int) hp_levels[indexg_neigh] - (int) hp_levels[indexg_sib]) > 0)
									check_levels_coarse(indexg_neigh,VInfo,'h');
						}
					}
				}
			}
		} else { // Check only neighbours
			hp_coarse_current[indexg] = 1;
		}
		break;
	}
}

void adapt_hp(void)
{
	// Initialize DB Parameters
	unsigned int DOF0       = DB.DOF0,
	             NV         = DB.NV,
	             NVglobal   = DB.NVglobal,
	             PMax       = DB.PMax,
	             LevelsMax  = DB.LevelsMax,
	             Adapt      = DB.Adapt;
	double       refine_frac = DB.refine_frac,
	             coarse_frac = DB.coarse_frac,
	             DOFcap_frac = DB.DOFcap_frac;

	// Standard datatypes
	unsigned int i, iMax, j, jMax, Nf, f, fh, Vf, fhMin, fhMax, vh, vhMin, vhMax,
	             iInd, DOF, NvnS, indexg, PMaxP1, LevelsMaxP1,
	             NFREFMAX_Total, refine_conflict, coarse_conflict, Indf,
	             *IndminRHS, *IndmaxRHS, *VNeigh, *VType_global, *vh_range, *fh_range,
	             *p_levels, *h_levels, *h_siblings, *h_forbid_coarse,
	             *hp_refine_current, *hp_coarse_current, *hp_coarse_current_local, *hp_coarse_current_err;
	double       minRHS, maxRHS, tmp_d, refine_frac_lim,
	             *RHS, *minRHS_Vec, *maxRHS_Vec, *minRHS_Vec_unsorted;

	struct S_ELEMENT *ELEMENT;
	struct S_VOLUME  *VOLUME, **VOLUME_Vec, *VOLUMEp, *VOLUMEc;

	// silence
	indexg = NVglobal;
	refine_frac_lim = 0.0;

	minRHS_Vec = malloc(NV * sizeof *minRHS_Vec); // free
	maxRHS_Vec = malloc(NV * sizeof *maxRHS_Vec); // free
	VOLUME_Vec = malloc(NV * sizeof *VOLUME_Vec); // free

// Change this from i to indexg (for MPI). Also then initialize min/maxRHS_Vec so that they can be sorted correctly.
	i = 0; DOF = 0;
	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		// Compute maxRHS in each VOLUME
		NvnS = VOLUME->NvnS;
		DOF += NvnS;

		RHS = VOLUME->RHS;

		RHS += DB.d*NvnS;

		minRHS = 1e10;
		maxRHS = 0.0;
		// Energy RHS
		for (iMax = NvnS; iMax--; ) {
			tmp_d = fabs(*RHS);
			if (tmp_d < minRHS)
				minRHS = tmp_d;
			if (tmp_d > maxRHS)
				maxRHS = tmp_d;
			RHS++;
		}
		minRHS_Vec[i] = minRHS;
		maxRHS_Vec[i] = maxRHS;
		VOLUME_Vec[i] = VOLUME;
		i++;
	}

	minRHS_Vec_unsorted = malloc(NV * sizeof *minRHS_Vec_unsorted); // free

	IndminRHS = malloc(NV * sizeof *IndminRHS); // free
	IndmaxRHS = malloc(NV * sizeof *IndmaxRHS); // free
	for (i = 0; i < NV; i++) {
		IndminRHS[i] = i;
		IndmaxRHS[i] = i;
		minRHS_Vec_unsorted[i] = minRHS_Vec[i];
	}

	array_sort_d(1,NV,minRHS_Vec,IndminRHS,'R','N');
	array_sort_d(1,NV,maxRHS_Vec,IndmaxRHS,'R','N');

// Add in MPI communication here for globally sorted array. (ToBeDeleted)

	// Make sure that the DOF cap is not exceeded
	refine_frac_lim = max(min((DOFcap_frac*DOF0)/DOF-1.0,refine_frac),0.0);
	printf("DOF, DOFcap: %d %d\n",DOF,(unsigned int) (DOFcap_frac*DOF0));
	if (refine_frac_lim < EPS) {
		printf("*** Warning: Consider raising DOFcap_frac. *** \n");
	}

	VInfo = malloc(sizeof *VInfo); // free
	VInfo->neigh       = malloc(NVglobal*NFMAX*NFREFMAX * sizeof *(VInfo->neigh));       // free
	VInfo->type        = calloc(NVglobal                , sizeof *(VInfo->type));        // free
	VInfo->vh_range    = malloc(NVglobal*2              * sizeof *(VInfo->vh_range));    // free
	VInfo->fh_range    = malloc(NVglobal*NFMAX*2        * sizeof *(VInfo->fh_range));    // free
	VInfo->href_type   = calloc(NVglobal                , sizeof *(VInfo->href_type));   // free
	VInfo->p_levels    = malloc(NVglobal                * sizeof *(VInfo->p_levels));    // free
	VInfo->h_levels    = malloc(NVglobal                * sizeof *(VInfo->h_levels));    // free
	VInfo->h_siblings  = malloc(NVglobal*NSIBMAX        * sizeof *(VInfo->h_siblings));  // free
	VInfo->h_forbid_c  = malloc(NVglobal                * sizeof *(VInfo->h_forbid_c));  // free
	VInfo->hp_refine   = calloc(NVglobal                , sizeof *(VInfo->hp_refine));   // free
	VInfo->hp_coarse   = calloc(NVglobal                , sizeof *(VInfo->hp_coarse));   // free
	VInfo->hp_coarse_l = calloc(NVglobal                , sizeof *(VInfo->hp_coarse_l)); // free

	VNeigh                  = VInfo->neigh;
	VType_global            = VInfo->type;
	vh_range                = VInfo->vh_range;
	vh_range                = VInfo->vh_range;
	fh_range                = VInfo->fh_range;
	p_levels                = VInfo->p_levels;
	h_levels                = VInfo->h_levels;
	h_siblings              = VInfo->h_siblings;
	h_forbid_coarse         = VInfo->h_forbid_c;
	hp_refine_current       = VInfo->hp_refine;
	hp_coarse_current       = VInfo->hp_coarse;
	hp_coarse_current_local = VInfo->hp_coarse_l;
	hp_coarse_current_err   = calloc(NVglobal , sizeof *hp_coarse_current_err); // free

	NFREFMAX_Total = NFMAX*NFREFMAX;

	for (i = 0; i < NVglobal; i++)
		h_forbid_coarse[i] = 1;

// Can remove initializations when the code is working (ToBeDeleted)
	PMaxP1      = PMax+1;
	LevelsMaxP1 = LevelsMax+1;
	for (i = 0; i < NVglobal; i++) {
		for (j = 0, jMax = NFREFMAX_Total; j < jMax; j++)
			VNeigh[i*jMax+j] = NVglobal;
		p_levels[i] = PMaxP1;
		h_levels[i] = LevelsMaxP1;
	}

	// Remove forbid for coarse_frac VOLUMES which satisfy tolerance condition
	for (i = 0, iMax = (unsigned int) (coarse_frac*NV); i < iMax; i++) {
		if (minRHS_Vec[i] < COARSE_TOL) {
			VOLUME = VOLUME_Vec[IndminRHS[i]];
			indexg = VOLUME->indexg;
			h_forbid_coarse[indexg] = 0;
		}
	}

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		indexg = VOLUME->indexg;
		ELEMENT = get_ELEMENT_type(VOLUME->type);

		Nf    = ELEMENT->Nf;
		for (f = 0; f < Nf; f++) {
			Indf = (indexg*NFMAX+f)*2;

			get_fh_range(VOLUME,f,&fhMin,&fhMax);
			fh_range[Indf  ] = fhMin;
			fh_range[Indf+1] = fhMax;
			for (fh = fhMin; fh <= fhMax; fh++) {
				Vf = f*NFREFMAX+fh;
				VNeigh[indexg*NFREFMAX_Total+Vf] = VOLUME->neigh[Vf];
			}
		}

		VType_global[indexg] = VOLUME->type;
		p_levels[indexg]     = VOLUME->P;
		h_levels[indexg]     = VOLUME->level;

		get_vh_range(VOLUME,&vhMin,&vhMax);
		vh_range[indexg*2  ] = vhMin;
		vh_range[indexg*2+1] = vhMax;

		VOLUMEp = VOLUME->parent;
		if (VOLUMEp) {
			VOLUMEc = VOLUMEp->child0;
			if (VOLUMEc->child0) {
				h_forbid_coarse[indexg] = 1;
			} else {
				for (i = 0, vh = vhMin; vh <= vhMax; vh++) {
					if (VOLUME->level != VOLUMEc->level || h_forbid_coarse[VOLUMEc->indexg]) {
						h_forbid_coarse[indexg] = 1;
						break;
					}
					h_siblings[indexg*NSIBMAX+i] = VOLUMEc->indexg;
					i++;
					VOLUMEc = VOLUMEc->next;
				}
			}
		} else {
			h_forbid_coarse[indexg] = 1;
		}
	}

	// Only VOLUMEs allowed to be coarsened.
	for (i = 0, iMax = (unsigned int) (coarse_frac*NV); i < iMax; i++) {
		if (minRHS_Vec[i] < COARSE_TOL) {
			VOLUME = VOLUME_Vec[IndminRHS[i]];
			indexg = VOLUME->indexg;
			if (Adapt == ADAPT_P || !h_forbid_coarse[indexg])
				hp_coarse_current_err[indexg] = 1;
		}
	}

	// Needs modification for MPI (ToBeDeleted)
	for (i = 0, iMax = NVglobal; i < iMax; i++) {
		if (VType_global[i] == 0) {
			printf("Error: Modifications needed for MPI in adapt_hp.\n"), EXIT_MSG;
			// Not all VOLUMEs are on this processor.
		}
	}

	// Mark refine_frac VOLUMEs for refinement
	for (i = 0, iMax = (unsigned int) (refine_frac_lim*NV); i < iMax; i++) {
		iInd = NV-i-1;
		if (maxRHS_Vec[iInd] > REFINE_TOL) {
			VOLUME = VOLUME_Vec[IndmaxRHS[iInd]];

			switch (Adapt) {
			default: // ADAPT_HP
				printf("Error: Code up the smoothness based indicator to choose h or p.\n"), EXIT_MSG;
				break;
			case ADAPT_P:
// No need to recompute VNeigh/VType_global for ADAPT_P at each time step (ToBeDeleted)
				check_levels_refine(VOLUME->indexg,VInfo,'p');
				break;
			case ADAPT_H:
				check_levels_refine(VOLUME->indexg,VInfo,'h');
				break;
			}
		}
	}

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		indexg = VOLUME->indexg;
		switch (Adapt) {
		default: // ADAPT_HP
			printf("Error: Fix case ADAPT_HP.\n"), EXIT_MSG;
			break;
		case ADAPT_P:
			if (hp_refine_current[indexg]) {
				VOLUME->adapt_type = PREFINE;
				if (p_levels[indexg] <= PMax)
					VOLUME->Vadapt = 1;
			}
			break;
		case ADAPT_H:
			if (hp_refine_current[indexg]) {
				VOLUME->adapt_type = HREFINE;
				if (h_levels[indexg] <= LevelsMax)
					VOLUME->Vadapt = 1;
// ToBeDeleted: Add in support for choosing which type of h-refinement to perform (isotropic (0)/several anisotropic)
// Perform isotropic refinement for elements refined due to refinement propagation in order to avoid unsupported
// match ups.
				VOLUME->hrefine_type = 0;
			}
			break;
		}
	}

	// Mark coarse_frac VOLUMEs for coarsening
	for (i = 0, iMax = (unsigned int) (coarse_frac*NV); i < iMax; i++) {
		if (minRHS_Vec[i] < COARSE_TOL) {
			VOLUME = VOLUME_Vec[IndminRHS[i]];

			for (j = 0; j < NVglobal; j++)
				hp_coarse_current_local[j] = 0;

			switch (Adapt) {
			default: // ADAPT_HP
				printf("Error: Code up the smoothness based indicator to choose h or p.\n"), EXIT_MSG;
				break;
			case ADAPT_P:
				check_levels_coarse(VOLUME->indexg,VInfo,'p');

				// Check for conflicts with elements flagged for refinement
				refine_conflict = 0;
				for (j = 0; j < NVglobal; j++) {
					if (hp_refine_current[j] && hp_coarse_current_local[j]) {
						refine_conflict = 1;
						break;
					}
				}
				if (!refine_conflict) {
					// Ensure that all elements to which the coarsening will propagate also have minRHS > COARSE_TOL
					coarse_conflict = 0;
					for (j = 0; j < NVglobal; j++) {
						if (hp_coarse_current_local[j] && minRHS_Vec_unsorted[j] > COARSE_TOL) {
							coarse_conflict = 1;
							break;
						}
					}

					if (!coarse_conflict) {
						for (j = 0; j < NVglobal; j++) {
							if (hp_coarse_current_local[j])
								hp_coarse_current[j] = 1;
						}
					}
				}
				break;
			case ADAPT_H:
				indexg = VOLUME->indexg;
				if (!h_forbid_coarse[indexg])
					check_levels_coarse(indexg,VInfo,'h');

				coarse_conflict = 0;
				for (j = 0; j < NVglobal; j++) {
					if (hp_coarse_current_local[j] && !hp_coarse_current_err[j]) {
						coarse_conflict = 1;
						break;
					}
				}
				if (!coarse_conflict) {
					for (j = 0; j < NVglobal; j++) {
						if (hp_coarse_current_local[j])
							hp_coarse_current[j] = 1;
					}
				}
				break;
			}
		}
	}

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		indexg = VOLUME->indexg;
		switch (Adapt) {
		default: // ADAPT_HP
			printf("Error: Fix case ADAPT_HP.\n"), EXIT_MSG;
			break;
		case ADAPT_P:
			if (hp_coarse_current[indexg]) {
				VOLUME->adapt_type = PCOARSE;
				if (p_levels[indexg] != 0)
					VOLUME->Vadapt = 1;
			}
			break;
		case ADAPT_H:
			if (hp_coarse_current[indexg]) {
				VOLUME->adapt_type = HCOARSE;
				if (h_levels[indexg] != 0)
					VOLUME->Vadapt = 1;
			}
			break;
		}
// ToBeDeleted
/*
if (hp_coarse_current_err[indexg]) {
	if (Adapt == ADAPT_P) {
		if (!VOLUME->adapt_type == PREFINE)
			VOLUME->adapt_type = PCOARSE;
	} else if (Adapt == ADAPT_H) {
		if (!VOLUME->adapt_type == HREFINE)
			VOLUME->adapt_type = HCOARSE;
	}
}
*/
	}

	free(minRHS_Vec);
	free(maxRHS_Vec);
	free(VOLUME_Vec);

	free(minRHS_Vec_unsorted);
	free(IndminRHS);
	free(IndmaxRHS);

	free(VInfo->neigh);
	free(VInfo->type);
	free(VInfo->vh_range);
	free(VInfo->fh_range);
	free(VInfo->href_type);
	free(VInfo->p_levels);
	free(VInfo->h_levels);
	free(VInfo->h_siblings);
	free(VInfo->h_forbid_c);
	free(VInfo->hp_refine);
	free(VInfo->hp_coarse);
	free(VInfo->hp_coarse_l);
	free(VInfo);

	free(hp_coarse_current_err);
}

void mesh_update(void)
{
	update_VOLUME_hp();
	update_FACE_hp();
	update_VOLUME_list();
	memory_free_children();
//	update_Vgrp();
	if (DB.Vectorized)
		printf("Error: update_Vgrp requires modifications when adaptation is enabled.\n"), EXIT_MSG;
	update_VOLUME_Ops();
	update_VOLUME_finalize();
}

void mesh_to_level(const unsigned int level)
{
	/*
	 *	Comments:
	 *		This function does not attempt to achieve the minimal L2 solution error while switching levels as it may
	 *		first project to level 0 even if the desired final level is above 0.
	 */

	unsigned int updated = 1, Vlevel, VlevelMax;
	struct S_VOLUME *VOLUME;

	if (level == 0) {
		while (updated) {
			updated = 0;
			VlevelMax = 0;
			for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
				Vlevel = VOLUME->level;
				if (Vlevel > VlevelMax)
					VlevelMax = Vlevel;
			}

			for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
				Vlevel = VOLUME->level;
				if (Vlevel != level && Vlevel == VlevelMax) {
					updated = 1;

					VOLUME->Vadapt = 1;
					VOLUME->adapt_type = HCOARSE;
				}
			}

			if (updated)
				mesh_update();
		}
	} else {
		// First project to level 0 if VOLUMEs are on different levels
		for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			if (VOLUME->level != DB.VOLUME->level) {
				mesh_to_level(0);
				break;
			}
		}

		while (updated) {
			updated = 0;
			for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
				if (VOLUME->level != level) {
					updated = 1;

					VOLUME->Vadapt = 1;
					if (VOLUME->level > level)
						VOLUME->adapt_type = HCOARSE;
					else
						VOLUME->adapt_type = HREFINE;
				}
			}

			if (updated)
				mesh_update();
		}
	}
}

void mesh_to_order(const unsigned int order)
{
	unsigned int updated = 1;
	struct S_VOLUME *VOLUME;

	while (updated) {
		updated = 0;

		for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			if (VOLUME->P != order) {
				updated = 1;

				VOLUME->Vadapt = 1;
				if (VOLUME->P > order)
					VOLUME->adapt_type = PCOARSE;
				else
					VOLUME->adapt_type = PREFINE;
			}
		}

		if (updated)
			mesh_update();
	}
}

void mesh_h_adapt(const unsigned int Nadapt, const char h_adapt_type)
{
	unsigned int n;
	struct S_VOLUME *VOLUME;

	for (n = 0; n < Nadapt; n++) {
		if (h_adapt_type == 'c') {
			for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
				VOLUME->Vadapt = 1;
				VOLUME->adapt_type = HCOARSE;
			}
		} else if (h_adapt_type == 'r') {
			for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
				VOLUME->Vadapt = 1;
				VOLUME->adapt_type = HREFINE;
			}
		}
		mesh_update();
	}
}

static bool check_conflict(const unsigned int indexg, const unsigned int AdaptType, const struct S_VInfo *VInfo)
{
	/*
	 *	Comments:
	 *		The conditions used to establish irregularity are necessary to allow h and p adaptation in neighbouring
	 *		elements while maintaining the 1-irregularity of the mesh. For example, taking the left element to have
	 *		[P,ML] = [3,2] and the right element [4,2], then a PREFINE in the left and HREFINE in the right element
	 *		would be acceptable but would have returned a conflict had the irregularity condition simply been set to
	 *		hp_levels[indexg]-hp_levels[indexg_neigh] != 0 because of the different AdaptType.
	 */

	unsigned int *hp_levels, *hp_update;

	hp_update = VInfo->hp_update;
	if (AdaptType == HREFINE || AdaptType == HCOARSE)
		hp_levels = VInfo->h_levels;
	else if (AdaptType == PREFINE || AdaptType == PCOARSE)
		hp_levels = VInfo->p_levels;
	else
		printf("Error: Unsupported.\n"), EXIT_MSG;

	unsigned int f, Nf, Indf, indexg_neigh;

	struct S_ELEMENT *ELEMENT;

	ELEMENT = get_ELEMENT_type(VInfo->type[indexg]);
	Nf = ELEMENT->Nf;

	Indf = indexg*NFREFMAX*NFMAX;
	if (AdaptType == HREFINE || AdaptType == PREFINE || AdaptType == PCOARSE) {
		unsigned int fh, fhMin, fhMax, Indfh, *fh_range;
		fh_range = VInfo->fh_range;
		for (f = 0; f < Nf; f++) {
			Indfh = Indf + f*NFREFMAX;
			fhMin = fh_range[(indexg*NFMAX+f)*2  ];
			fhMax = fh_range[(indexg*NFMAX+f)*2+1];
			for (fh = fhMin; fh <= fhMax; fh++) {
				indexg_neigh = VInfo->neigh[Indfh+fh];

				bool irregular;
				if (AdaptType == HREFINE || AdaptType == PREFINE)
					irregular = ((int) hp_levels[indexg] - (int) hp_levels[indexg_neigh]) > 0;
				else
					irregular = ((int) hp_levels[indexg] - (int) hp_levels[indexg_neigh]) < 0;

				if (irregular) {
					if (!hp_update[indexg_neigh]) {
						if (check_conflict(indexg_neigh,AdaptType,VInfo))
							return 1;
					} else {
						// Problem if neighbour is flagged to be updated with a different type of adaptation.
						if (hp_update[indexg_neigh] != AdaptType)
							return 1;
					}
				}
			}
		}
	} else {
		printf("Add support for HCOARSE.\n"), EXIT_MSG;
	}
	return 0;
}

static void update_levels(const unsigned int indexg, const unsigned AdaptType, const struct S_VInfo *VInfo)
{
	/*
	 *	Comments:
	 *		Remove redundant code when working (ToBeDeleted).
	 */

	bool         conflict;
	unsigned int *hp_update, *hp_levels;

	if (AdaptType == HREFINE)
		conflict = 0;
	else
		// This call is redundant after the first for the element under consideration. (ToBeModified)
		conflict = check_conflict(indexg,AdaptType,VInfo);

	hp_update = VInfo->hp_update;
	if (conflict) {
		hp_update[indexg] = 0;
		return;
	}

	hp_update[indexg] = AdaptType;
	if (AdaptType == HREFINE || AdaptType == HCOARSE)
		hp_levels = VInfo->h_levels;
	else if (AdaptType == PREFINE || AdaptType == PCOARSE)
		hp_levels = VInfo->p_levels;
	else
		printf("Error: Unsupported.\n"), EXIT_MSG;

	unsigned int f, Nf, Indf, indexg_neigh;

	struct S_ELEMENT *ELEMENT;

	ELEMENT = get_ELEMENT_type(VInfo->type[indexg]);
	Nf = ELEMENT->Nf;

	Indf = indexg*NFREFMAX*NFMAX;
	if (AdaptType == HREFINE || AdaptType == PREFINE || AdaptType == PCOARSE) {
		unsigned int fh, fhMin, fhMax, Indfh, *fh_range;
		fh_range = VInfo->fh_range;
		for (f = 0; f < Nf; f++) {
			Indfh = Indf + f*NFREFMAX;
			fhMin = fh_range[(indexg*NFMAX+f)*2  ];
			fhMax = fh_range[(indexg*NFMAX+f)*2+1];
			for (fh = fhMin; fh <= fhMax; fh++) {
				indexg_neigh = VInfo->neigh[Indfh+fh];

				bool irregular;
				if (AdaptType == HREFINE || AdaptType == PREFINE)
					irregular = ((int) hp_levels[indexg] - (int) hp_levels[indexg_neigh]) > 0;
				else
					irregular = ((int) hp_levels[indexg] - (int) hp_levels[indexg_neigh]) < 0;

				if (!hp_update[indexg_neigh] && irregular)
					update_levels(indexg_neigh,AdaptType,VInfo);
			}
		}
	} else {
		printf("Add support for HCOARSE.\n"), EXIT_MSG;
	}
}

void ensure_1irregular(unsigned int *hp_update)
{
	/*
	 *	Purpose:
	 *		Update hp_update to ensure that the mesh resulting from the refinement/coarsening as specified in hp_update
	 *		is not more than 1-irregular in h and p.
	 *
	 *	Comments:
	 *		Marking VOLUMEs for h and p adaptation is performed sequentially to ensure that only one option is provided
	 *		for each VOLUME. The order of precedence is: HREFINE, PREFINE, HCOARSE, PCOARSE.
	 *		To enable parallel treatment (MPI) reliance on local structures (VOLUME, FACE) is avoided.
	 *		Must be updated for MPI support.
	 */

	// Initialize DB Parameters
	unsigned int NV       = DB.NV,
	             NVglobal = DB.NVglobal,
	             Adapt    = DB.Adapt;

	// Standard datatypes
	unsigned int v, NFREFMAX_Total;

	struct S_VOLUME  *VOLUME, **VOLUME_Vec;

	NFREFMAX_Total = NFMAX*NFREFMAX;

	VOLUME_Vec = malloc(NV * sizeof *VOLUME_Vec); // free

	for (VOLUME = DB.VOLUME, v = 0; VOLUME; VOLUME = VOLUME->next, v++) {
		VOLUME_Vec[v] = VOLUME;
	}

	VInfo = malloc(sizeof *VInfo); // free
	VInfo->hp_update = hp_update;

	VInfo->neigh       = malloc(NVglobal*NFMAX*NFREFMAX * sizeof *(VInfo->neigh));       // free
	VInfo->type        = calloc(NVglobal                , sizeof *(VInfo->type));        // free
	VInfo->fh_range    = malloc(NVglobal*NFMAX*2        * sizeof *(VInfo->fh_range));    // free
	VInfo->p_levels    = malloc(NVglobal                * sizeof *(VInfo->p_levels));    // free
	VInfo->h_levels    = malloc(NVglobal                * sizeof *(VInfo->h_levels));    // free

	unsigned int *VNeigh, *VType_global, *fh_range, *p_levels, *h_levels;

	VNeigh                  = VInfo->neigh;
	VType_global            = VInfo->type;
	fh_range                = VInfo->fh_range;
	p_levels                = VInfo->p_levels;
	h_levels                = VInfo->h_levels;
	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		unsigned int indexg, f, Nf;

		struct S_ELEMENT *ELEMENT;

		indexg = VOLUME->indexg;
		ELEMENT = get_ELEMENT_type(VOLUME->type);

		Nf = ELEMENT->Nf;
		for (f = 0; f < Nf; f++) {
			unsigned int Indf, fh, fhMin, fhMax, Vf;

			Indf = (indexg*NFMAX+f)*2;

			get_fh_range(VOLUME,f,&fhMin,&fhMax);
			fh_range[Indf  ] = fhMin;
			fh_range[Indf+1] = fhMax;
			for (fh = fhMin; fh <= fhMax; fh++) {
				Vf = f*NFREFMAX+fh;
				VNeigh[indexg*NFREFMAX_Total+Vf] = VOLUME->neigh[Vf];
			}
		}

		VType_global[indexg] = VOLUME->type;
		p_levels[indexg]     = VOLUME->P;
		h_levels[indexg]     = VOLUME->level;
	}

	// Needs modification for MPI (ToBeDeleted)
	for (v = 0; v < NVglobal; v++) {
		if (VType_global[v] == 0)
			printf("Error: Add support.\n"), EXIT_MSG; // Not all VOLUMEs are on this processor.
	}

	// HREFINE
	if (Adapt == ADAPT_H || Adapt == ADAPT_HP) {
		for (v = 0; v < NV; v++) {
			VOLUME = VOLUME_Vec[v];
			if (hp_update[v] == HREFINE)
				update_levels(VOLUME->indexg,HREFINE,VInfo);
		}
	}

	// PREFINE
	if (Adapt == ADAPT_P || Adapt == ADAPT_HP) {
		for (v = 0; v < NV; v++) {
			VOLUME = VOLUME_Vec[v];
			if (hp_update[v] == PREFINE)
				update_levels(VOLUME->indexg,PREFINE,VInfo);
		}
	}

	// HCOARSE
	if (Adapt == ADAPT_H || Adapt == ADAPT_HP) {
		for (v = 0; v < NV; v++) {
			VOLUME = VOLUME_Vec[v];
			if (hp_update[v] == HCOARSE)
				update_levels(VOLUME->indexg,HCOARSE,VInfo);
		}
	}

	// PCOARSE
	if (Adapt == ADAPT_P || Adapt == ADAPT_HP) {
		for (v = 0; v < NV; v++) {
			VOLUME = VOLUME_Vec[v];
			if (hp_update[v] == PCOARSE)
				update_levels(VOLUME->indexg,PCOARSE,VInfo);
		}
	}

	free(VOLUME_Vec);
	free(VNeigh);
	free(VType_global);
	free(fh_range);
	free(p_levels);
	free(h_levels);
	free(VInfo);
}
