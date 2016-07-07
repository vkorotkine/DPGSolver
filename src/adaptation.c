// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

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

void get_PS_range(unsigned int *PSMin, unsigned int *PSMax)
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

void get_Pb_range(const unsigned int P, unsigned int *PbMin, unsigned int *PbMax)
{
	/*
	 *	Purpose:
	 *		Return the range of orders used for operators which interpolate between different orders.
	 */

	// Initialize DB Parameters
	unsigned int PMax  = DB.PMax,
	             Adapt = DB.Adapt;

	switch (Adapt) {
		default: // ADAPT_P or ADAPT_HP
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

static void check_levels_refine(const unsigned int indexg, const unsigned int *VNeigh, const unsigned int *VType,
                                unsigned int *hp_levels, unsigned int *hp_refine_current, const char hp_type)
{
	// Standard datatypes
	unsigned int Nf, *Nfref, indexg_neigh, f, fh, fhMax;

	struct S_ELEMENT *ELEMENT;

	ELEMENT = get_ELEMENT_type(VType[indexg]);
	Nf = ELEMENT->Nf;
	Nfref = ELEMENT->Nfref;

	if (!hp_refine_current[indexg]) {
		hp_refine_current[indexg] = 1;

		switch (hp_type) {
		default: // p
			hp_levels[indexg]++;
			for (f = 0; f < Nf; f++) {
				indexg_neigh = VNeigh[indexg*NFREFMAX*NFMAX+f*NFREFMAX];
//if (indexg == 24)
//	printf("indexg_neigh: %d\n",indexg_neigh);
				if (!hp_refine_current[indexg_neigh] && ((int) hp_levels[indexg] - (int) hp_levels[indexg_neigh]) > 1)
					check_levels_refine(indexg_neigh,VNeigh,VType,hp_levels,hp_refine_current,'p');
			}
			break;
		case 'h':
			for (f = 0; f < Nf; f++) {
			for (fh = 0, fhMax = Nfref[f]; fh < fhMax; fh++) {
			}}
			printf("Error: Write check_levels_refine for hp_type = 'h'.\n"), exit(1);
			break;
		}
	}
}

static void check_levels_coarse(const unsigned int indexg, const unsigned int *VNeigh, const unsigned int *VType,
                                const unsigned int *hp_levels, unsigned int *hp_coarse_current, const char hp_type)
{
	// Standard datatypes
	unsigned int Nf, *Nfref, indexg_neigh, f, fh, fhMax;

	struct S_ELEMENT *ELEMENT;

	ELEMENT = get_ELEMENT_type(VType[indexg]);
	Nf    = ELEMENT->Nf;
	Nfref = ELEMENT->Nfref;

	hp_coarse_current[indexg] = 1;

	switch (hp_type) {
	default: // p
		for (f = 0; f < Nf; f++) {
			indexg_neigh = VNeigh[indexg*NFREFMAX*NFMAX+f*NFREFMAX];
			if (!hp_coarse_current[indexg_neigh] && ((int) hp_levels[indexg_neigh] - (int) hp_levels[indexg]) > 0)
				check_levels_coarse(indexg_neigh,VNeigh,VType,hp_levels,hp_coarse_current,'p');
		}
		break;
	case 'h':
		for (f = 0; f < Nf; f++) {
		for (fh = 0, fhMax = Nfref[f]; fh < fhMax; fh++) {
		}}
		printf("Error: Write check_levels_refine for hp_type = 'h'.\n"), exit(1);
		break;
	}
}

void adapt_hp(void)
{
	// Initialize DB Parameters
	unsigned int DOF0     = DB.DOF0,
	             NV       = DB.NV,
				 NVglobal = DB.NVglobal,
				 PMax     = DB.PMax,
	             Adapt    = DB.Adapt;
	double       refine_frac = DB.refine_frac,
	             coarse_frac = DB.coarse_frac,
	             DOFcap_frac = DB.DOFcap_frac;

	// Standard datatypes
	unsigned int i, iMax, j, jMax, iInd, DOF, NvnS, indexg, PMaxP1, Nf, *Nfref, f, fh, Vf, fhMax,
	             NFREFMAX_Total, refine_conflict, coarse_conflict,
	             *IndminRHS, *IndmaxRHS, *VNeigh, *VType_global,
	             *p_levels, *h_levels, *hp_refine_current, *hp_coarse_current, *hp_coarse_current_local;
	double       minRHS, maxRHS, tmp_d, refine_frac_lim,
	             *RHS, *minRHS_Vec, *maxRHS_Vec, *minRHS_Vec_unsorted;

	struct S_ELEMENT *ELEMENT;
	struct S_VOLUME  *VOLUME, **VOLUME_Vec;

	// silence
	indexg = NVglobal;
	refine_frac_lim = 0.0;

	minRHS_Vec = malloc(NV * sizeof *minRHS_Vec); // free
	maxRHS_Vec = malloc(NV * sizeof *maxRHS_Vec); // free
	VOLUME_Vec = malloc(NV * sizeof *VOLUME_Vec); // free

// Change this from i to indexg (for MPI). Also then initialize min/maxRHS_Vec so that they can be sorted correctly.
	i = 0; DOF = 0;
	for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
		// Compute maxRHS in each VOLUME
		NvnS = VOLUME->NvnS;
		DOF += NvnS;

		RHS = VOLUME->RHS;

		minRHS = 1e10;
		maxRHS = 0.0;
		// Density RHS
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
printf("ref_frac_lim: %f\n",refine_frac_lim);
	if (refine_frac_lim < EPS) {
		printf("*** Warning: Consider raising DOFcap_frac. *** \n");
	}

	VNeigh            = malloc(NVglobal*NFMAX*NFREFMAX * sizeof *VNeigh); // free

	VType_global      = calloc(NVglobal , sizeof *VType_global);      // free
	p_levels          = malloc(NVglobal * sizeof *p_levels);          // free
	h_levels          = malloc(NVglobal * sizeof *h_levels);          // free
	hp_refine_current = calloc(NVglobal , sizeof *hp_refine_current); // free
	hp_coarse_current = calloc(NVglobal , sizeof *hp_coarse_current); // free
	
	NFREFMAX_Total = NFMAX*NFREFMAX;

	for (i = 0, PMaxP1 = PMax+1; i < NVglobal; i++) {
		for (j = 0, jMax = NFREFMAX_Total; j < jMax; j++)
			VNeigh[i*jMax+j] = NVglobal;
		p_levels[i] = PMaxP1;
		h_levels[i] = 100; // Can initialize to h_levels to hlevelsMax+1 in future (ToBeDeleted)
	}

	for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
		indexg = VOLUME->indexg;
		ELEMENT = get_ELEMENT_type(VOLUME->type);

		Nf    = ELEMENT->Nf;
		Nfref = ELEMENT->Nfref;
		for (f = 0; f < Nf; f++) {
		for (fh = 0, fhMax = Nfref[f]; fh < fhMax; fh++) {
			Vf = f*NFREFMAX+fh;
			VNeigh[indexg*NFREFMAX_Total+Vf] = VOLUME->neigh[Vf];
		}}

		VType_global[indexg] = VOLUME->type;
		p_levels[indexg]     = VOLUME->P;
		h_levels[indexg]     = VOLUME->hlevel;
	}
	// Needs modification for MPI (ToBeDeleted)
	for (i = 0, iMax = NVglobal; i < iMax; i++) {
		if (VType_global[i] == 0) {
			printf("Error: Modifications needed for MPI in adapt_hp.\n"), exit(1);
			// Not all VOLUMEs are on this processor.
		}
	}

//array_print_ui(NVglobal,1,p_levels,'R');
//array_print_ui(NVglobal,NFREFMAX_Total,VNeigh,'R');

	// Mark refine_frac VOLUMEs for refinement
//printf("\n\n\nrefinement marking\n\n");
	for (i = 0, iMax = (unsigned int) (refine_frac_lim*NV); i < iMax; i++) {
		iInd = NV-i-1;
//printf("i, indexg, maxRHS: %d, %d, % .3e\n",i,VOLUME_Vec[IndmaxRHS[iInd]]->indexg,maxRHS_Vec[iInd]);
		if (maxRHS_Vec[iInd] > REFINE_TOL) {
			VOLUME = VOLUME_Vec[IndmaxRHS[iInd]];

			switch (Adapt) {
			default: // ADAPT_HP
				printf("Error: Code up the smoothness based indicator to choose h or p.\n"), exit(1);
				break;
			case ADAPT_P:
// No need to recompute VNeigh/VType_global for ADAPT_P
				check_levels_refine(VOLUME->indexg,VNeigh,VType_global,p_levels,hp_refine_current,'p');
/*
for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next){
	printf("%d %d\n",VOLUME->indexg,hp_refine_current[VOLUME->indexg]);
}
printf("\n\n\n");
*/
				break;
			case ADAPT_H:
				check_levels_refine(VOLUME->indexg,VNeigh,VType_global,h_levels,hp_refine_current,'h');
for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next){
	printf("%d %d\n",VOLUME->indexg,hp_refine_current[VOLUME->indexg]);
}
				break;
			}
		}
	}

	for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
		indexg = VOLUME->indexg;
		switch (Adapt) {
		default: // ADAPT_HP
			printf("Error: Fix case ADAPT_HP.\n"), exit(1);
			break;
		case ADAPT_P:
			if (hp_refine_current[indexg] && p_levels[indexg] <= PMax) {
				VOLUME->Vadapt = 1;
				VOLUME->adapt_type = PREFINE;
			}
			break;
		case ADAPT_H:
// ToBeDeleted: Add condition for level Max
			if (hp_refine_current[indexg]) {
				VOLUME->Vadapt = 1;
				VOLUME->adapt_type = HREFINE;
			}
			break;
		}
	}

	// Mark coarse_frac VOLUMEs for coarsening
	for (i = 0, iMax = (unsigned int) (coarse_frac*NV); i < iMax; i++) {
		if (minRHS_Vec[i] < COARSE_TOL) {
			VOLUME = VOLUME_Vec[IndminRHS[i]];
//printf("indexg, minRHS: %d, % .3e\n",VOLUME->indexg,minRHS_Vec[i]);

			hp_coarse_current_local = calloc(NVglobal , sizeof *hp_coarse_current_local); // free

			switch (Adapt) {
			default: // ADAPT_HP
				printf("Error: Code up the smoothness based indicator to choose h or p.\n"), exit(1);
				break;
			case ADAPT_P:
				check_levels_coarse(VOLUME->indexg,VNeigh,VType_global,p_levels,hp_coarse_current_local,'p');
/*
for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next){
	printf("%d %d %d\n",VOLUME->indexg,hp_coarse_current_local[VOLUME->indexg],p_levels[VOLUME->indexg]);
}
printf("\n\n\n");
exit(1);
*/

				break;
			case ADAPT_H:
				check_levels_coarse(VOLUME->indexg,VNeigh,VType_global,h_levels,hp_coarse_current_local,'h');
				break;
			}

			// Check for conflicts with elements flagged for refinement
			refine_conflict = 0;
			for (j = 0; j < NVglobal; j++) {
				if (hp_refine_current[j] && hp_coarse_current_local[j]) {
					refine_conflict = 1;
//printf("refine conflict found!\n");
					break;
				}
			}
			if (!refine_conflict) {
				// Ensure that all elements to which the coarsening will propagate also have minRHS < COARSE_TOL
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
			free(hp_coarse_current_local);
		}
	}

	for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
		indexg = VOLUME->indexg;
		switch (Adapt) {
		default: // ADAPT_HP
			printf("Error: Fix case ADAPT_HP.\n"), exit(1);
			break;
		case ADAPT_P:
			if (hp_coarse_current[indexg] && p_levels[indexg] != 0) {
				VOLUME->Vadapt = 1;
				VOLUME->adapt_type = PCOARSE;
			}
			break;
		case ADAPT_H:
			if (hp_coarse_current[indexg] && h_levels[indexg] != 0) {
				VOLUME->Vadapt = 1;
				VOLUME->adapt_type = HCOARSE;
			}
			break;
		}
	}
/*
for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next){
	if (VOLUME->Vadapt)
		printf("\t");
	printf("indexg, P, Vadapt, adapt_type: %d %d %d %d\n",VOLUME->indexg,VOLUME->P,VOLUME->Vadapt,VOLUME->adapt_type);
}
//array_print_ui(NVglobal,1,hp_refine_current,'R');
*/

// ToBeDeleted: Note that there is no conflict after the max DOF is reached as refinement is not attempted and all
//elements are marked as if they need not be refined. Probably add a flag indicating whether elements wanted to be refined
//but could not be.

//exit(1);

	free(minRHS_Vec);
	free(maxRHS_Vec);
	free(VOLUME_Vec);

	free(minRHS_Vec_unsorted);
	free(IndminRHS);
	free(IndmaxRHS);

	free(VNeigh);
	free(VType_global);
	free(p_levels);
	free(h_levels);
	free(hp_refine_current);
	free(hp_coarse_current);
}
