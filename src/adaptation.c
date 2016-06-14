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

void adapt_hp(void)
{
	// Initialize DB Parameters
	unsigned int Nvar  = DB.Nvar,
	             DOF0  = DB.DOF0,
	             Adapt = DB.Adapt;
	double       refine_frac = DB.refine_frac,
	             coarse_frac = DB.coarse_frac,
	             DOFcap_frac = DB.DOFcap_frac;

	// Standard datatypes
	unsigned int i, iMax, NV, tmp_ui, DOF, NvnS,
	             *IndminRES, *IndmaxRES;
	double       minRES, maxRES, tmp_d,
	             *RES, *minRES_Vec, *maxRES_Vec;

	struct S_VOLUME *VOLUME, **VOLUME_Vec;

	NV = 0;
	for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next)
		NV++;

	minRES_Vec = malloc(NV * sizeof *minRES_Vec); // free
	maxRES_Vec = malloc(NV * sizeof *maxRES_Vec); // free
	VOLUME_Vec = malloc(NV * sizeof *VOLUME_Vec); // free

	i = 0; DOF = 0;
	for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
		// Compute maxRES in each VOLUME
		NvnS = VOLUME->NvnS;
		DOF += NvnS;

		RES = VOLUME->RES;

		minRES = 1e10;
		maxRES = 0.0;
		for (iMax = NvnS*Nvar; iMax--; ) {
			tmp_d = fabs(*RES);
			if (tmp_d < minRES)
				minRES = tmp_d;
			if (tmp_d > maxRES)
				maxRES = tmp_d;
			RES++;
		}
		minRES_Vec[i] = minRES;
		maxRES_Vec[i] = maxRES;
		VOLUME_Vec[i] = VOLUME;
		i++;
	}

	IndminRES = malloc(NV * sizeof *IndminRES); // free
	IndmaxRES = malloc(NV * sizeof *IndmaxRES); // free
	for (i = 0; i < NV; i++) {
		IndminRES[i] = i;
		IndmaxRES[i] = i;
	}

	array_sort_d(1,NV,minRES_Vec,IndminRES,'R','N');
	array_sort_d(1,NV,maxRES_Vec,IndmaxRES,'R','N');
	// reverse IndmaxRES
	for (i = 0, iMax = NV/2; i < iMax; i++) {
// Write an array_reverse function if this is used again (ToBeDeleted)
		tmp_ui            = IndmaxRES[i];
		IndmaxRES[i]      = IndmaxRES[NV-i-1];
		IndmaxRES[NV-i-1] = tmp_ui;
	}

// Add in MPI communication here for globally sorted array. (ToBeDeleted)
	// Make sure that the DOF cap is not exceeded
	if (refine_frac*DOF > DOFcap_frac*DOF0) {
		printf("*** Warning: Consider raising DOFcap_frac.\n ***");
		// refine_frac = (DOFcap_frac*DOF0)/DOF;
		/* enabling this may result in endless loop while not enabling may result in infinite refinement */
	}

	// Mark refine_frac VOLUMEs for refinement
	for (i = 0, iMax = (unsigned int) refine_frac*NV; i < iMax; i++) {
		VOLUME_Vec[IndmaxRES[i]]->Vadapt = 1;
		switch (Adapt) {
		default: // ADAPT_HP
			printf("Error: Code up the smoothness based indicator to choose h or p.\n"), exit(1);
			break;
		case ADAPT_P:
			VOLUME_Vec[IndmaxRES[i]]->adapt_type = PREFINE;
			break;
		case ADAPT_H:
			VOLUME_Vec[IndmaxRES[i]]->adapt_type = HREFINE;
			break;
		}
	}

	// Mark coarse_frac VOLUMEs for coarsening
	for (i = 0, iMax = (unsigned int) coarse_frac*NV; i < iMax; i++) {
		VOLUME_Vec[IndminRES[i]]->Vadapt = 1;
		switch (Adapt) {
		default: // ADAPT_HP
			printf("Error: Code up the smoothness based indicator to choose h or p.\n"), exit(1);
			break;
		case ADAPT_P:
			VOLUME_Vec[IndminRES[i]]->adapt_type = PCOARSE;
			break;
		case ADAPT_H:
			VOLUME_Vec[IndminRES[i]]->adapt_type = HCOARSE;
			break;
		}
	}

	free(minRES_Vec);
	free(maxRES_Vec);
	free(VOLUME_Vec);

	free(IndminRES);
	free(IndmaxRES);
}
