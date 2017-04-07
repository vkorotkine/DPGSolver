// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "explicit_VOLUME_info_c.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"

#include "element_functions.h"
#include "solver_functions.h"
#include "solver_functions_c.h"
#include "matrix_functions.h"
#include "fluxes_inviscid_c.h"
#include "array_print.h"

/*
 *	Purpose:
 *		Identical to explicit_VOLUME_info using complex variables (for complex step verification).
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

static void compute_VOLUME_RHS_EFE(void);

void explicit_VOLUME_info_c(void)
{
	// Initialize DB Parameters
	unsigned int EFE        = DB.EFE,
	             Vectorized = DB.Vectorized;

	if (EFE) {
		switch (Vectorized) {
		case 0:
			compute_VOLUME_RHS_EFE();
			break;
		default:
			EXIT_UNSUPPORTED;
//			compute_VOLUMEVec_RHS_EFE();
			break;
		}
	} else {
		;
	}
}

static void compute_VOLUME_RHS_EFE(void)
{
	// Initialize DB Parameters
	char         *Form = DB.Form;
	unsigned int d          = DB.d,
	             Collocated = DB.Collocated,
				 Nvar       = DB.Nvar,
				 Neq        = DB.Neq;

	// Standard datatypes
	struct S_OPERATORS_V *OPS[2];
	struct S_VDATA       *VDATA;

	VDATA = malloc(sizeof *VDATA); // free
	VDATA->OPS = OPS;

	for (size_t i = 0; i < 2; i++)
		OPS[i] = malloc(sizeof *OPS[i]); // free

	if (strstr(Form,"Weak")) {
		for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			init_VDATA(VDATA,VOLUME);

			// Obtain W_vI
			unsigned int NvnI = VDATA->OPS[0]->NvnI;
			if (Collocated) {
				VDATA->W_vI_c = VOLUME->What_c;
			} else {
				VDATA->W_vI_c = malloc(NvnI*Nvar * sizeof *(VDATA->W_vI_c)); // free
				coef_to_values_vI_c(VDATA,'W');
			}

			// Compute Flux in reference space
			double complex *F_vI = malloc(NvnI*d*Neq * sizeof *F_vI); // free

			flux_inviscid_c(NvnI,1,VDATA->W_vI_c,F_vI,d,Neq);

			if (!Collocated)
				free(VDATA->W_vI_c);

			// Convert to reference space
			double complex *Fr_vI = malloc(NvnI*d*Neq * sizeof *Fr_vI); // free
			convert_between_rp_c(NvnI,Neq,VOLUME->C_vI,F_vI,Fr_vI,"FluxToRef");
			free(F_vI);

			// RHS
			unsigned int NvnS = VDATA->OPS[0]->NvnS;

			if (VOLUME->RHS_c)
				free(VOLUME->RHS_c);
			double complex *RHS = calloc(NvnS*Neq , sizeof *RHS); // keep
			VOLUME->RHS_c = RHS;

			finalize_VOLUME_Inviscid_Weak_c(Neq,Fr_vI,RHS,'E',VDATA);
			free(Fr_vI);
		}
	} else if (strstr(Form,"Strong")) {
		EXIT_UNSUPPORTED;
	}

	free(VDATA);
	for (size_t i = 0; i < 2; i++)
		free(OPS[i]);
}
