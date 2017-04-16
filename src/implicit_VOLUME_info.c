// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "implicit_VOLUME_info.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_VOLUME.h"

#include "solver_functions.h"
#include "fluxes_inviscid.h"
#include "jacobian_fluxes_inviscid.h"
#include "array_print.h"

/*
 *	Purpose:
 *		Evaluate the VOLUME contributions to the RHS and LHS terms.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

static void compute_VOLUME_EFE(void);

void implicit_VOLUME_info(void)
{
	// Initialize DB Parameters
	unsigned int EFE        = DB.EFE,
	             Vectorized = DB.Vectorized;

	if (EFE) {
		switch (Vectorized) {
		case 0:
			compute_VOLUME_EFE();
			break;
		default:
			printf("Error: Unsupported Vectorized (%d).\n",Vectorized), EXIT_MSG;
			break;
		}
	} else {
		;
	}
}

static void compute_VOLUME_EFE(void)
{
	// Initialize DB Parameters
	char         *Form = DB.Form;
	unsigned int d          = DB.d,
	             Collocated = DB.Collocated,
				 Nvar       = DB.Nvar,
				 Neq        = DB.Neq;

	struct S_OPERATORS_V *OPS[2];
	struct S_VDATA       *VDATA;

	VDATA = malloc(sizeof *VDATA); // free
	VDATA->OPS = (struct S_OPERATORS_V const *const *) OPS;

	for (size_t i = 0; i < 2; i++)
		OPS[i] = malloc(sizeof *OPS[i]); // free

	if (strstr(Form,"Weak")) {
		for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			init_VDATA(VDATA,VOLUME);

			// Obtain W_vI
			unsigned int NvnI = VDATA->OPS[0]->NvnI;
			if (Collocated) {
				VDATA->W_vI = VOLUME->What;
			} else {
				VDATA->W_vI = malloc(NvnI*Nvar * sizeof *(VDATA->W_vI)); // free
				coef_to_values_vI(VDATA,'W');
			}

			// Compute Flux and its Jacobian in reference space
			double *F_vI    = malloc(NvnI*d*Neq      * sizeof *F_vI);    // free
			double *dFdW_vI = malloc(NvnI*d*Nvar*Neq * sizeof *dFdW_vI); // free

			flux_inviscid(NvnI,1,VDATA->W_vI,F_vI,d,Neq);
			jacobian_flux_inviscid(NvnI,1,VDATA->W_vI,dFdW_vI,d,Neq);

			if (!Collocated)
				free(VDATA->W_vI);

			// Convert to reference space
			double *Fr_vI = malloc(NvnI*d*Neq * sizeof *Fr_vI); // free
			convert_between_rp(NvnI,Neq,VOLUME->C_vI,F_vI,Fr_vI,"FluxToRef");
			free(F_vI);

			double *dFrdW_vI = malloc(NvnI*d*Nvar*Neq * sizeof *dFrdW_vI); // free
			convert_between_rp(NvnI,Nvar*Neq,VOLUME->C_vI,dFdW_vI,dFrdW_vI,"FluxToRef");
			free(dFdW_vI);

			// Compute RHS and LHS terms
			unsigned int NvnS = VDATA->OPS[0]->NvnS;

			// RHS
			memset(VOLUME->RHS,0.0,NvnS*Neq * sizeof *(VOLUME->RHS));
			finalize_VOLUME_Inviscid_Weak(Neq,Fr_vI,VOLUME->RHS,'E',VDATA);
			free(Fr_vI);

			// LHS
			memset(VOLUME->LHS,0.0,NvnS*NvnS*Neq*Nvar * sizeof *(VOLUME->LHS));
			finalize_VOLUME_Inviscid_Weak(NvnS*Neq*Nvar,dFrdW_vI,VOLUME->LHS,'I',VDATA);
			free(dFrdW_vI);
		}
	} else if (strstr(Form,"Strong")) {
		EXIT_UNSUPPORTED;
	}

	free(VDATA);
	for (size_t i = 0; i < 2; i++)
		free(OPS[i]);
}
