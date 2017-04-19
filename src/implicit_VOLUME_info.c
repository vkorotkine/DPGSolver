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
#include "fluxes_viscous.h"
#include "jacobian_fluxes_inviscid.h"
#include "jacobian_fluxes_viscous.h"
#include "array_free.h"
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

static void compute_Inviscid_VOLUME_EFE(void);
static void compute_Viscous_VOLUME_EFE(void);

void implicit_VOLUME_info(void)
{
	unsigned int EFE        = DB.EFE,
	             Vectorized = DB.Vectorized;

	if (EFE) {
		switch (Vectorized) {
		case 0:
			compute_Inviscid_VOLUME_EFE();
			compute_Viscous_VOLUME_EFE();
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
		}
	} else {
		;
	}
}

static void compute_Inviscid_VOLUME_EFE(void)
{
	unsigned int d    = DB.d,
				 Nvar = d+2,
				 Neq  = d+2;

	struct S_OPERATORS_V *OPS[2];
	struct S_VDATA       *VDATA;

	VDATA = malloc(sizeof *VDATA); // free
	VDATA->OPS = (struct S_OPERATORS_V const *const *) OPS;

	for (size_t i = 0; i < 2; i++)
		OPS[i] = malloc(sizeof *OPS[i]); // free

	if (strstr(DB.Form,"Weak")) {
		for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			init_VDATA(VDATA,VOLUME);

			// Obtain W_vI
			unsigned int NvnI = VDATA->OPS[0]->NvnI;
			if (DB.Collocated) {
				VDATA->W_vI = VOLUME->What;
			} else {
				VDATA->W_vI = malloc(NvnI*Nvar * sizeof *(VDATA->W_vI)); // free
				coef_to_values_vI(VDATA,'W');
			}

			// Compute Flux and its Jacobian in reference space
			double *const F_vI    = malloc(NvnI*d*Neq      * sizeof *F_vI),    // free
			       *const dFdW_vI = malloc(NvnI*d*Nvar*Neq * sizeof *dFdW_vI); // free

// Consider combining flux_inviscid and jacobian_flux_inviscid (ToBeDeleted)
			flux_inviscid(NvnI,1,VDATA->W_vI,F_vI,d,Neq);
			jacobian_flux_inviscid(NvnI,1,VDATA->W_vI,dFdW_vI,d,Neq);

			if (!DB.Collocated)
				free(VDATA->W_vI);

			// Convert to reference space
			double *const Fr_vI = malloc(NvnI*d*Neq * sizeof *Fr_vI); // free
			convert_between_rp(NvnI,Neq,VOLUME->C_vI,F_vI,Fr_vI,"FluxToRef");
			free(F_vI);

			double *const dFrdW_vI = malloc(NvnI*d*Nvar*Neq * sizeof *dFrdW_vI); // free
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
	} else if (strstr(DB.Form,"Strong")) {
		EXIT_UNSUPPORTED;
	}

	free(VDATA);
	for (size_t i = 0; i < 2; i++)
		free(OPS[i]);
}

static void compute_Viscous_VOLUME_EFE(void)
{
	unsigned int d    = DB.d,
				 Nvar = d+2,
				 Neq  = d+2;

	struct S_OPERATORS_V *OPS[2];
	struct S_VDATA       *VDATA;

	VDATA = malloc(sizeof *VDATA); // free
	VDATA->OPS = (struct S_OPERATORS_V const *const *) OPS;

	for (size_t i = 0; i < 2; i++)
		OPS[i] = malloc(sizeof *OPS[i]); // free

	if (strstr(DB.Form,"Weak")) {
		for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			init_VDATA(VDATA,VOLUME);

			// Obtain W_vI and Q_vI
			unsigned int const NvnI = VDATA->OPS[0]->NvnI;
			if (DB.Collocated) {
				VDATA->W_vI = VOLUME->What;
				VDATA->Q_vI = VOLUME->Qhat;
			} else {
				VDATA->W_vI = malloc(NvnI*Nvar * sizeof *(VDATA->W_vI)); // free
				VDATA->Q_vI = malloc(d         * sizeof *(VDATA->Q_vI)); // free
				for (size_t dim = 0; dim < d; dim++)
					VDATA->Q_vI[dim] = malloc(NvnI*Nvar * sizeof *(VDATA->Q_vI[dim])); // free

				coef_to_values_vI(VDATA,'W');
				coef_to_values_vI(VDATA,'Q');
			}

			// Compute negated Flux and its Jacobian in reference space
			double * const F_vI    = malloc(NvnI*d*Neq      * sizeof *F_vI),    // free
			       * const dFdW_vI = malloc(NvnI*d*Nvar*Neq * sizeof *dFdW_vI), // free
			       **const dFdQ_vI = malloc(d               * sizeof *dFdQ_vI); // free
			for (size_t dim = 0; dim < d; dim++)
				dFdQ_vI[dim] = malloc(NvnI*d*Nvar*Neq * sizeof *dFdQ_vI[dim]); // free

// Consider combining flux_viscous and jacobian_flux_viscous (ToBeDeleted)
			flux_viscous(NvnI,1,VDATA->W_vI,(const double *const *const) VDATA->Q_vI,F_vI);
			jacobian_flux_viscous(NvnI,1,VDATA->W_vI,(const double *const *const) VDATA->Q_vI,dFdW_vI,dFdQ_vI);

			if (!DB.Collocated) {
				free(VDATA->W_vI);
				array_free2_d(d,VDATA->Q_vI);
			}

			// Convert to reference space
			double *const Fr_vI = malloc(NvnI*d*Neq * sizeof *Fr_vI); // free
			convert_between_rp(NvnI,Neq,VOLUME->C_vI,F_vI,Fr_vI,"FluxToRef");
			free(F_vI);

			double *const dFrdW_vI = malloc(NvnI*d*Nvar*Neq * sizeof *dFrdW_vI); // free
			convert_between_rp(NvnI,Nvar*Neq,VOLUME->C_vI,dFdW_vI,dFrdW_vI,"FluxToRef");
			free(dFdW_vI);

			double **const dFrdQ_vI = malloc(d * sizeof *dFrdQ_vI); // free
			for (size_t dim = 0; dim < d; dim++) {
				dFrdQ_vI[dim] = malloc(NvnI*d*Nvar*Neq * sizeof *dFrdW_vI); // free
				convert_between_rp(NvnI,Nvar*Neq,VOLUME->C_vI,dFdQ_vI[dim],dFrdQ_vI[dim],"FluxToRef");
			}
			array_free2_d(d,dFdQ_vI);

			// Compute RHS and LHS terms
			unsigned int NvnS = VDATA->OPS[0]->NvnS;

			// RHS
			memset(VOLUME->RHS,0.0,NvnS*Neq * sizeof *(VOLUME->RHS));
			finalize_VOLUME_Inviscid_Weak(Neq,Fr_vI,VOLUME->RHS,'E',VDATA);
			free(Fr_vI);

			// LHS
// Add contribution of dFrdQ_vI below.
EXIT_UNSUPPORTED;
			memset(VOLUME->LHS,0.0,NvnS*NvnS*Neq*Nvar * sizeof *(VOLUME->LHS));
			finalize_VOLUME_Inviscid_Weak(NvnS*Neq*Nvar,dFrdW_vI,VOLUME->LHS,'I',VDATA);
			free(dFrdW_vI);
			array_free2_d(d,dFrdQ_vI);
		}
	} else if (strstr(DB.Form,"Strong")) {
		EXIT_UNSUPPORTED;
	}

	free(VDATA);
	for (size_t i = 0; i < 2; i++)
		free(OPS[i]);
}
