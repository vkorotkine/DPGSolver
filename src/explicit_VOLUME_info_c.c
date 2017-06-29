// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "explicit_VOLUME_info_c.h"
#include "S_VOLUME.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <stdbool.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"

#include "explicit_info_c.h"
#include "solver_functions.h"
#include "solver_functions_c.h"
#include "fluxes_structs.h"
#include "fluxes_inviscid_c.h"
#include "fluxes_viscous_c.h"
#include "array_free.h"
#include "array_print.h"

/*
 *	Purpose:
 *		Identical to explicit_VOLUME_info using complex variables (for complex step verification).
 *
 *	Comments:
 *		VOLUME is passed to the functions such that only the necessary RHS terms are evaluated when applicable. When
 *		compute_all is passed as TRUE, the RHS of every VOLUME is computed.
 *
 *	Notation:
 *
 *	References:
 */

static void compute_Inviscid_VOLUME_RHS_EFE (struct S_VOLUME *const VOLUME_perturbed, bool const compute_all);
static void compute_Viscous_VOLUME_RHS_EFE  (struct S_VOLUME *const VOLUME_perturbed, bool const compute_all);

void explicit_VOLUME_info_c (struct S_VOLUME *const VOLUME_perturbed, bool const compute_all)
{
	// Initialize DB Parameters
	unsigned int EFE        = DB.EFE,
	             Vectorized = DB.Vectorized;

	if (EFE) {
		switch (Vectorized) {
		case 0:
			compute_Inviscid_VOLUME_RHS_EFE(VOLUME_perturbed,compute_all);
			compute_Viscous_VOLUME_RHS_EFE(VOLUME_perturbed,compute_all);
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
		}
	} else {
		;
	}
}

static void compute_Inviscid_VOLUME_RHS_EFE (struct S_VOLUME *const VOLUME_perturbed, bool const compute_all)
{
	// Initialize DB Parameters
	unsigned int d    = DB.d,
				 Nvar = d+2,
				 Neq  = d+2;

	// Standard datatypes
	struct S_OPERATORS_V *OPS[2];

	struct S_VDATA *VDATA = malloc(sizeof *VDATA); // free
	VDATA->OPS = (struct S_OPERATORS_V const *const *) OPS;

	struct S_FLUX *const FLUXDATA = malloc(sizeof *FLUXDATA); // free
	FLUXDATA->PDE_index = DB.PDE_index;
	FLUXDATA->d   = d;
	FLUXDATA->Nel = 1;

	for (size_t i = 0; i < 2; i++)
		OPS[i] = malloc(sizeof *OPS[i]); // free

	struct S_LOCAL_MESH_ELEMENTS local_ELEMENTs;
	if (!compute_all)
		local_ELEMENTs = compute_local_ELEMENT_list(VOLUME_perturbed,'V');

	if (strstr(DB.Form,"Weak")) {
		for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			if (!compute_all) {
				if (!is_VOLUME_in_local_list(VOLUME,&local_ELEMENTs)) {
					continue;
				} else if (VOLUME != VOLUME_perturbed) {
					if (VOLUME->RHS_c)
						free(VOLUME->RHS_c);
					VOLUME->RHS_c = calloc((VOLUME->NvnS)*Neq , sizeof *(VOLUME->RHS_c)); // keep
					continue;
				}
			}

			init_VDATA(VDATA,VOLUME);

			// Obtain W_vI
			unsigned int NvnI = VDATA->OPS[0]->NvnI;
			if (DB.Collocated) {
				VDATA->W_vI_c = VOLUME->What_c;
			} else {
				VDATA->W_vI_c = malloc(NvnI*Nvar * sizeof *(VDATA->W_vI_c)); // free
				coef_to_values_vI_c(VDATA,'W');
			}

			// Compute Flux in reference space
			double complex *const F_vI = malloc(NvnI*d*Neq * sizeof *F_vI); // free

			FLUXDATA->Nn  = NvnI;
			FLUXDATA->W_c = VDATA->W_vI_c;
			FLUXDATA->F_c = F_vI;

			flux_inviscid_c(FLUXDATA);

			if (!DB.Collocated)
				free(VDATA->W_vI_c);

			// Convert to reference space
			double complex *const Fr_vI = malloc(NvnI*d*Neq * sizeof *Fr_vI); // free
			convert_between_rp_c(NvnI,Neq,VOLUME->C_vI,F_vI,Fr_vI,"FluxToRef");
			free(F_vI);

			// Compute RHS term
			unsigned int NvnS = VDATA->OPS[0]->NvnS;

			if (VOLUME->RHS_c)
				free(VOLUME->RHS_c);
			VOLUME->RHS_c = calloc(NvnS*Neq , sizeof *(VOLUME->RHS_c)); // keep

			finalize_VOLUME_Inviscid_Weak_c(Neq,Fr_vI,VOLUME->RHS_c,'E',VDATA);
			free(Fr_vI);
		}
	} else if (strstr(DB.Form,"Strong")) {
		EXIT_UNSUPPORTED;
	}

	free(VDATA);
	free(FLUXDATA);
	for (size_t i = 0; i < 2; i++)
		free(OPS[i]);
}

static void compute_Viscous_VOLUME_RHS_EFE (struct S_VOLUME *const VOLUME_perturbed, bool const compute_all)
{
	if (!DB.Viscous)
		return;

	unsigned int const d    = DB.d,
	                   Nvar = d+2,
	                   Neq  = d+2;

	struct S_OPERATORS_V *OPS[2];

	struct S_VDATA *VDATA = malloc(sizeof *VDATA); // free
	VDATA->OPS = (struct S_OPERATORS_V const *const *) OPS;

	struct S_FLUX *const FLUXDATA = malloc(sizeof *FLUXDATA); // free
	FLUXDATA->d   = d;
	FLUXDATA->Nel = 1;

	for (size_t i = 0; i < 2; i++)
		OPS[i] = malloc(sizeof *OPS[i]); // free

	struct S_LOCAL_MESH_ELEMENTS local_ELEMENTs;
	if (!compute_all)
		local_ELEMENTs = compute_local_ELEMENT_list(VOLUME_perturbed,'V');

	if (strstr(DB.Form,"Weak")) {
		for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			if (!compute_all) {
				// Note: Can't skip neighbouring VOLUMEs to VOLUME_perturbed as Qhat depends on What of VNeigh.
				if (!is_VOLUME_in_local_list(VOLUME,&local_ELEMENTs))
					continue;
			}

			init_VDATA(VDATA,VOLUME);

			// Obtain W_vI and Q_vI
			unsigned int const NvnI = VDATA->OPS[0]->NvnI;
			if (DB.Collocated) {
				VDATA->W_vI_c = VOLUME->What_c;
				VDATA->Q_vI_c = VOLUME->Qhat_c;
			} else {
				VDATA->W_vI_c = malloc(NvnI*Nvar * sizeof *(VDATA->W_vI_c)); // free
				VDATA->Q_vI_c = malloc(d         * sizeof *(VDATA->Q_vI_c)); // free
				for (size_t dim = 0; dim < d; dim++)
					VDATA->Q_vI_c[dim] = malloc(NvnI*Nvar * sizeof *(VDATA->Q_vI_c[dim])); // free

				coef_to_values_vI_c(VDATA,'W');
				coef_to_values_vI_c(VDATA,'Q');
			}

			// Compute negated Flux in reference space
			double complex *const F_vI = malloc(NvnI*d*Neq * sizeof *F_vI);

			FLUXDATA->Nn  = NvnI;
			FLUXDATA->W_c = VDATA->W_vI_c;
			FLUXDATA->Q_c = (double complex const *const *const) VDATA->Q_vI_c;
			FLUXDATA->F_c = F_vI;
if (VOLUME->indexg == 0) {
printf("eVi_c\n");
size_t dim = 0;
array_print_cmplx(VOLUME->NvnS,Neq,VDATA->W_vI_c,'C');
array_print_cmplx(VOLUME->NvnS,Neq,VDATA->Q_vI_c[dim],'C');
//EXIT_UNSUPPORTED;
}

			flux_viscous_c(FLUXDATA);

			if (!DB.Collocated) {
				free(VDATA->W_vI_c);
				array_free2_cmplx(d,VDATA->Q_vI_c);
			}

			// Convert to reference space
			double complex *const Fr_vI = malloc(NvnI*Neq*d * sizeof *Fr_vI); // free
			convert_between_rp_c(NvnI,Neq,VOLUME->C_vI,F_vI,Fr_vI,"FluxToRef");
			free(F_vI);

			// Compute RHS term
			finalize_VOLUME_Viscous_Weak_c(Neq,Fr_vI,VOLUME->RHS_c,'E',VDATA);
			free(Fr_vI);
		}
	} else if (strstr(DB.Form,"Strong")) {
		EXIT_UNSUPPORTED;
	}

	free(VDATA);
	free(FLUXDATA);
	for (size_t i = 0; i < 2; i++)
		free(OPS[i]);
}
