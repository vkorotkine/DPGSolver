// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "implicit_VOLUME_info_DG.h"

#include <stdlib.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_VOLUME.h"

#include "solver_functions.h"
#include "fluxes_structs.h"
#include "jacobian_fluxes_viscous.h"
#include "support.h"

/*
 *	Purpose:
 *		Evaluate the VOLUME contributions to the RHS and LHS terms.
 *
 *	Comments:
 *		For the Navier-Stokes solver, the VOLUME term flux has a dependence on the solution (W) and gradient (Q). As,
 *		the gradient has a dependence on the solution in adjacent VOLUMEs, this VOLUME term must receive contribution
 *		from both the VOLUME and the FACEs when performing the linearization. The 'compute_Viscous_VOLUME_VOLUME_EFE'
 *		function provides the VOLUME contribution to this part of the linearization and the associated function in
 *		implicit_FACE_info adds the remaining part.
 *
 *		Computing fluxes and jacobians require many of the same operations and could be combined for speed-up in the
 *		future. A similar comment applies to the inviscid and viscous contributions. Profile and revisit (ToBeModified).
 *
 *	Notation:
 *
 *	References:
 */

static void compute_Inviscid_VOLUME_EFE(void);
static void compute_Viscous_VOLUME_EFE(void);

void implicit_VOLUME_info_DG (const bool PrintEnabled)
{
	if (PrintEnabled)
		printf("V");

	if (DB.EFE) {
		switch (DB.Vectorized) {
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
				 Nvar = DB.Nvar,
				 Neq  = DB.Neq;

	if (!DB.Inviscid) {
		for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			unsigned int NvnS = VOLUME->NvnS;
			set_to_zero_d(NvnS*Neq,VOLUME->RHS);
			set_to_zero_d(NvnS*NvnS*Neq*Nvar,VOLUME->LHS);
		}
		return;
	}
	struct S_OPERATORS_V *OPS[2];

	struct S_VDATA *const VDATA = malloc(sizeof *VDATA); // free
	VDATA->OPS = (struct S_OPERATORS_V const *const *) OPS;

	struct S_FLUX *const FLUXDATA = malloc(sizeof *FLUXDATA); // free
	FLUXDATA->d   = d;
	FLUXDATA->Nel = 1;
	FLUXDATA->PDE_index = DB.PDE_index;

	struct S_DATA *const DATA = malloc(sizeof *DATA); // free
	DATA->VDATA     = VDATA;
	DATA->FLUXDATA  = FLUXDATA;
	DATA->feature   = 'V';
	DATA->imex_type = 'I';

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
				manage_solver_memory(DATA,'A','W'); // free
				coef_to_values_vI(VDATA,'W');
			}

			// Compute Flux and its Jacobian in reference space
			manage_solver_memory(DATA,'A','I'); // free

			if (DB.PDE_index == PDE_ADVECTION)
				manage_solver_memory(DATA,'A','X'); // free

			compute_flux_inviscid(VDATA,FLUXDATA,'I');

			if (!DB.Collocated)
				manage_solver_memory(DATA,'F','W');

			if (DB.PDE_index == PDE_ADVECTION)
				manage_solver_memory(DATA,'F','X');

			// Convert to reference space
			convert_between_rp(NvnI,Neq,VOLUME->C_vI,FLUXDATA->F,FLUXDATA->Fr,"FluxToRef");
			convert_between_rp(NvnI,Nvar*Neq,VOLUME->C_vI,FLUXDATA->dFdW,FLUXDATA->dFrdW,"FluxToRef");

			// Compute RHS and LHS terms
			unsigned int NvnS = VDATA->OPS[0]->NvnS;

			// RHS
			set_to_zero_d(NvnS*Neq,VOLUME->RHS);
			finalize_VOLUME_Inviscid_Weak(Neq,FLUXDATA->Fr,VOLUME->RHS,'E',VDATA);

			// LHS
			set_to_zero_d(NvnS*NvnS*Neq*Nvar,VOLUME->LHS);
			finalize_VOLUME_Inviscid_Weak(NvnS*Neq*Nvar,FLUXDATA->dFrdW,VOLUME->LHS,'I',VDATA);

			manage_solver_memory(DATA,'F','I');
		}
	} else if (strstr(DB.Form,"Strong")) {
		EXIT_UNSUPPORTED;
	}

	free(VDATA);
	free(FLUXDATA);
	for (size_t i = 0; i < 2; i++)
		free(OPS[i]);
	free(DATA);
}

static void compute_Viscous_VOLUME_EFE(void)
{
	if (!DB.Viscous)
		return;

	unsigned int d    = DB.d,
				 Nvar = DB.Nvar,
				 Neq  = DB.Neq;

	struct S_OPERATORS_V *OPS[2];

	struct S_VDATA *const VDATA = malloc(sizeof *VDATA); // free
	VDATA->OPS = (struct S_OPERATORS_V const *const *) OPS;

	struct S_FLUX *const FLUXDATA = malloc(sizeof *FLUXDATA); // free
	FLUXDATA->d   = d;
	FLUXDATA->Nel = 1;
	FLUXDATA->PDE_index = DB.PDE_index;

	struct S_DATA *const DATA = malloc(sizeof *DATA); // free
	DATA->VDATA     = VDATA;
	DATA->FLUXDATA  = FLUXDATA;
	DATA->feature   = 'V';
	DATA->imex_type = 'I';

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
				manage_solver_memory(DATA,'A','W'); // free
				manage_solver_memory(DATA,'A','Q'); // free

				coef_to_values_vI(VDATA,'W');
				coef_to_values_vI(VDATA,'Q');
			}

			// Compute negated Flux and its Jacobian in reference space
			manage_solver_memory(DATA,'A','V'); // free

			FLUXDATA->Nn = NvnI;
			FLUXDATA->W  = VDATA->W_vI;
			FLUXDATA->Q  = (double const *const *const) VDATA->Q_vI;
			if (!DB.Fv_func_of_W)
				FLUXDATA->dFdW = NULL;

			jacobian_flux_viscous(FLUXDATA);

			if (!DB.Collocated) {
				manage_solver_memory(DATA,'F','W');
				manage_solver_memory(DATA,'F','Q');
			}

			// Convert to reference space
			convert_between_rp(NvnI,Neq,VOLUME->C_vI,FLUXDATA->F,FLUXDATA->Fr,"FluxToRef");
			if (DB.Fv_func_of_W)
				convert_between_rp(NvnI,Nvar*Neq,VOLUME->C_vI,FLUXDATA->dFdW,FLUXDATA->dFrdW,"FluxToRef");
			for (size_t dim = 0; dim < d; dim++)
				convert_between_rp(NvnI,Nvar*Neq,VOLUME->C_vI,FLUXDATA->dFdQ[dim],FLUXDATA->dFrdQ[dim],"FluxToRef");

			// Compute RHS and LHS terms
			unsigned int NvnS = VDATA->OPS[0]->NvnS;

			// RHS
			finalize_VOLUME_Viscous_Weak(Neq,FLUXDATA->Fr,VOLUME->RHS,'E',VDATA);

			// LHS
			if (DB.Fv_func_of_W)
				finalize_VOLUME_Viscous_Weak(NvnS*Neq*Nvar,FLUXDATA->dFrdW,VOLUME->LHS,'I',VDATA);

			for (size_t dim = 0; dim < d; dim++)
				set_to_zero_d(NvnS*NvnS*Neq*Nvar,VOLUME->LHSQ[dim]);
			initialize_VOLUME_LHSQ_Weak(NvnS*Neq*Nvar,(double const *const *const) FLUXDATA->dFrdQ,VOLUME->LHSQ,VDATA);
			finalize_VOLUME_LHSQV_Weak(VOLUME);

			manage_solver_memory(DATA,'F','V');
		}
	} else if (strstr(DB.Form,"Strong")) {
		EXIT_UNSUPPORTED;
	}

	free(VDATA);
	free(FLUXDATA);
	for (size_t i = 0; i < 2; i++)
		free(OPS[i]);
	free(DATA);
}
