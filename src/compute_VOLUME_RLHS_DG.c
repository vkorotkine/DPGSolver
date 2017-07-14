// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "compute_VOLUME_RLHS_DG.h"
#include "solver.h"

#include <stdlib.h>
#include <string.h>

#include "petscmat.h"

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_VOLUME.h"

#include "solver_functions.h"
#include "fluxes_structs.h"
#include "fluxes_viscous.h"
#include "jacobian_fluxes_viscous.h"
#include "support.h"

/*
 *	Purpose:
 *		Evaluate the VOLUME contributions to the RHS and (optionally) LHS terms.
 *
 *	Comments:
 *		EFE stands for (E)xact (F)lux (E)valuation, meaning that the flux is not represented as a polynomial and then
 *		interpolated to the cubature nodes. EFE is the analogue of the (C)hain(R)ule approach for the strong form of the
 *		scheme. See Zwanenburg(2016) for additional discussion.
 *
 *		The inviscid and viscous contributions could be combined as they both require identical operations after the
 *		initial flux evaluation. Consider implementing this in the future based on profiling results. (ToBeModified)
 *
 *		For the Navier-Stokes solver, the VOLUME term flux has a dependence on the solution (W) and gradient (Q). As,
 *		the gradient has a dependence on the solution in adjacent VOLUMEs, this VOLUME term must receive contribution
 *		from both the VOLUME and the FACEs when performing the linearization. The 'compute_Viscous_VOLUME_EFE' function
 *		provides the VOLUME contribution to this part of the linearization and the associated function in
 *		compute_FACE_RLHS_DG adds the remaining part.
 *
 *		Computing fluxes and jacobians require many of the same operations and could be combined for speed-up in the
 *		future. A similar comment applies to the inviscid and viscous contributions. Profile and revisit (ToBeModified).
 *
 *	Notation:
 *
 *	References:
 *		Zwanenburg(2016)-Equivalence_between_the_Energy_Stable_Flux_Reconstruction_and_Discontinuous_Galerkin_Schemes
 */

struct S_LHS_info {
	size_t IndA[2], Nn[2];
	const double* LHS;

	InsertMode addv;
};
static struct S_LHS_info constructor_LHS_info (const double*const LHS, const struct S_VOLUME*const V0,
                                               const struct S_VOLUME*const V1, const InsertMode addv);
static void fill_PetscMat (const struct S_solver_info*const solver_info, const struct S_LHS_info*const LHS_info);

static void set_memory_to_zero_VOLUMEs  (const char imex_type);
static void compute_Inviscid_VOLUME_EFE (const struct S_solver_info*const solver_info);
static void compute_Viscous_VOLUME_EFE  (const char imex_type);

void compute_VOLUME_RLHS_DG (const struct S_solver_info*const solver_info)
{
	if (solver_info->display)
		printf("V");

	if (DB.EFE) {
		switch (DB.Vectorized) {
		case 0:
			set_memory_to_zero_VOLUMEs(solver_info->imex_type);
			compute_Inviscid_VOLUME_EFE(solver_info);
			compute_Viscous_VOLUME_EFE(solver_info->imex_type);
// Add source contribution here as well (ToBeDeleted)
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
		}
	} else {
		;
	}
}

static void set_memory_to_zero_VOLUMEs (const char imex_type)
{
	unsigned int const d    = DB.d,
	                   Nvar = DB.Nvar,
	                   Neq  = DB.Neq;

	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		unsigned int NvnS = VOLUME->NvnS;
		set_to_zero_d(NvnS*Neq,VOLUME->RHS);
		if (imex_type == 'I') {
			set_to_zero_d(NvnS*NvnS*Neq*Nvar,VOLUME->LHS);
			if (DB.Viscous) {
				for (size_t dim = 0; dim < d; dim++)
					set_to_zero_d(NvnS*NvnS*Neq*Nvar,VOLUME->LHSQ[dim]);
			}
		}
	}
}

static void compute_Inviscid_VOLUME_EFE (const struct S_solver_info*const solver_info)
{
	const char imex_type = solver_info->imex_type;

	if (!DB.Inviscid)
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
	DATA->imex_type = imex_type;

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

			compute_flux_inviscid(VDATA,FLUXDATA,imex_type);

			if (!DB.Collocated)
				manage_solver_memory(DATA,'F','W');

			if (DB.PDE_index == PDE_ADVECTION)
				manage_solver_memory(DATA,'F','X');

			// Convert to reference space
			convert_between_rp(NvnI,Neq,VOLUME->C_vI,FLUXDATA->F,FLUXDATA->Fr,"FluxToRef");
			if (imex_type == 'I')
				convert_between_rp(NvnI,Nvar*Neq,VOLUME->C_vI,FLUXDATA->dFdW,FLUXDATA->dFrdW,"FluxToRef");

			// Compute RHS terms
			finalize_VOLUME_Inviscid_Weak(Neq,FLUXDATA->Fr,VOLUME->RHS,'E',VDATA);

			// Compute LHS terms
			if (imex_type == 'I') {
				unsigned int NvnS = VDATA->OPS[0]->NvnS;
// Initialize memory for LHS, add to PetscMat, then free (ToBeModified)
// Initialize memory for LHS in manage_solver_memory for imex_type == 'I'
				double*const LHS = calloc(NvnS*NvnS*Neq*Nvar , sizeof *LHS); // free

// Remove addition of LHS in finalize_LHS to test equivalence of this treatment.
//				finalize_VOLUME_Inviscid_Weak(NvnS*Neq*Nvar,FLUXDATA->dFrdW,VOLUME->LHS,'I',VDATA);
				finalize_VOLUME_Inviscid_Weak(NvnS*Neq*Nvar,FLUXDATA->dFrdW,LHS,'I',VDATA);

				struct S_LHS_info LHS_info = constructor_LHS_info(LHS,VOLUME,VOLUME,ADD_VALUES);
				fill_PetscMat(solver_info,&LHS_info);
				free(LHS);
			}
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

static void compute_Viscous_VOLUME_EFE (const char imex_type)
{
	/*
	 *	Comments:
	 *		The viscous VOLUME contributions have a nearly identical form to those of the inviscid contributions.
	 *		Consider combining the two functions in the future. (ToBeModified)
	 */

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
	DATA->imex_type = imex_type;

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

			if (imex_type == 'E') {
				flux_viscous(FLUXDATA);
			} else if (imex_type == 'I') {
				if (!DB.Fv_func_of_W)
					FLUXDATA->dFdW = NULL;
				jacobian_flux_viscous(FLUXDATA);
			}

			if (!DB.Collocated) {
				manage_solver_memory(DATA,'F','W');
				manage_solver_memory(DATA,'F','Q');
			}

			// Convert to reference space
			convert_between_rp(NvnI,Neq,VOLUME->C_vI,FLUXDATA->F,FLUXDATA->Fr,"FluxToRef");
			if (imex_type == 'I') {
				if (DB.Fv_func_of_W)
					convert_between_rp(NvnI,Nvar*Neq,VOLUME->C_vI,FLUXDATA->dFdW,FLUXDATA->dFrdW,"FluxToRef");
				for (size_t dim = 0; dim < d; dim++)
					convert_between_rp(NvnI,Nvar*Neq,VOLUME->C_vI,FLUXDATA->dFdQ[dim],FLUXDATA->dFrdQ[dim],"FluxToRef");
			}


			// Compute RHS terms
			finalize_VOLUME_Viscous_Weak(Neq,FLUXDATA->Fr,VOLUME->RHS,'E',VDATA);

			// Compute LHS terms
			if (imex_type == 'I') {
				unsigned int NvnS = VDATA->OPS[0]->NvnS;
				if (DB.Fv_func_of_W)
					finalize_VOLUME_Viscous_Weak(NvnS*Neq*Nvar,FLUXDATA->dFrdW,VOLUME->LHS,'I',VDATA);

				initialize_VOLUME_LHSQ_Weak(NvnS*Neq*Nvar,(double const *const *const) FLUXDATA->dFrdQ,VOLUME->LHSQ,VDATA);
				finalize_VOLUME_LHSQV_Weak(VOLUME);
			}

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

static struct S_LHS_info constructor_LHS_info (const double*const LHS, const struct S_VOLUME*const V0,
                                               const struct S_VOLUME*const V1, const InsertMode addv)
{
	struct S_LHS_info LHS_info;

	LHS_info.LHS     = LHS;
	LHS_info.IndA[0] = V0->IndA;
	LHS_info.IndA[1] = V1->IndA;
	LHS_info.Nn[0]   = V0->NvnS;
	LHS_info.Nn[1]   = V1->NvnS;
	LHS_info.addv    = addv;

	return LHS_info;
}

static void fill_PetscMat (const struct S_solver_info*const solver_info, const struct S_LHS_info*const LHS_info)
{
	/*
	 *	Purpose:
	 *		Add individual LHS contributions to the global system matrix.
	 */

// need to call compute_dof before calling this function to set VOLUME->IndA (ToBeDeleted)
	Mat *A = solver_info->A;

	const double*const LHS = LHS_info->LHS;

	const size_t*const IndA = LHS_info->IndA,
	            *const Nn   = LHS_info->Nn;

	PetscInt idxm[Nn[0]],
	         idxn[Nn[1]];
	const InsertMode addv = LHS_info->addv;

	const unsigned int Neq  = DB.Neq,
	                   Nvar = DB.Nvar;

	for (size_t eq = 0; eq < Neq; eq++) {
		const size_t Indm = IndA[0] + eq*Nn[0];
		for (size_t i = 0; i < Nn[0]; i++)
			idxm[i] = Indm+i;

		for (size_t var = 0; var < Nvar; var++) {
			const size_t Indn = IndA[1] + var*Nn[1];
			for (size_t i = 0; i < Nn[1]; i++)
				idxn[i] = Indn+i;

			const PetscScalar*const vv = &LHS[(eq*Nvar+var)*Nn[0]*Nn[1]];

			MatSetValues(*A,Nn[0],idxm,Nn[1],idxn,vv,addv);
		}
	}
}

