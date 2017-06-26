// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "solver_functions_c.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_VOLUME.h"
#include "S_FACE.h"
#include "Test.h"

#include "solver_functions.h"
#include "matrix_functions.h"

#include "array_swap.h"
#include "boundary_conditions.h"
#include "boundary_conditions_c.h"
#include "fluxes_structs.h"
#include "fluxes_inviscid_c.h"
#include "fluxes_viscous_c.h"

#include "array_free.h"
#include "array_print.h"

#include "support.h"

/*
 *	Purpose:
 *		Provide solver related functions for complex variables (used for complex step verification).
 *
 *	Comments:
 *		Only supports imex_type = 'E' as the complex solver functions are not used to compute linearizations.
 *
 *	Notation:
 *
 *	References:
 */

// **************************************************************************************************** //
// General functions
// **************************************************************************************************** //

void manage_solver_memory_c(struct S_DATA *const DATA, char const mem_op, char const mem_type)
{
	char const feature   = DATA->feature,
	           imex_type = DATA->imex_type;

	if (!(imex_type == 'E'))
		EXIT_UNSUPPORTED;

	struct S_VDATA         *const VDATA     = DATA->VDATA;
	struct S_FDATA         *const FDATAL    = DATA->FDATAL,
	                       *const FDATAR    = DATA->FDATAR;
	struct S_FLUX          *const FLUXDATA  = DATA->FLUXDATA;
	struct S_NUMERICALFLUX *const NFLUXDATA = DATA->NFLUXDATA;


	struct S_OPERATORS_F const *const *const OPSF = FDATAL->OPS;

	unsigned int const d        = DB.d,
	                   Nvar     = DB.Nvar,
	                   Neq      = DB.Neq,
	                   NvnI     = VDATA->OPS[0]->NvnI,
	                   IndFType = FDATAL->IndFType,
	                   NfnI     = OPSF[IndFType]->NfnI;

	if (mem_op == 'A') {
		if (feature == 'V') {
			if (mem_type == 'W') {
				VDATA->W_vI_c = malloc(NvnI*Nvar * sizeof *(VDATA->W_vI_c)); // keep
			} else if (mem_type == 'Q') {
				VDATA->Q_vI_c = malloc(d * sizeof *(VDATA->Q_vI_c)); // keep
				for (size_t dim = 0; dim < d; dim++)
					VDATA->Q_vI_c[dim] = malloc(NvnI*Nvar * sizeof *(VDATA->Q_vI_c[dim])); // keep
			} else if (mem_type == 'I' || mem_type == 'V') {
				FLUXDATA->F_c  = malloc(NvnI*d*Neq * sizeof *(FLUXDATA->F_c));  // keep
				FLUXDATA->Fr_c = malloc(NvnI*Neq*d * sizeof *(FLUXDATA->Fr_c)); // keep
			} else {
				EXIT_UNSUPPORTED;
			}
		} else if (feature == 'F') {
			if (mem_type == 'W') {
				FDATAL->W_fIL_c = malloc(NfnI*Nvar * sizeof *(FDATAL->W_fIL_c)), // keep
				FDATAR->W_fIL_c = malloc(NfnI*Nvar * sizeof *(FDATAR->W_fIL_c)); // keep
			} else if (mem_type == 'Q') {
				FDATAL->Qp_fIL_c = malloc(d * sizeof *(FDATAL->Qp_fIL_c)), // keep
				FDATAR->Qp_fIL_c = malloc(d * sizeof *(FDATAR->Qp_fIL_c)); // keep
				for (size_t dim = 0; dim < d; dim++) {
					FDATAL->Qp_fIL_c[dim] = malloc(NfnI*Nvar * sizeof *(FDATAL->Qp_fIL_c[dim])); // keep
					FDATAR->Qp_fIL_c[dim] = malloc(NfnI*Nvar * sizeof *(FDATAR->Qp_fIL_c[dim])); // keep
				}
			} else if (mem_type == 'S') {
				NFLUXDATA->nSolNum_c = malloc(d * sizeof *(NFLUXDATA->nSolNum_c)); // keep
				for (size_t dim = 0; dim < d; dim++)
					NFLUXDATA->nSolNum_c[dim] = malloc(NfnI*Neq * sizeof *(NFLUXDATA->nSolNum_c[dim])); // keep
			} else if (mem_type == 'I') {
				NFLUXDATA->nFluxNum_c = malloc(NfnI*Neq * sizeof *(NFLUXDATA->nFluxNum_c)); // keep
			} else if (mem_type == 'V') {
				NFLUXDATA->nFluxNum_c = malloc(NfnI*Neq * sizeof *(NFLUXDATA->nFluxNum_c)); // keep
			} else if (mem_type == 'P') {
				NFLUXDATA->nFluxNum_c = malloc(NfnI*Neq * sizeof *(NFLUXDATA->nFluxNum_c)); // keep
			} else {
				EXIT_UNSUPPORTED;
			}
		} else {
			EXIT_UNSUPPORTED;
		}
	} else if (mem_op == 'F') {
		if (feature == 'V') {
			if (mem_type == 'W') {
				free(VDATA->W_vI_c);
			} else if (mem_type == 'Q') {
				array_free2_cmplx(d,VDATA->Q_vI_c);
			} else if (mem_type == 'I' || mem_type == 'V') {
				free(FLUXDATA->F_c);
				free(FLUXDATA->Fr_c);
			} else {
				EXIT_UNSUPPORTED;
			}
		} else if (feature == 'F') {
			if (mem_type == 'W') {
				free(FDATAL->W_fIL_c);
				free(FDATAR->W_fIL_c);
			} else if (mem_type == 'Q') {
				array_free2_cmplx(d,FDATAL->Qp_fIL_c);
				array_free2_cmplx(d,FDATAR->Qp_fIL_c);
			} else if (mem_type == 'S') {
				array_free2_cmplx(d,NFLUXDATA->nSolNum_c);
			} else if (mem_type == 'I') {
				free(NFLUXDATA->nFluxNum_c);
			} else if (mem_type == 'V') {
				free(NFLUXDATA->nFluxNum_c);
			} else if (mem_type == 'P') {
				free(NFLUXDATA->nFluxNum_c);
			}
		} else {
			EXIT_UNSUPPORTED;
		}
	} else {
		EXIT_UNSUPPORTED;
	}
}

// **************************************************************************************************** //
// VOLUME functions
// **************************************************************************************************** //

void coef_to_values_vI_c(struct S_VDATA const *const VDATA, char const coef_type)
{
	if (DB.Collocated) {
		EXIT_UNSUPPORTED; // Set A_vI = VOLUME->Ahat in calling function.
	} else {
		if (!(coef_type == 'W' || coef_type == 'Q'))
			EXIT_UNSUPPORTED;

		unsigned int const d    = DB.d,
		                   Nvar = DB.Nvar;

		struct S_OPERATORS_V const *const *const OPS    = (struct S_OPERATORS_V const *const *const) VDATA->OPS;
		struct S_VOLUME      const *const        VOLUME = VDATA->VOLUME;

		if (coef_type == 'W') {
			mm_dcc(CBCM,CBT,CBNT,OPS[0]->NvnI,Nvar,OPS[0]->NvnS,1.0,0.0,OPS[0]->ChiS_vI,VOLUME->What_c,VDATA->W_vI_c);
		} else if (coef_type == 'Q') {
			for (size_t dim = 0; dim < d; dim++)
				mm_dcc(CBCM,CBT,CBNT,OPS[0]->NvnI,Nvar,OPS[0]->NvnS,1.0,0.0,OPS[0]->ChiS_vI,VOLUME->Qhat_c[dim],VDATA->Q_vI_c[dim]);
		}
	}
}

void convert_between_rp_c(unsigned int const Nn, unsigned int const Nrc, double const *const C,
                          double complex *const Ap, double complex *const Ar, char const *const conv_type)
{
	unsigned int d = DB.d;

	if (strstr(conv_type,"FluxToRef")) {
		set_to_zero_cmplx(Nn*Nrc*d,Ar);
		for (size_t col = 0; col < Nrc; col++) {
		for (size_t dim1 = 0; dim1 < d; dim1++) {
			size_t const IndAr = (Nrc*dim1+col)*Nn;
			for (size_t dim2 = 0; dim2 < d; dim2++) {
				size_t const IndAp = (col*d+dim2)*Nn,
				             IndC  = (dim1*d+dim2)*Nn;
				for (size_t n = 0; n < Nn; n++)
					Ar[IndAr+n] += Ap[IndAp+n]*C[IndC+n];
			}
		}}
	} else {
		EXIT_UNSUPPORTED;
	}
}

void finalize_VOLUME_Inviscid_Weak_c(unsigned int const Nrc, double complex const *const Ar_vI,
                                     double complex *const RLHS, char const imex_type,
                                     struct S_VDATA const *const VDATA)
{
	struct S_OPERATORS_V const *const *const OPS = (struct S_OPERATORS_V const *const *const) VDATA->OPS;

	unsigned int const d    = DB.d,
	                   NvnS = OPS[0]->NvnS,
	                   NvnI = OPS[0]->NvnI;

	if (imex_type == 'E') {
		for (size_t dim = 0; dim < d; dim++)
			mm_dcc(CBCM,CBT,CBNT,NvnS,Nrc,NvnI,1.0,1.0,OPS[0]->D_Weak[dim],&Ar_vI[NvnI*Nrc*dim],RLHS);
	} else {
		EXIT_UNSUPPORTED;
	}
}

void finalize_VOLUME_Viscous_Weak_c(unsigned int const Nrc, double complex *const Ar_vI, double complex *const RLHS,
                                    char const imex_type, struct S_VDATA const *const VDATA)
{
	finalize_VOLUME_Inviscid_Weak_c(Nrc,Ar_vI,RLHS,imex_type,VDATA);
}

// **************************************************************************************************** //
// FACE functions
// **************************************************************************************************** //

void coef_to_values_fI_c(struct S_FDATA *const FDATA, char const coef_type, char const imex_type)
{
	if (!(imex_type == 'E'))
		EXIT_UNSUPPORTED;

	if (!(coef_type == 'W' || coef_type == 'Q'))
		EXIT_UNSUPPORTED;

	struct S_OPERATORS_F const *const *const OPS    = (struct S_OPERATORS_F const *const *const) FDATA->OPS;
	struct S_VOLUME      const *const        VOLUME = FDATA->VOLUME;
	struct S_FACE              *const        FACE   = (struct S_FACE *const) FDATA->FACE;

	unsigned int const d        = DB.d,
	                   Nvar     = DB.Nvar,
	                   Vf       = FDATA->Vf,
	                   IndFType = FDATA->IndFType,
	                   NfnI     = OPS[IndFType]->NfnI,
	                   NvnS     = OPS[0]->NvnS;

	double complex **Qphat_c = NULL;
	if (coef_type == 'Q') {
		Qphat_c = malloc(d * sizeof *Qphat_c);

		// Add VOLUME contribution
		for (size_t dim = 0; dim < d; dim++) {
			Qphat_c[dim] = malloc(NvnS*Nvar * sizeof *Qphat_c[dim]); // free
			for (size_t i = 0; i < NvnS*Nvar; i++)
				Qphat_c[dim][i] = VOLUME->QhatV_c[dim][i];
		}

		if (DB.ViscousFluxType == FLUX_CDG2) {
			// Evaluate from which side scaling should be computed based on area switch (Brdar(2012), eq. (4.5))
			if (FDATA->side == 'L') {
				if (FACE->Boundary) {
					FACE->CDG2_side = 'L';
				} else {
					struct S_VOLUME const *const VL = FACE->VL,
					                      *const VR = FACE->VR;

					unsigned int const NvnSL = VL->NvnS,
					                   NvnSR = VR->NvnS;

					double const *const detJVL_vIL = VL->detJV_vI,
					             *const detJVR_vIR = VR->detJV_vI;

// Potentially scale by volume of reference element (necessary for mixed meshes?) (ToBeDeleted)
					double VolL = 0.0;
					for (size_t n = 0; n < NvnSL; n++)
						VolL = max(VolL,detJVL_vIL[n]);

					double VolR = 0.0;
					for (size_t n = 0; n < NvnSR; n++)
						VolR = max(VolR,detJVR_vIR[n]);

					if (VolL <= VolR)
						FACE->CDG2_side = 'L';
					else
						FACE->CDG2_side = 'R';
				}
			}
		}

		unsigned int const Boundary = FACE->Boundary;

		double chi = 0.0;
		if (DB.ViscousFluxType == FLUX_CDG2) {
			TestDB.EnteredViscousFlux[0]++;
			if (FDATA->side == FACE->CDG2_side) {
				// Multiplied by 2.0 as the viscous flux terms from each side are subsequently averaged.
				if (Boundary)
					chi =     (DB.NfMax);
				else
					chi = 2.0*(DB.NfMax);
			}
		} else if (DB.ViscousFluxType == FLUX_BR2) {
			TestDB.EnteredViscousFlux[1]++;
			if (Boundary)
				chi =     (DB.NfMax);
			else
				chi =     (DB.NfMax);
		} else {
			EXIT_UNSUPPORTED;
		}
		chi *= PENALIZATION_SCALING;

		// Add partial FACE contribution
		if (chi != 0.0) {
			for (size_t dim = 0; dim < d; dim++) {
				for (size_t i = 0; i < NvnS*Nvar; i++)
					Qphat_c[dim][i] += chi*FDATA->QhatF_c[dim][i];
			}
		}
	}

	if (coef_type == 'W') {
		mm_dcc(CBCM,CBT,CBNT,NfnI,Nvar,OPS[0]->NvnS,1.0,0.0,OPS[0]->ChiS_fI[Vf],VOLUME->What_c,FDATA->W_fIL_c);
	} else if (coef_type == 'Q') {
		for (size_t dim = 0; dim < d; dim++)
			mm_dcc(CBCM,CBT,CBNT,NfnI,Nvar,OPS[0]->NvnS,1.0,0.0,OPS[0]->ChiS_fI[Vf],Qphat_c[dim],FDATA->Qp_fIL_c[dim]);
	}

	if (coef_type == 'Q')
		array_free2_cmplx(d,Qphat_c);
}

static void evaluate_boundary_c(struct S_BC *const BCdata, bool const ComputeGradient)
{
	unsigned int const BC   = BCdata->BC;

	if (ComputeGradient) {
		if (!(BC % BC_STEP_SC == BC_NOSLIP_T         ||
		      BC % BC_STEP_SC == BC_NOSLIP_ADIABATIC ||
		      BC % BC_STEP_SC == BC_DIRICHLET        ||
		      BC % BC_STEP_SC == BC_NEUMANN)) {
			EXIT_UNSUPPORTED;
		}
	}

	compute_boundary_values_c(BCdata);
}

void compute_WR_fIL_c(struct S_FDATA const *const FDATA, double complex const *const WL_fIL,
                      double complex *const WR_fIL)
{
	struct S_OPERATORS_F const *const *const OPS  = FDATA->OPS;
	struct S_FACE        const *const        FACE = FDATA->FACE;

	unsigned int const Nvar     = DB.Nvar,
	                   IndFType = FDATA->IndFType,
	                   BC       = FACE->BC,
	                   NfnI     = OPS[IndFType]->NfnI,
	                   *const nOrdRL = OPS[IndFType]->nOrdRL;

	double const *const XYZ_fIL = FACE->XYZ_fI,
	             *const n_fIL   = FACE->n_fI;

	if (BC == 0 || (BC % BC_STEP_SC > 50)) { // Internal/Periodic FACE
		coef_to_values_fI_c((struct S_FDATA *const) FDATA,'W','E');
		array_rearrange_cmplx(NfnI,Nvar,nOrdRL,'C',WR_fIL);
	} else { // Boundary FACE
		struct S_BC *const BCdata = malloc(sizeof *BCdata); // free

		BCdata->d   = DB.d;
		BCdata->Nn  = NfnI;
		BCdata->Nel = 1;
		BCdata->BC  = BC;

		BCdata->XYZ  = XYZ_fIL;
		BCdata->nL   = n_fIL;
		BCdata->WL_c = WL_fIL;
		BCdata->WB_c = WR_fIL;
		BCdata->QL_c = NULL;
		BCdata->QB_c = NULL;

		evaluate_boundary_c(BCdata,0);
		free(BCdata);
	}
}

void compute_WR_QpR_fIL_c(struct S_FDATA const *const FDATA, double complex const *const WL_fIL,
                             double complex *const WR_fIL, double complex const *const *const QpL_fIL,
                             double complex *const *const QpR_fIL, char const imex_type)
{
	struct S_OPERATORS_F const *const *const OPS  = FDATA->OPS;
	struct S_FACE        const *const        FACE = FDATA->FACE;

	unsigned int const d        = DB.d,
	                   Nvar     = DB.Nvar,
	                   IndFType = FDATA->IndFType,
	                   BC       = FACE->BC,
	                   NfnI     = OPS[IndFType]->NfnI,
	                   *const nOrdRL = OPS[IndFType]->nOrdRL;

	if (!FACE->Boundary) { // Internal/Periodic FACE
		compute_WR_fIL_c(FDATA,WL_fIL,WR_fIL);

		coef_to_values_fI_c((struct S_FDATA *const) FDATA,'Q',imex_type);
		for (size_t dim = 0; dim < d; dim++)
			array_rearrange_cmplx(NfnI,Nvar,nOrdRL,'C',QpR_fIL[dim]);
	} else { // Boundary FACE
		struct S_BC *const BCdata = malloc(sizeof *BCdata); // free

		BCdata->d   = DB.d;
		BCdata->Nn  = NfnI;
		BCdata->Nel = 1;
		BCdata->BC  = BC;

		BCdata->XYZ  = FACE->XYZ_fI;
		BCdata->nL   = FACE->n_fI;
		BCdata->WL_c = WL_fIL;
		BCdata->WB_c = WR_fIL;
		BCdata->QL_c = QpL_fIL;
		BCdata->QB_c = QpR_fIL;

		evaluate_boundary_c(BCdata,1);
		free(BCdata);
	}
}

void compute_numerical_flux_c(struct S_FDATA const *const FDATA, char const imex_type)
{
	if (imex_type != 'E')
		EXIT_UNSUPPORTED;

	struct S_OPERATORS_F const *const *const OPS  = (struct S_OPERATORS_F const *const *const) FDATA->OPS;
	struct S_FACE        const *const        FACE = FDATA->FACE;

	unsigned int const d        = DB.d,
	                   IndFType = FDATA->IndFType,
	                   NfnI     = OPS[IndFType]->NfnI;

	double const *const n_fIL = FACE->n_fI;

	double complex const *const WL_fIL       = FDATA->NFLUXDATA->WL_c,
	                     *const WR_fIL       = FDATA->NFLUXDATA->WR_c;
	double complex       *const nFluxNum_fIL = FDATA->NFLUXDATA->nFluxNum_c;

	struct S_NUMERICALFLUX *const NUMFLUXDATA = malloc(sizeof *NUMFLUXDATA); // free
	NUMFLUXDATA->NumFluxInviscid_index = DB.InviscidFluxType;
	NUMFLUXDATA->d   = d;
	NUMFLUXDATA->Nn  = NfnI;
	NUMFLUXDATA->Nel = 1;

	NUMFLUXDATA->nL       = n_fIL;
	NUMFLUXDATA->XYZ      = FACE->XYZ_fI;
	NUMFLUXDATA->WL_c       = WL_fIL;
	NUMFLUXDATA->WR_c       = WR_fIL;
	NUMFLUXDATA->nFluxNum_c = nFluxNum_fIL;

	flux_num_inviscid_c(NUMFLUXDATA);
	free(NUMFLUXDATA);
}

void compute_numerical_solution_c(struct S_FDATA const *const FDATA, char const imex_type)
{
	if (!(imex_type == 'E'))
		EXIT_UNSUPPORTED;

	struct S_OPERATORS_F const *const *const OPS  = (struct S_OPERATORS_F const *const *const) FDATA->OPS;
	struct S_FACE        const *const        FACE = FDATA->FACE;

	unsigned int const d        = DB.d,
	                   Nvar     = DB.Nvar,
	                   IndFType = FDATA->IndFType,
	                   NfnI     = OPS[IndFType]->NfnI;

	double const *const n_fIL = FACE->n_fI;

	double complex const *const WL_fIL = FDATA->NFLUXDATA->WL_c,
	                     *const WR_fIL = FDATA->NFLUXDATA->WR_c;

	double complex       *const *const nSolNum_fIL = FDATA->NFLUXDATA->nSolNum_c;

	for (size_t dim = 0; dim < d; dim++) {
	for (size_t var = 0; var < Nvar; var++) {
	for (size_t n = 0; n < NfnI; n++) {
		nSolNum_fIL[dim][var*NfnI+n] = n_fIL[n*d+dim]*0.5*(WL_fIL[var*NfnI+n]+WR_fIL[var*NfnI+n]);
	}}}
}

static void correct_numerical_solution_strong_c(struct S_FDATA const *const FDATA, char const imex_type, char const side)
{
	if (!(imex_type == 'E'))
		EXIT_UNSUPPORTED;

	struct S_OPERATORS_F const *const *const OPS  = (struct S_OPERATORS_F const *const *const) FDATA->OPS;
	struct S_FACE        const *const        FACE = FDATA->FACE;

	unsigned int const d        = DB.d,
	                   Nvar     = DB.Nvar,
	                   IndFType = FDATA->IndFType,
	                   NfnI     = OPS[IndFType]->NfnI;

	double const *const n_fIL     = FACE->n_fI,
	             *const detJF_fIL = FACE->detJF_fI;

	double complex const *const WL_fIL = FDATA->NFLUXDATA->WL_c,
	                     *const WR_fIL = FDATA->NFLUXDATA->WR_c;

	double complex       *const *const nSolNum_fIL = FDATA->NFLUXDATA->nSolNum_c;

	if (side == 'L') {
		for (size_t dim = 0; dim < d; dim++) {
		for (size_t var = 0; var < Nvar; var++) {
		for (size_t n = 0; n < NfnI; n++) {
			nSolNum_fIL[dim][var*NfnI+n] -= n_fIL[n*d+dim]*detJF_fIL[n]*WL_fIL[var*NfnI+n];
		}}}
	} else if (side == 'R') {
		for (size_t dim = 0; dim < d; dim++) {
		for (size_t var = 0; var < Nvar; var++) {
		for (size_t n = 0; n < NfnI; n++) {
			nSolNum_fIL[dim][var*NfnI+n] += n_fIL[n*d+dim]*detJF_fIL[n]*(WL_fIL[var*NfnI+n]-WR_fIL[var*NfnI+n]);
		}}}
	} else {
		EXIT_UNSUPPORTED;
	}
}

static void dot_with_normal_c(unsigned int const Nn, unsigned int const NCol, double const *const n_fIL,
                            double complex const *const ANum_fIL, double complex *const nANum_fIL)
{
	unsigned int const d = DB.d;

	set_to_zero_cmplx(Nn*NCol,nANum_fIL);
	for (size_t col = 0; col < NCol; col++) {
		double complex const *const ANum_ptr = &ANum_fIL[Nn*d*col];
		for (size_t dim = 0; dim < d; dim++) {
		for (size_t n = 0; n < Nn; n++) {
			nANum_fIL[n+Nn*col] += n_fIL[n*d+dim]*ANum_ptr[Nn*dim+n];
		}}
	}
}

void compute_numerical_flux_viscous_c(struct S_FDATA const *const FDATAL, struct S_FDATA const *const FDATAR,
                                      char const imex_type)
{
	if (!(imex_type == 'E'))
		EXIT_UNSUPPORTED;

	struct S_OPERATORS_F const *const *const OPSL = (struct S_OPERATORS_F const *const *const) FDATAL->OPS;
	struct S_FACE        const *const        FACE = FDATAL->FACE;

	char const *PDE = DB.PDE;

	unsigned int const d        = DB.d,
	                   Nvar     = DB.Nvar,
	                   Neq      = DB.Neq,
	                   IndFType = FDATAL->IndFType,
	                   Boundary = FACE->Boundary,
	                   NfnI     = OPSL[IndFType]->NfnI;

	double const *const n_fIL                  = FACE->n_fI;
	double complex const *const WL_fIL         = FDATAL->NFLUXDATA->WL_c,
	                     *const WR_fIL         = FDATAL->NFLUXDATA->WR_c,
	                     *const *const QpL_fIL = (double complex const *const *const) FDATAL->Qp_fIL_c,
	                     *const *const QpR_fIL = (double complex const *const *const) FDATAR->Qp_fIL_c;

	double complex *const nFluxViscNum_fIL = FDATAL->NFLUXDATA->nFluxNum_c;

	double complex *const FluxViscNum_fIL = malloc(NfnI*d*Neq * sizeof *FluxViscNum_fIL); // free

	struct S_FLUX *const FLUXDATA = malloc(sizeof *FLUXDATA); // free
	FLUXDATA->d   = d;
	FLUXDATA->Nn  = NfnI;
	FLUXDATA->Nel = 1;

	if (Boundary) {
		double complex *const W_fIL = malloc(NfnI*Nvar * sizeof *W_fIL); // free
		for (size_t i = 0; i < NfnI*Nvar; i++)
			W_fIL[i] = 0.5*(WL_fIL[i]+WR_fIL[i]);

		double complex **Qp_fIL = malloc(d * sizeof *Qp_fIL); // free
		for (size_t dim = 0; dim < d; dim++) {
			Qp_fIL[dim] = malloc(NfnI*Nvar * sizeof *Qp_fIL[dim]); // free
			for (size_t i = 0; i < NfnI*Nvar; i++)
				Qp_fIL[dim][i] = 0.5*(QpL_fIL[dim][i]+QpR_fIL[dim][i]);
		}

		if (strstr(PDE,"NavierStokes")) {
			FLUXDATA->W_c = W_fIL;
			FLUXDATA->Q_c = (double complex const *const *const) Qp_fIL;
			FLUXDATA->F_c = FluxViscNum_fIL;
			flux_viscous_c(FLUXDATA);
		} else if (strstr(PDE,"Poisson")) {
			for (size_t dim = 0; dim < d; dim++) {
				for (size_t i = 0, iMax = NfnI*Nvar; i < iMax; i++)
					FluxViscNum_fIL[iMax*dim+i] = Qp_fIL[dim][i];
			}
		} else {
			EXIT_UNSUPPORTED;
		}

		free(W_fIL);
		array_free2_cmplx(d,Qp_fIL);
	} else {
		if (strstr(PDE,"NavierStokes")) {
			double complex *const FluxViscL_fIL = malloc(NfnI*d*Neq * sizeof *FluxViscL_fIL); // free
			double complex *const FluxViscR_fIL = malloc(NfnI*d*Neq * sizeof *FluxViscR_fIL); // free

			FLUXDATA->W_c = WL_fIL;
			FLUXDATA->Q_c = (double complex const *const *const) QpL_fIL;
			FLUXDATA->F_c = FluxViscL_fIL;
			flux_viscous_c(FLUXDATA);

			FLUXDATA->W_c = WR_fIL;
			FLUXDATA->Q_c = (double complex const *const *const) QpR_fIL;
			FLUXDATA->F_c = FluxViscR_fIL;
			flux_viscous_c(FLUXDATA);

			for (size_t i = 0; i < NfnI*d*Neq; i++)
				FluxViscNum_fIL[i] = 0.5*(FluxViscL_fIL[i]+FluxViscR_fIL[i]);

			free(FluxViscL_fIL);
			free(FluxViscR_fIL);
		} else if (strstr(PDE,"Poisson")) {
			for (size_t dim = 0; dim < d; dim++) {
				for (size_t i = 0, iMax = NfnI*Nvar; i < iMax; i++)
					FluxViscNum_fIL[iMax*dim+i] = 0.5*(QpL_fIL[dim][i] + QpR_fIL[dim][i]);
			}
		} else {
			EXIT_UNSUPPORTED;
		}
	}
	free(FLUXDATA);

	dot_with_normal_c(NfnI,Neq,n_fIL,FluxViscNum_fIL,nFluxViscNum_fIL);
	free(FluxViscNum_fIL);

	// Modify nFluxViscNum_fIL to account for boundary conditions (if necessary)
	unsigned int const BC = FACE->BC;

	if (BC % BC_STEP_SC == BC_NOSLIP_ADIABATIC) {
		// Set last component of nFluxViscNum to zero
		size_t const eq = Neq-1;
		for (size_t n = 0; n < NfnI; n++)
			nFluxViscNum_fIL[n+NfnI*eq] = 0.0;
	}
}

void add_Jacobian_scaling_FACE_c(struct S_FDATA const *const FDATA, char const imex_type, char const coef_type)
{
	if (imex_type != 'E')
		EXIT_UNSUPPORTED;

	struct S_OPERATORS_F const *const *const OPS  = (struct S_OPERATORS_F const *const *const) FDATA->OPS;
	struct S_FACE        const *const        FACE = FDATA->FACE;

	unsigned int const d        = DB.d,
	                   Neq      = DB.Neq,
	                   IndFType = FDATA->IndFType,
	                   NfnI     = OPS[IndFType]->NfnI;

	double const *const detJF_fIL = FACE->detJF_fI;

	if (coef_type == 'W' || coef_type == 'V') {
		double complex *const nFluxNum_fIL = FDATA->NFLUXDATA->nFluxNum_c;

		for (size_t eq = 0; eq < Neq; eq++) {
			size_t const IndnF = eq*NfnI;
			for (size_t n = 0; n < NfnI; n++)
				nFluxNum_fIL[IndnF+n] *= detJF_fIL[n];
		}
	} else if (coef_type == 'Q') {
		double complex *const *const nSolNum_fIL = FDATA->NFLUXDATA->nSolNum_c;

		for (size_t dim = 0; dim < d; dim++) {
			for (size_t eq = 0; eq < Neq; eq++) {
				size_t const IndnF = eq*NfnI;
				for (size_t n = 0; n < NfnI; n++)
					nSolNum_fIL[dim][IndnF+n] *= detJF_fIL[n];
			}
		}
	} else {
		EXIT_UNSUPPORTED;
	}

}

static void swap_FACE_orientation_c(struct S_FDATA const *const FDATA, char const imex_type, char const coef_type)
{
	if (!(imex_type == 'E'))
		EXIT_UNSUPPORTED;

	struct S_OPERATORS_F const *const *const OPS  = (struct S_OPERATORS_F const *const *const) FDATA->OPS;

	unsigned int const d             = DB.d,
	                   Neq           = DB.Neq,
	                   IndFType      = FDATA->IndFType,
	                   NfnI          = OPS[IndFType]->NfnI,
	                   *const nOrdLR = OPS[IndFType]->nOrdLR;


	if (coef_type == 'W' || coef_type == 'V') {
		double complex *const nFluxNum_fI = FDATA->NFLUXDATA->nFluxNum_c;

		for (size_t i = 0, iMax = Neq*NfnI; i < iMax; i++)
			nFluxNum_fI[i] *= -1.0;

		array_rearrange_cmplx(NfnI,Neq,nOrdLR,'C',nFluxNum_fI);
	} else if (coef_type == 'Q') {
		double complex *const *const nSolNum_fI = FDATA->NFLUXDATA->nSolNum_c;

		for (size_t dim = 0; dim < d; dim++) {
			for (size_t i = 0, iMax = Neq*NfnI; i < iMax; i++)
				nSolNum_fI[dim][i] *= -1.0;

			array_rearrange_cmplx(NfnI,Neq,nOrdLR,'C',nSolNum_fI[dim]);
		}
	} else {
		EXIT_UNSUPPORTED;
	}
}

void finalize_FACE_Inviscid_Weak_c(struct S_FDATA const *const FDATAL, struct S_FDATA const *const FDATAR,
                                   double complex const *const nANumL_fI, char const side, char const imex_type,
                                   char const coef_type)
{
	if (!(imex_type == 'E'))
		EXIT_UNSUPPORTED;

	struct S_FDATA const *FDATA = NULL;
	double complex       *RHS   = NULL;

	if (side == 'L') {
		FDATA = FDATAL;
		RHS   = FDATA->FACE->RHSL_c;
	} else if (side == 'R') {
		FDATA = FDATAR;
		RHS   = FDATA->FACE->RHSR_c;

		swap_FACE_orientation_c(FDATAR,imex_type,coef_type);
	} else {
		EXIT_UNSUPPORTED;
	}

	struct S_OPERATORS_F const *const *const OPS = (struct S_OPERATORS_F const *const *const) FDATA->OPS;

	unsigned int const Neq      = DB.Neq,
	                   Vf       = FDATA->Vf,
	                   IndFType = FDATA->IndFType,
	                   NfnI     = OPS[IndFType]->NfnI;

	mm_dcc(CBCM,CBT,CBNT,OPS[0]->NvnS,Neq,NfnI,1.0,1.0,OPS[0]->I_Weak_FV[Vf],nANumL_fI,RHS);
}

void finalize_QhatF_Weak_c(struct S_FDATA const *const FDATAL, struct S_FDATA const *const FDATAR, char const side,
                           char const imex_type)
{
	if (!(imex_type == 'E'))
		EXIT_UNSUPPORTED;

	struct S_FACE const *const FACE = FDATAL->FACE;

	correct_numerical_solution_strong_c(FDATAL,'E',side);

	struct S_FDATA const *FDATA = NULL;

	double complex *const *QhatF;
	if (side == 'L') {
		FDATA = FDATAL;
		QhatF = (double complex *const *const) FACE->QhatL_c;
	} else if (side == 'R') {
		FDATA = FDATAR;
		QhatF = (double complex *const *const) FACE->QhatR_c;
		swap_FACE_orientation_c(FDATAR,imex_type,'Q');
	} else {
		EXIT_UNSUPPORTED;
	}

	struct S_OPERATORS_F const *const *const OPS = (struct S_OPERATORS_F const *const *const) FDATA->OPS;

	unsigned int const d        = DB.d,
	                   Neq      = DB.Neq,
	                   Nvar     = DB.Nvar,
	                   Vf       = FDATA->Vf,
	                   IndFType = FDATA->IndFType,
	                   NfnI     = OPS[IndFType]->NfnI,
	                   NvnS     = OPS[0]->NvnS;

	double complex const *const *const nSolNum_fI = (double complex const *const *const) FDATA->NFLUXDATA->nSolNum_c;

	double complex *const InSNum = malloc(NvnS*Nvar * sizeof *InSNum); // free
	for (size_t dim = 0; dim < d; dim++) {
		// Note that there is a minus sign included in the definition of I_Weak_FV.
		mm_dcc(CBCM,CBT,CBNT,OPS[0]->NvnS,Neq,NfnI,1.0,0.0,OPS[0]->I_Weak_FV[Vf],nSolNum_fI[dim],InSNum);
		for (size_t i = 0, iMax = NvnS*Nvar; i < iMax; i++)
			QhatF[dim][i] = -InSNum[i];
	}
	free(InSNum);
}

void finalize_FACE_Viscous_Weak_c(struct S_FDATA const *const FDATAL, struct S_FDATA const *const FDATAR,
                                  double complex const *const nANumL_fI, char const side, char const imex_type,
                                  char const coef_type)
{
	finalize_FACE_Inviscid_Weak_c(FDATAL,FDATAR,nANumL_fI,side,imex_type,coef_type);
}
