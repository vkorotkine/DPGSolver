// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#include "Macros.h"
#include "Parameters.h"
#include "S_DB.h"
#include "S_VOLUME.h"
#include "S_FACE.h"

#include "solver_functions.h"
#include "matrix_functions.h"
#include "array_free.h"
#include "array_swap.h"
#include "boundary_conditions.h"
#include "update_VOLUMEs.h"

/*
 *	Purpose:
 *		Compute weak gradients required for the computation of viscous fluxes.
 *
 *	Comments:
 *		The weak form of the equation is used (i.e. integrated by parts once).
 *
 *	Notation:
 *
 *	References:
 */

static void explicit_GradW_VOLUME   (void);
static void explicit_GradW_FACE     (void);
static void explicit_GradW_finalize (void);

void explicit_GradW(void)
{
	explicit_GradW_VOLUME();
	explicit_GradW_FACE();
	explicit_GradW_finalize();
}

struct S_Dxyz {
	unsigned int dim, Nbf, Nn;
	double const *const *D, *C;
};


static double *compute_Dxyz(struct S_Dxyz *DxyzInfo, unsigned int d)
{
	/*
	 *	Purpose:
	 *		Compute physical derivative operator matrices using the chain rule.
	 *
	 *	Comments:
	 *		Note the ordering of C specified in the comments of setup_geom_factors.
	 *
	 *	References:
	 *		Zwanenburg(2016)-Equivalence_between_the_Energy_Stable_Flux_Reconstruction_and_Discontinuous_Galerkin_
	 *		                 Schemes (eq. B.2)
	 */

	unsigned int Nbf, Nn, dim1, IndC;
	double       *Dxyz;

	Nbf  = DxyzInfo->Nbf;
	Nn   = DxyzInfo->Nn;
	dim1 = DxyzInfo->dim;

	double const *const *const D = DxyzInfo->D;
	double const *const        C = DxyzInfo->C;

	Dxyz = calloc(Nbf*Nn, sizeof *Dxyz); // keep
	for (size_t dim2 = 0; dim2 < d; dim2++) {
		IndC = (dim1+dim2*d)*Nn;
		mm_diag_d(Nbf,Nn,&C[IndC],D[dim2],Dxyz,1.0,1.0,'R','R');
	}

	return Dxyz;
}

static void explicit_GradW_VOLUME(void)
{
	/*
	 *	Purpose:
	 *		Compute intermediate VOLUME contribution to Qhat.
	 *
	 *	Comments:
	 *		This is an intermediate contribution because the multiplication by MInv is not included.
	 *		The contribution from this function is duplicated in QhatV as the local contribution is required for the
	 *		numerical flux in the second equation of the mixed form.
	 *		It is currently hard-coded that GradW is of the same order as the solution.
	 *		Note, if Collocation is enable, that D_Weak includes the inverse cubature weights.
	 *
	 *	References:
	 *		Zwanenburg(2016)-Equivalence_between_the_Energy_Stable_Flux_Reconstruction_and_Discontinuous_Galerkin_Schemes
	 */

	// Initialize DB Parameters
	unsigned int d    = DB.d,
	             Nvar = d+2;

	// Standard datatypes
	struct S_OPERATORS_V *OPS[2];

	struct S_VDATA *VDATA = malloc(sizeof *VDATA); // free
	VDATA->OPS = (struct S_OPERATORS_V const *const *) OPS;

	for (size_t i = 0; i < 2; i++)
		OPS[i] = malloc(sizeof *OPS[i]); // free

	struct S_Dxyz *const DxyzInfo = malloc(sizeof *DxyzInfo); // free

	if (0) { // Weak form for 1st equation in the mixed form
		for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			init_VDATA(VDATA,VOLUME);

			unsigned int const NvnS = VDATA->OPS[0]->NvnS,
			                   NvnI = VDATA->OPS[0]->NvnI;

			DxyzInfo->Nbf = VDATA->OPS[0]->NvnS;
			DxyzInfo->Nn  = VDATA->OPS[0]->NvnI;
			DxyzInfo->D   = (double const *const *const) VDATA->OPS[0]->D_Weak;
			DxyzInfo->C   = VOLUME->C_vI;

			double const *const        ChiS_vI  = VDATA->OPS[0]->ChiS_vI;
			double       *      *const DxyzChiS = VOLUME->DxyzChiS;
			for (size_t dim = 0; dim < d; dim++) {
				DxyzInfo->dim = dim;
				double *const Dxyz = compute_Dxyz(DxyzInfo,d); // free

				// Note: The detJ_vI term cancels with the gradient operator (Zwanenburg(2016), eq. (B.2))
				if (DB.Collocated) { // ChiS_vI == I
					DxyzChiS[dim] = Dxyz;
				} else {
					DxyzChiS[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS,NvnS,NvnI,1.0,Dxyz,ChiS_vI); // free
					free(Dxyz);
				}

				// Compute intermediate Qhat contribution
				if (VOLUME->QhatV[dim])
					free(VOLUME->QhatV[dim]);
				VOLUME->QhatV[dim] = malloc(NvnS*Nvar * sizeof *(VOLUME->QhatV[dim])); // keep
				mm_d(CBCM,CBT,CBNT,NvnS,Nvar,NvnS,-1.0,0.0,DxyzChiS[dim],VOLUME->What,VOLUME->QhatV[dim]);

				if (VOLUME->Qhat[dim])
					free(VOLUME->Qhat[dim]);
				VOLUME->Qhat[dim]  = malloc(NvnS*Nvar * sizeof *(VOLUME->Qhat[dim]));  // keep (used below)

				for (size_t i = 0; i < NvnS*Nvar; i++)
					VOLUME->Qhat[dim][i] = VOLUME->QhatV[dim][i];


				// Need to store DxyzChiS for implicit runs. (Don't forget -ve sign) (ToBeDeleted)
				free(DxyzChiS[dim]);
			}
		}
	} else { // Strong form
		for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			init_VDATA(VDATA,VOLUME);

			unsigned int const NvnS = VDATA->OPS[0]->NvnS,
			                   NvnI = VDATA->OPS[0]->NvnI;

			DxyzInfo->Nbf = VDATA->OPS[0]->NvnS;
			DxyzInfo->Nn  = VDATA->OPS[0]->NvnI;
			DxyzInfo->D   = (double const *const *const) VDATA->OPS[0]->D_Weak;
			DxyzInfo->C   = VOLUME->C_vI;

			double const *const        ChiS_vI  = VDATA->OPS[0]->ChiS_vI;
			double       *      *const DxyzChiS = VOLUME->DxyzChiS;
			for (size_t dim = 0; dim < d; dim++) {
				DxyzInfo->dim = dim;
				double *const Dxyz = compute_Dxyz(DxyzInfo,d); // free

				// Note: The detJ_vI term cancels with the gradient operator (Zwanenburg(2016), eq. (B.2))
				if (DB.Collocated) { // ChiS_vI == I
					DxyzChiS[dim] = Dxyz;
				} else {
					DxyzChiS[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS,NvnS,NvnI,1.0,Dxyz,ChiS_vI); // free
					free(Dxyz);
				}

				// Compute intermediate Qhat contribution
				if (VOLUME->QhatV[dim])
					free(VOLUME->QhatV[dim]);
				VOLUME->QhatV[dim] = malloc(NvnS*Nvar * sizeof *(VOLUME->QhatV[dim])); // keep
				// Note: Using CBCM with CBNT for DxyzChiS (stored in row-major ordering) gives DxyzChiS' in the
				//       operation below.
				mm_d(CBCM,CBNT,CBNT,NvnS,Nvar,NvnS,1.0,0.0,DxyzChiS[dim],VOLUME->What,VOLUME->QhatV[dim]);

				if (VOLUME->Qhat[dim])
					free(VOLUME->Qhat[dim]);
				VOLUME->Qhat[dim]  = malloc(NvnS*Nvar * sizeof *(VOLUME->Qhat[dim]));  // keep (used below)

				for (size_t i = 0; i < NvnS*Nvar; i++)
					VOLUME->Qhat[dim][i] = VOLUME->QhatV[dim][i];

				// Need to store DxyzChiS for implicit runs. (Don't forget -ve sign) (ToBeDeleted)
				free(DxyzChiS[dim]);
			}
		}
	}

	free(VDATA);
	for (size_t i = 0; i < 2; i++)
		free(OPS[i]);
	free(DxyzInfo);
}

static void explicit_GradW_FACE(void)
{
	/*
	 *	Purpose:
	 *		Compute intermediate FACE contribution to Qhat.
	 *
	 *	Comments:
	 *		L/R is used in place of In/Out (the previous convention). (ToBeDeleted)
	 *		This is an intermediate contribution because the multiplication by MInv is not included.
	 *		It is currently hard-coded that GradW is of the same order as the solution and that a central numerical flux
	 *		is used.
	 *		Note, if Collocation is enable, that I_Weak includes the inverse cubature weights.
	 */

	// Initialize DB Parameters
	unsigned int d    = DB.d,
	             Nvar = d+2,
	             Neq  = d+2;

	// Standard datatypes
	struct S_OPERATORS_F *OPSL[2], *OPSR[2];
	struct S_FDATA       *FDATAL = malloc(sizeof *FDATAL), // free
	                     *FDATAR = malloc(sizeof *FDATAR); // free
	FDATAL->OPS = (struct S_OPERATORS_F const *const *) OPSL;
	FDATAR->OPS = (struct S_OPERATORS_F const *const *) OPSR;

	struct S_NumericalFlux *const NFluxData = malloc(sizeof *NFluxData); // free
	FDATAL->NFluxData = NFluxData;
	FDATAR->NFluxData = NFluxData;

	for (size_t i = 0; i < 2; i++) {
		OPSL[i]  = malloc(sizeof *OPSL[i]);  // free
		OPSR[i]  = malloc(sizeof *OPSR[i]);  // free
	}

	if (0) { // Weak form for 1st equation in the mixed form
		for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
			init_FDATA(FDATAL,FACE,'L');
			init_FDATA(FDATAR,FACE,'R');

			// Compute WL_fIL and WR_fIL (i.e. as seen from the (L)eft VOLUME)
			unsigned int const IndFType = FDATAL->IndFType,
			                   NfnI     = OPSL[IndFType]->NfnI;

			FDATAL->W_fIL = malloc(NfnI*Nvar * sizeof *(FDATAL->W_fIL)), // free
			FDATAR->W_fIL = malloc(NfnI*Nvar * sizeof *(FDATAR->W_fIL)); // free

			coef_to_values_fI(FDATAL,'W');
			compute_WR_fIL(FDATAR,FDATAL->W_fIL,FDATAR->W_fIL);

			// Compute numerical flux as seen from the left VOLUME
			NFluxData->WL_fIL     = FDATAL->W_fIL;
			NFluxData->WR_fIL     = FDATAR->W_fIL;
			NFluxData->nSolNum_fI = malloc(d * sizeof *(NFluxData->nSolNum_fI)); // free
			for (size_t dim = 0; dim < d; dim++)
				NFluxData->nSolNum_fI[dim] = malloc(NfnI*Neq * sizeof *(NFluxData->nSolNum_fI[dim])); // free

			compute_numerical_solution(FDATAL,'E','W');
			add_Jacobian_scaling_FACE(FDATAL,'E','Q');

			free(FDATAL->W_fIL);
			free(FDATAR->W_fIL);

			// Add in memory allocation here. (ToBeDeleted)
			EXIT_UNSUPPORTED;

			finalize_QhatF_Weak(FDATAL,FDATAR,'L','E');
			if (!FACE->Boundary)
				finalize_QhatF_Weak(FDATAL,FDATAR,'R','E');

			array_free2_d(d,NFluxData->nSolNum_fI);
		}
	} else { // Strong form
		for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
			init_FDATA(FDATAL,FACE,'L');
			init_FDATA(FDATAR,FACE,'R');

			// Compute WL_fIL and WR_fIL (i.e. as seen from the (L)eft VOLUME)
			unsigned int const IndFType = FDATAL->IndFType,
			                   NfnI     = OPSL[IndFType]->NfnI;

			FDATAL->W_fIL = malloc(NfnI*Nvar * sizeof *(FDATAL->W_fIL)), // free
			FDATAR->W_fIL = malloc(NfnI*Nvar * sizeof *(FDATAR->W_fIL)); // free

			coef_to_values_fI(FDATAL,'W');
			compute_WR_fIL(FDATAR,FDATAL->W_fIL,FDATAR->W_fIL);

			// Compute numerical flux as seen from the left VOLUME
			NFluxData->WL_fIL     = FDATAL->W_fIL;
			NFluxData->WR_fIL     = FDATAR->W_fIL;
			NFluxData->nSolNum_fI = malloc(d * sizeof *(NFluxData->nSolNum_fI)); // free
			for (size_t dim = 0; dim < d; dim++)
				NFluxData->nSolNum_fI[dim] = malloc(NfnI*Neq * sizeof *(NFluxData->nSolNum_fI[dim])); // free

			compute_numerical_solution(FDATAL,'E','S');
			add_Jacobian_scaling_FACE(FDATAL,'E','Q');

			free(FDATAL->W_fIL);
			free(FDATAR->W_fIL);

			for (size_t dim = 0; dim < d; dim++) {
				if (FACE->QhatL[dim])
					free(FACE->QhatL[dim]);
				FACE->QhatL[dim] = malloc(OPSL[0]->NvnS*Nvar * sizeof *(FACE->QhatL[dim])); // keep
			}
			finalize_QhatF_Weak(FDATAL,FDATAR,'L','E');

			if (!FACE->Boundary) {
				for (size_t dim = 0; dim < d; dim++) {
					if (FACE->QhatR[dim])
						free(FACE->QhatR[dim]);
					FACE->QhatR[dim] = malloc(OPSR[0]->NvnS*Nvar * sizeof *(FACE->QhatR[dim])); // keep
				}
				finalize_QhatF_Weak(FDATAL,FDATAR,'R','E');
			}

			array_free2_d(d,NFluxData->nSolNum_fI);
		}
	}

	for (size_t i = 0; i < 2; i++) {
		free(OPSL[i]);
		free(OPSR[i]);
	}
	free(NFluxData);
	free(FDATAL);
	free(FDATAR);
}

static void finalize_Qhat(struct S_VOLUME const *const VOLUME, unsigned int const NvnS, double *const *const Qhat)
{
	unsigned int const d    = DB.d,
	                   Nvar = d+2;

	if (DB.Collocated) {
		double const *const detJV_vI = VOLUME->detJV_vI;
		for (size_t dim = 0; dim < d; dim++) {
			for (size_t var = 0; var < Nvar; var++) {
			for (size_t n = 0; n < NvnS; n++) {
				Qhat[dim][var*NvnS+n] /= detJV_vI[n];
			}}
		}
	} else {
		double *Qhat_tmp = malloc(NvnS*Nvar * sizeof *Qhat_tmp); // free
		for (size_t dim = 0; dim < d; dim++) {
			mm_CTN_d(NvnS,Nvar,NvnS,VOLUME->MInv,Qhat[dim],Qhat_tmp);
			for (size_t var = 0; var < Nvar; var++) {
			for (size_t n = 0; n < NvnS; n++) {
				Qhat[dim][var*NvnS+n] = Qhat_tmp[var*NvnS+n];
			}}
		}
		free(Qhat_tmp);
	}
}

static void explicit_GradW_finalize(void)
{
	/*
	 *	Purpose:
	 *		Add in inverse mass matrix contribution to VOLUME->Qhat, VOLUME->QhatV.
	 *
	 *	Comments:
	 *		All contributions continue to be stored individually as they must be used as such for the computation of the
	 *		viscous numerical flux. The FACE Qhat contributions are added to the VOLUME Qhat contribution to store the
	 *		entire weak gradient in VOLUME->Qhat (used for VOLUME terms).
	 *
	 *		The FACE contribution were included direcly in VL/VR->Qhat while the VOLUME contributions were stored in
	 *		QhatV. The two are now summed in VL/VR->Qhat and QhatV is retained for use in the numerical flux for the
	 *		second equation of the mixed formulation.
	 *		If Collocation is enable, only the inverse Jacobian determinant is missing.
	 */

	unsigned int const d    = DB.d,
	                   Nvar = d+2;

	// Add FACE contributions to VOLUME->Qhat then multiply by MInv
	for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
		struct S_VOLUME const *const VL = FACE->VIn,
		                      *const VR = FACE->VOut;

		unsigned int const NvnSL = VL->NvnS,
		                   NvnSR = VR->NvnS;

		for (size_t dim = 0; dim < d; dim++) {
			for (size_t i = 0; i < NvnSL*Nvar; i++)
				VL->Qhat[dim][i] += FACE->QhatL[dim][i];

			for (size_t i = 0; i < NvnSR*Nvar; i++)
				VR->Qhat[dim][i] += FACE->QhatR[dim][i];
		}

		finalize_Qhat(VL,NvnSL,FACE->QhatL);
		finalize_Qhat(VR,NvnSR,FACE->QhatR);
	}

	// Multiply VOLUME Qhat terms by MInv
	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		unsigned int const NvnS = VOLUME->NvnS;

		finalize_Qhat(VOLUME,NvnS,VOLUME->Qhat);
		finalize_Qhat(VOLUME,NvnS,VOLUME->QhatV);
	}
}
