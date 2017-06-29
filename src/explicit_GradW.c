// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "explicit_GradW.h"

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
#include "fluxes_structs.h"
#include "array_free.h"
#include "array_print.h"

/*
 *	Purpose:
 *		Compute weak gradients required for the computation of viscous fluxes.
 *
 *	Comments:
 *		The use of the doubly integrated by parts weak form for the first equation in the mixed formulation is much more
 *		consistent with the theoretical formulation of the stabilized mixed form as presented in Brezzi(2000) (and in
 *		general); the stabilization is a penalization on the solution jumps across the elements. This also gives a much
 *		more natural symmetry to the diffusive operator when compared with using the weak form for the first equation.
 *		Furthermore, the partially corrected gradient contributions required for the numerical viscous flux can be
 *		directly computed from the summation of the VOLUME and FACE terms to Qhat in this case. Note that the doubly
 *		integrated by parts weak form was previously called the strong form.
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
	if (!DB.Viscous)
		return;

	explicit_GradW_VOLUME();
	explicit_GradW_FACE();
	explicit_GradW_finalize();
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

	unsigned int const d    = DB.d,
	                   Nvar = d+2;

	struct S_OPERATORS_V *OPS[2];

	struct S_VDATA *const VDATA = malloc(sizeof *VDATA); // free
	VDATA->OPS = (struct S_OPERATORS_V const *const *) OPS;

	for (size_t i = 0; i < 2; i++)
		OPS[i] = malloc(sizeof *OPS[i]); // free

	struct S_Dxyz *const DxyzInfo = malloc(sizeof *DxyzInfo); // free

	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		init_VDATA(VDATA,VOLUME);

		unsigned int const NvnS = VDATA->OPS[0]->NvnS,
		                   NvnI = VDATA->OPS[0]->NvnI;

		DxyzInfo->Nbf = VDATA->OPS[0]->NvnS;
		DxyzInfo->Nn  = VDATA->OPS[0]->NvnI;
		DxyzInfo->D   = (double const *const *const) VDATA->OPS[0]->D_Weak;
		DxyzInfo->C   = VOLUME->C_vI;

		double const *const        ChiS_vI  = VDATA->OPS[0]->ChiS_vI;
		for (size_t dim = 0; dim < d; dim++) {
			DxyzInfo->dim = dim;
			double *const Dxyz = compute_Dxyz(DxyzInfo,d); // free

			// Note: The detJ_vI term cancels with the gradient operator (Zwanenburg(2016), eq. (B.2))
			double *DxyzChiS = NULL;
			if (DB.Collocated) { // ChiS_vI == I
				DxyzChiS = Dxyz;
			} else {
				DxyzChiS = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS,NvnS,NvnI,1.0,Dxyz,ChiS_vI); // free
				free(Dxyz);
			}

			// Compute intermediate Qhat contribution

			// Note: Using CBCM with CBNT for DxyzChiS (stored in row-major ordering) gives DxyzChiS transposed in the
			//       operation below.
			mm_d(CBCM,CBNT,CBNT,NvnS,Nvar,NvnS,1.0,0.0,DxyzChiS,VOLUME->What,VOLUME->QhatV[dim]);
			free(DxyzChiS);

			for (size_t i = 0; i < NvnS*Nvar; i++)
				VOLUME->Qhat[dim][i] = VOLUME->QhatV[dim][i];
if (VOLUME->indexg == 0 && dim == 0)
array_print_d(NvnS,Nvar,VOLUME->Qhat[dim],'C');
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

	struct S_OPERATORS_F *OPSL[2], *OPSR[2];
	struct S_FDATA       *const FDATAL = malloc(sizeof *FDATAL), // free
	                     *const FDATAR = malloc(sizeof *FDATAR); // free
	FDATAL->OPS = (struct S_OPERATORS_F const *const *) OPSL;
	FDATAR->OPS = (struct S_OPERATORS_F const *const *) OPSR;

	struct S_NUMERICALFLUX *const NFLUXDATA = malloc(sizeof *NFLUXDATA); // free
	FDATAL->NFLUXDATA = NFLUXDATA;
	FDATAR->NFLUXDATA = NFLUXDATA;

	struct S_DATA *const DATA = malloc(sizeof *DATA); // free
	DATA->FDATAL    = FDATAL;
	DATA->FDATAR    = FDATAR;
	DATA->NFLUXDATA = NFLUXDATA;
	DATA->feature   = 'F';
	DATA->imex_type = 'E';

	for (size_t i = 0; i < 2; i++) {
		OPSL[i]  = malloc(sizeof *OPSL[i]);  // free
		OPSR[i]  = malloc(sizeof *OPSR[i]);  // free
	}

	for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
		init_FDATA(FDATAL,FACE,'L');
		init_FDATA(FDATAR,FACE,'R');

		// Compute WL_fIL and WR_fIL (i.e. as seen from the (L)eft VOLUME)
		manage_solver_memory(DATA,'A','W'); // free

		coef_to_values_fI(FDATAL,'W','E');
		compute_WR_fIL(FDATAR,FDATAL->W_fIL,FDATAR->W_fIL);

		// Compute numerical flux as seen from the left VOLUME
		NFLUXDATA->WL = FDATAL->W_fIL;
		NFLUXDATA->WR = FDATAR->W_fIL;
		manage_solver_memory(DATA,'A','S'); // free

		compute_numerical_solution(FDATAL,'E');
		add_Jacobian_scaling_FACE(FDATAL,'E','Q');

		finalize_QhatF_Weak(FDATAL,FDATAR,'L','E');
		if (!FACE->Boundary)
			finalize_QhatF_Weak(FDATAL,FDATAR,'R','E');

		manage_solver_memory(DATA,'F','W');
		manage_solver_memory(DATA,'F','S');
	}

	for (size_t i = 0; i < 2; i++) {
		free(OPSL[i]);
		free(OPSR[i]);
	}
	free(NFLUXDATA);
	free(FDATAL);
	free(FDATAR);
	free(DATA);
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
		double *const Qhat_tmp = malloc(NvnS*Nvar * sizeof *Qhat_tmp); // free
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
		struct S_VOLUME const *const VL = FACE->VL,
		                      *const VR = FACE->VR;

		unsigned int const NvnSL = VL->NvnS,
		                   NvnSR = VR->NvnS;

		for (size_t dim = 0; dim < d; dim++) {
			for (size_t i = 0; i < NvnSL*Nvar; i++)
				VL->Qhat[dim][i] += FACE->QhatL[dim][i];

			if (!FACE->Boundary) {
				for (size_t i = 0; i < NvnSR*Nvar; i++)
					VR->Qhat[dim][i] += FACE->QhatR[dim][i];
			}
		}

		finalize_Qhat(VL,NvnSL,FACE->QhatL);
		if (!FACE->Boundary)
			finalize_Qhat(VR,NvnSR,FACE->QhatR);
	}

	// Multiply VOLUME Qhat terms by MInv
	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		unsigned int const NvnS = VOLUME->NvnS;

		finalize_Qhat(VOLUME,NvnS,VOLUME->Qhat);
		finalize_Qhat(VOLUME,NvnS,VOLUME->QhatV);
	}
}
