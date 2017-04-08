// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

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
	 *		It is currently hard-coded that GradW is of the same order as the solution.
	 *		Note, if Collocation is enable, that D_Weak includes the inverse cubature weights.
	 *
	 *	References:
	 *		Zwanenburg(2016)-Equivalence_between_the_Energy_Stable_Flux_Reconstruction_and_Discontinuous_Galerkin_Schemes
	 */

	// Initialize DB Parameters
	unsigned int d          = DB.d,
	             Neq        = DB.Neq,
	             Collocated = DB.Collocated;

	// Standard datatypes
	struct S_OPERATORS_V *OPS;
	struct S_Dxyz        *DxyzInfo;

	OPS      = malloc(sizeof *OPS);      // free
	DxyzInfo = malloc(sizeof *DxyzInfo); // free

	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		unsigned int NvnS, NvnI;
		double       **DxyzChiS;

		init_ops_VOLUME(OPS,VOLUME,0);
//		init_ops_VOLUME(OPS[0],VOLUME,0);
//		if (VOLUME->type == WEDGE)
//			init_ops_VOLUME(OPS[1],VOLUME,1);

		NvnS = OPS->NvnS;
		NvnI = OPS->NvnI;

		DxyzInfo->Nbf = OPS->NvnS;
		DxyzInfo->Nn  = OPS->NvnI;
		DxyzInfo->D   = (double const *const *const) OPS->D_Weak;
		DxyzInfo->C   = VOLUME->C_vI;

		double const *const ChiS_vI = OPS->ChiS_vI;

		DxyzChiS = VOLUME->DxyzChiS;
		for (size_t dim = 0; dim < d; dim++) {
			double *Dxyz;

			DxyzInfo->dim = dim;
			Dxyz = compute_Dxyz(DxyzInfo,d); // free/keep

			// Note: The detJ_vI term cancels with the gradient operator (Zwanenburg(2016), eq. (B.2))
			if (Collocated) { // ChiS_vI == I
				DxyzChiS[dim] = Dxyz;
			} else {
				DxyzChiS[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS,NvnS,NvnI,1.0,Dxyz,ChiS_vI); // keep
				free(Dxyz);
			}
// Might not need to store DxyzChiS for explicit (ToBeDeleted)

			// Compute intermediate (see comments) Qhat contribution
			VOLUME->Qhat[dim] = mm_Alloc_d(CBCM,CBT,CBNT,NvnS,Neq,NvnS,1.0,DxyzChiS[dim],VOLUME->What); // keep
		}
	}

	free(OPS);
	free(DxyzInfo);
}

static void boundary_NavierStokes(const unsigned int Nn, const unsigned int Nel, const double *XYZ, const double *nL,
                                  const double *WL, double *WB, const unsigned int d, const unsigned int Nvar,
                                  const unsigned int BC)
{
	if (BC % BC_STEP_SC == BC_DIRICHLET) {
		boundary_NoSlip_Dirichlet(Nn,Nel,XYZ,WL,WB,nL,d,Nvar);
	} else if (BC % BC_STEP_SC == BC_NOSLIP_ADIABATIC) {
		boundary_NoSlip_Adiabatic(Nn,Nel,XYZ,WL,WB,nL,d,Nvar);
	} else {
		EXIT_UNSUPPORTED;
	}
}

static void explicit_GradW_FACE(void)
{
	/*
	 *	Purpose:
	 *		Compute intermediate FACE contribution to Qhat.
	 *
	 *	Comments:
	 *		L/R is used in place of In/Out (the previous convention).
	 *		This is an intermediate contribution because the multiplication by MInv is not included.
	 *		It is currently hard-coded that GradW is of the same order as the solution and that a central numerical flux
	 *		is used.
	 *		Note, if Collocation is enable, that I_Weak includes the inverse cubature weights.
	 */

	// Initialize DB Parameters
	unsigned int d    = DB.d,
	             Nvar = DB.Nvar;

	// Standard datatypes
	struct S_OPERATORS_F *OPSL[2], *OPSR[2];
	struct S_FDATA       *FDATAL, *FDATAR;
	FDATAL = malloc(sizeof *FDATAL); // free
	FDATAR = malloc(sizeof *FDATAR); // free
	FDATAL->OPS = OPSL;
	FDATAR->OPS = OPSR;

	for (size_t i = 0; i < 2; i++) {
		OPSL[i]  = malloc(sizeof *OPSL[i]);  // free
		OPSR[i]  = malloc(sizeof *OPSR[i]);  // free
	}

	for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
		// Load required FACE operators (Left and Right FACEs)
		init_FDATA(FDATAL,FACE,'L');
		init_FDATA(FDATAR,FACE,'R');

		struct S_VOLUME const *const VL = FDATAL->VOLUME,
		                      *const VR = FDATAR->VOLUME;

		unsigned int NfnI, NvnSL, NvnSR, VfL, VfR;
		NfnI  = OPSL[0]->NfnI;
		NvnSL = OPSL[0]->NvnS;
		NvnSR = OPSR[0]->NvnS;

		VfL   = FDATAL->Vf;
		VfR   = FDATAR->Vf;

		unsigned int Boundary = FACE->Boundary;
		unsigned int BC       = FACE->BC;

		// Compute WL_fIL
		double *WL_fIL = malloc(NfnI*Nvar * sizeof *WL_fIL); // free
		FDATAL->W_fIL = WL_fIL;
		coef_to_values_fI(FDATAL,'W');

		// Compute WR_fIL (Taking BCs into account if applicable)
		double *ChiSR_fIL = NULL;
		double *n_fI = FACE->n_fI;

		double *WR_fIL = malloc(NfnI*Nvar * sizeof *WR_fIL); // free
		if (!Boundary) {
			unsigned int const *nOrdRL = OPSL[FDATAL->IndFType]->nOrdRL;

			// Use reordered operator (as opposed to solution) for ease of linearization.
			ChiSR_fIL = malloc(NfnI*NvnSR * sizeof *ChiSR_fIL); // free

			double const *ChiS_fI = OPSR[0]->ChiS_fI[VfR];
			for (size_t i = 0; i < NfnI; i++) {
			for (size_t j = 0; j < NvnSR; j++) {
				ChiSR_fIL[i*NvnSR+j] = ChiS_fI[nOrdRL[i]*NvnSR+j];
			}}

			// Could be performed with sum factorization using a reordering of the sum factorized FACE operator.
			mm_CTN_d(NfnI,Nvar,NvnSR,ChiSR_fIL,VR->What,WR_fIL);
		} else {
			boundary_NavierStokes(NfnI,1,FACE->XYZ_fI,n_fI,WL_fIL,WR_fIL,d,Nvar,BC);
		}

		// Compute normal numerical solution (nWnum_fI == n_fI[dim] (dot) 0.5*(WL_fI+WR_fI))
		double **nWnum_fI = malloc(d * sizeof *nWnum_fI); // free
		for (size_t dim = 0; dim < d; dim++)
			nWnum_fI[dim] = malloc(NfnI*Nvar * sizeof *nWnum_fI[dim]); // free

		for (size_t dim = 0; dim < d; dim++) {
		for (size_t var = 0; var < Nvar; var++) {
		for (size_t n = 0; n < NfnI; n++) {
			nWnum_fI[dim][var*NfnI+n] = n_fI[dim*NfnI+n]*0.5*(WL_fIL[var*NfnI+n]+WR_fIL[var*NfnI+n]);
		}}}
		free(WL_fIL);
		free(WR_fIL);

		// Add in FACE Jacobian determinant term
		double *detJF_fI = FACE->detJF_fI;

		double **JnWnum_fI = nWnum_fI;
		for (size_t dim = 0; dim < d; dim++) {
		for (size_t var = 0; var < Nvar; var++) {
		for (size_t n = 0; n < NfnI; n++) {
			JnWnum_fI[dim][n+var*NfnI] *= detJF_fI[n];
		}}}

		// Compute intermediate Qhat contributions

		// Interior VOLUME
		for (size_t dim = 0; dim < d; dim++) {
			// Note that there is a minus sign included in the definition of I_Weak_FF.
			mm_d(CBCM,CBT,CBNT,NvnSL,Nvar,NfnI,-1.0,1.0,OPSL[0]->I_Weak_FF[VfL],JnWnum_fI[dim],VL->Qhat[dim]);
		}

		// Exterior VOLUME
		if (!Boundary) {
			unsigned int const *nOrdLR = OPSL[FDATAL->IndFType]->nOrdLR;
			for (size_t dim = 0; dim < d; dim++) {
				array_rearrange_d(NfnI,Nvar,nOrdLR,'C',JnWnum_fI[dim]);

				// minus sign from using negative normal for the opposite VOLUME cancels with minus sign in I_Weak_FF.
				mm_d(CBCM,CBT,CBNT,NvnSR,Nvar,NfnI,1.0,1.0,OPSR[0]->I_Weak_FF[VfR],nWnum_fI[dim],VR->Qhat[dim]);
			}
			free(ChiSR_fIL);
		}
		array_free2_d(d,nWnum_fI);
	}

	for (size_t i = 0; i < 2; i++) {
		free(OPSL[i]);
		free(OPSR[i]);
	}
	free(FDATAL);
	free(FDATAR);
}

static void explicit_GradW_finalize(void)
{
	/*
	 *	Purpose:
	 *		Add in inverse mass matrix contribution to VOLUME->Qhat.
	 *
	 *	Comments:
	 *		The FACE contribution were included direcly in VL/VR->Qhat.
	 *		If Collocation is enable, only the inverse Jacobian determinant is missing.
	 *
	 */

	unsigned int d          = DB.d,
	             Nvar       = DB.Nvar,
	             Collocated = DB.Collocated;

	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		unsigned int NvnS = VOLUME->NvnS;
		double       *Qhat;

		if (Collocated) {
			double *detJV_vI = VOLUME->detJV_vI;

			for (size_t dim = 0; dim < d; dim++) {
				Qhat = VOLUME->Qhat[dim];
				for (size_t var = 0; var < Nvar; var++) {
				for (size_t n = 0; n < NvnS; n++) {
					Qhat[var*NvnS+n] /= detJV_vI[n];
				}}
			}
		} else {

			if (VOLUME->MInv == NULL)
				compute_inverse_mass(VOLUME);

			for (size_t dim = 0; dim < d; dim++) {
				Qhat = malloc(NvnS*Nvar * sizeof *Qhat); // keep (free previously stored Qhat)

				mm_CTN_d(NvnS,Nvar,NvnS,VOLUME->MInv,VOLUME->Qhat[dim],Qhat);
				free(VOLUME->Qhat[dim]);
				VOLUME->Qhat[dim] = Qhat;
			}
		}
	}
}
