// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#include "Macros.h"
#include "Parameters.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "S_FACE.h"

#include "element_functions.h"
#include "matrix_functions.h"
#include "sum_factorization.h"
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


struct S_OPERATORS {
	// VOLUME
	unsigned int NvnS, NvnI;
	double       **ChiS_vI, **D_Weak;

	// FACE
	unsigned int NfnI, NvnS_SF, NvnI_SF, NfnI_SF, *nOrdLR, *nOrdRL;
	double       **ChiS_fI, **ChiS_fI_SF, **ChiS_vI_SF,
	             **I_Weak_FF;

	struct S_OpCSR **ChiS_fI_sp;
};

struct S_FDATA {
	unsigned int P, Vf, f, SpOp, Eclass, IndFType;

	struct S_OPERATORS **OPS;
	struct S_VOLUME    *VOLUME;
};

struct S_Dxyz {
	unsigned int dim, Nbf, Nn;
	double       **D, *C;
};

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME)
{
	// Standard datatypes
	unsigned int P, type, curved;

	struct S_ELEMENT *ELEMENT;

	P      = VOLUME->P;
	type   = VOLUME->type;
	curved = VOLUME->curved;

	ELEMENT = get_ELEMENT_type(type);

	OPS->NvnS = ELEMENT->NvnS[P];
	if (!curved) {
		OPS->NvnI = ELEMENT->NvnIs[P];

		OPS->ChiS_vI = ELEMENT->ChiS_vIs[P][P];
		OPS->D_Weak  = ELEMENT->Ds_Weak_VV[P][P][0];
	} else {
		OPS->NvnI = ELEMENT->NvnIc[P];

		OPS->ChiS_vI = ELEMENT->ChiS_vIc[P][P];
		OPS->D_Weak  = ELEMENT->Dc_Weak_VV[P][P][0];
	}
}

static void init_opsF(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const struct S_FACE *FACE,
                      const unsigned int IndFType)
{
	// Initialize DB Parameters
	unsigned int ***SF_BE = DB.SF_BE;

	// Standard datatypes
	unsigned int PV, PF, Vtype, Eclass, FtypeInt, IndOrdRL, IndOrdLR;

	struct S_ELEMENT *ELEMENT, *ELEMENT_FACE, *ELEMENT_SF;

	// silence
	ELEMENT_SF = NULL;

	PV      = VOLUME->P;
	PF      = FACE->P;
	Vtype   = VOLUME->type;
	Eclass  = VOLUME->Eclass;

	FtypeInt = FACE->typeInt;
	IndOrdRL = FACE->IndOrdOutIn;
	IndOrdLR = FACE->IndOrdInOut;

	ELEMENT      = get_ELEMENT_type(Vtype);
	ELEMENT_FACE = get_ELEMENT_FACE(Vtype,IndFType);
	if ((Eclass == C_TP && SF_BE[PF][0][1]) || (Eclass == C_WEDGE && SF_BE[PF][1][1]))
		ELEMENT_SF = ELEMENT->ELEMENTclass[IndFType];
	else
		ELEMENT_SF = ELEMENT;

	OPS->NvnS    = ELEMENT->NvnS[PV];
	OPS->NvnS_SF = ELEMENT_SF->NvnS[PV];
	if (FtypeInt == 's') {
		// Straight FACE Integration
		OPS->NfnI    = ELEMENT->NfnIs[PF][IndFType];
		OPS->NfnI_SF = ELEMENT_SF->NfnIs[PF][0];
		OPS->NvnI_SF = ELEMENT_SF->NvnIs[PF];

		OPS->nOrdLR = ELEMENT_FACE->nOrd_fIs[PF][IndOrdLR];
		OPS->nOrdRL = ELEMENT_FACE->nOrd_fIs[PF][IndOrdRL];

		OPS->ChiS_fI    = ELEMENT->ChiS_fIs[PV][PF];
		OPS->ChiS_fI_sp = ELEMENT->ChiS_fIs_sp[PV][PF];
		OPS->ChiS_fI_SF = ELEMENT_SF->ChiS_fIs[PV][PF];
		OPS->ChiS_vI_SF = ELEMENT_SF->ChiS_vIs[PV][PF];

		OPS->I_Weak_FF = ELEMENT_SF->Is_Weak_FF[PV][PF];
	} else {
		// Curved FACE Integration
		OPS->NfnI    = ELEMENT->NfnIc[PF][IndFType];
		OPS->NfnI_SF = ELEMENT_SF->NfnIc[PF][0];
		OPS->NvnI_SF = ELEMENT_SF->NvnIc[PF];

		OPS->nOrdLR = ELEMENT_FACE->nOrd_fIc[PF][IndOrdLR];
		OPS->nOrdRL = ELEMENT_FACE->nOrd_fIc[PF][IndOrdRL];

		OPS->ChiS_fI    = ELEMENT->ChiS_fIc[PV][PF];
		OPS->ChiS_fI_sp = ELEMENT->ChiS_fIc_sp[PV][PF];
		OPS->ChiS_fI_SF = ELEMENT_SF->ChiS_fIc[PV][PF];
		OPS->ChiS_vI_SF = ELEMENT_SF->ChiS_vIc[PV][PF];

		OPS->I_Weak_FF = ELEMENT_SF->Ic_Weak_FF[PV][PF];
	}
}


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
	double       **D, *C, *Dxyz;

	Nbf  = DxyzInfo->Nbf;
	Nn   = DxyzInfo->Nn;
	dim1 = DxyzInfo->dim;

	D    = DxyzInfo->D;
	C    = DxyzInfo->C;

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
	struct S_OPERATORS *OPS;
	struct S_Dxyz      *DxyzInfo;

	OPS      = malloc(sizeof *OPS);      // free
	DxyzInfo = malloc(sizeof *DxyzInfo); // free

	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		unsigned int NvnS, NvnI;
		double       *ChiS_vI, **DxyzChiS;

		init_ops(OPS,VOLUME);

		NvnS = OPS->NvnS;
		NvnI = OPS->NvnI;

		DxyzInfo->Nbf = OPS->NvnS;
		DxyzInfo->Nn  = OPS->NvnI;
		DxyzInfo->D   = OPS->D_Weak;
		DxyzInfo->C   = VOLUME->C_vI;

		ChiS_vI = OPS->ChiS_vI[0];

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

static void compute_W_fI(const struct S_FDATA *FDATA, double *W_fI)
{
	/*
	 *	Purpose:
	 *		Interpolate VOLUME What coefficients to FACE cubature nodes.
	 *
	 *	Comments:
	 *		Various options are available for computing W_fI:
	 *			1) Using sum factorized operators for TP and WEDGE elements;
	 *			2) Using the standard operator but exploiting sparsity;
	 *			3) Using the standard approach.
	 */

	// Initialize DB Parameters
	unsigned int d            = DB.d,
	             Nvar         = DB.Nvar,
	             *VFPartUnity = DB.VFPartUnity,
	             ***SF_BE     = DB.SF_BE;

	// Standard datatypes
	unsigned int P, dim, Eclass, Vf, f, SpOp, IndFType, NfnI, NvnS,
	             NIn[DMAX], NOut[DMAX], Diag[DMAX], NOut0, NOut1;
	double       *OP[DMAX], **OPF0, **OPF1;

	struct S_OPERATORS **OPS;
	struct S_VOLUME    *VOLUME;

	OPS    = FDATA->OPS;
	VOLUME = FDATA->VOLUME;

	P        = FDATA->P;
	Eclass   = FDATA->Eclass;
	Vf       = FDATA->Vf;
	f        = FDATA->f;
	SpOp     = FDATA->SpOp;
	IndFType = FDATA->IndFType;

	NfnI = OPS[IndFType]->NfnI;
	NvnS = OPS[0]->NvnS;

	if (Eclass == C_TP && SF_BE[P][0][1]) {
		get_sf_parametersF(OPS[0]->NvnS_SF,OPS[0]->NvnI_SF,OPS[0]->ChiS_vI_SF,
		                   OPS[0]->NvnS_SF,OPS[0]->NfnI_SF,OPS[0]->ChiS_fI_SF,NIn,NOut,OP,d,Vf,C_TP);

		if (SpOp) {
			for (dim = 0; dim < d; dim++)
				Diag[dim] = 2;
			Diag[f/2] = 0;
		} else {
			for (dim = 0; dim < d; dim++)
				Diag[dim] = 0;
		}

		sf_apply_d(VOLUME->What,W_fI,NIn,NOut,Nvar,OP,Diag,d);
	} else if (Eclass == C_WEDGE && SF_BE[P][1][1]) {
		if (f < 3) { OPF0  = OPS[0]->ChiS_fI_SF, OPF1  = OPS[1]->ChiS_vI_SF;
		             NOut0 = OPS[0]->NfnI_SF,    NOut1 = OPS[1]->NvnI_SF;
		} else {     OPF0  = OPS[0]->ChiS_vI_SF, OPF1  = OPS[1]->ChiS_fI_SF;
		             NOut0 = OPS[0]->NvnI_SF,    NOut1 = OPS[1]->NfnI_SF; }
		get_sf_parametersF(OPS[0]->NvnS_SF,NOut0,OPF0,OPS[1]->NvnS_SF,NOut1,OPF1,NIn,NOut,OP,d,Vf,C_WEDGE);

		if (SpOp) {
			for (dim = 0; dim < d; dim++)
				Diag[dim] = 2;
			if (f < 3)
				Diag[0] = 0;
			else
				Diag[2] = 0;
		} else {
			for (dim = 0; dim < d; dim++)
				Diag[dim] = 0;
			Diag[1] = 2;
		}

		sf_apply_d(VOLUME->What,W_fI,NIn,NOut,Nvar,OP,Diag,d);
	} else if ((SpOp && (Eclass == C_TP || Eclass == C_WEDGE)) || (VFPartUnity[Eclass])) {
		mm_CTN_CSR_d(NfnI,Nvar,NvnS,OPS[0]->ChiS_fI_sp[Vf],VOLUME->What,W_fI);
	} else {
		mm_CTN_d(NfnI,Nvar,NvnS,OPS[0]->ChiS_fI[Vf],VOLUME->What,W_fI);
	}
}

static void init_FDATA(struct S_FDATA *FDATA, const struct S_FACE *FACE, const unsigned int side)
{
	// Initialize DB Parameters
	unsigned int Collocated = DB.Collocated;

	FDATA->P = FACE->P;

	if (side == 'L') {
		FDATA->VOLUME = FACE->VIn;
		FDATA->Vf     = FACE->VfIn;
	} else if (side == 'R') {
		FDATA->VOLUME = FACE->VOut;
		FDATA->Vf     = FACE->VfOut;
	} else {
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}
	FDATA->f    = (FDATA->Vf)/NFREFMAX;
	FDATA->SpOp = Collocated && ((FDATA->Vf) % NFREFMAX == 0 && FDATA->VOLUME->P == FDATA->P);

	FDATA->Eclass = FDATA->VOLUME->Eclass;
	FDATA->IndFType = get_IndFType(FDATA->Eclass,FDATA->f);

	init_opsF(FDATA->OPS[0],FDATA->VOLUME,FACE,0);
	if (FDATA->VOLUME->type == WEDGE || FDATA->VOLUME->type == PYR)
		// Needed for sum factorized operators and alternate FACE operators (TRIs/QUADs)
		init_opsF(FDATA->OPS[1],FDATA->VOLUME,FACE,1);
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
	struct S_OPERATORS *OPSL[2], *OPSR[2];
	struct S_FDATA     *FDATAL, *FDATAR;

	for (size_t i = 0; i < 2; i++) {
		OPSL[i]  = malloc(sizeof *OPSL[i]);  // free
		OPSR[i]  = malloc(sizeof *OPSR[i]);  // free
	}
	FDATAL = malloc(sizeof *FDATAL); // free
	FDATAR = malloc(sizeof *FDATAR); // free

	for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
		// Load required FACE operators (Left and Right FACEs)
		struct S_VOLUME *VL, *VR;

		FDATAL->OPS = OPSL;
		FDATAR->OPS = OPSR;

		init_FDATA(FDATAL,FACE,'L');
		init_FDATA(FDATAR,FACE,'R');

		VL = FDATAL->VOLUME;
		VR = FDATAR->VOLUME;

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
		compute_W_fI(FDATAL,WL_fIL);

		// Compute WR_fIL (Taking BCs into account if applicable)
		double *ChiSR_fIL = NULL;
		double *n_fI = FACE->n_fI;

		double *WR_fIL = malloc(NfnI*Nvar * sizeof *WR_fIL); // free
		if (!Boundary) {
			unsigned int *nOrdRL = OPSL[FDATAL->IndFType]->nOrdRL;

			// Use reordered operator (as opposed to solution) for ease of linearization.
			ChiSR_fIL = malloc(NfnI*NvnSR * sizeof *ChiSR_fIL); // free

			double *ChiS_fI = OPSR[0]->ChiS_fI[VfR];
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
			unsigned int *nOrdLR = OPSL[FDATAL->IndFType]->nOrdLR;
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
