// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "solver_Poisson_c.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "S_FACET.h"

#include "solver_Poisson.h"
#include "element_functions.h"
#include "matrix_functions.h"
#include "exact_solutions.h"
#include "array_swap.h"
#include "array_free.h"

/*
 *	Purpose:
 *		Compute RHS for Poisson solver using complex variables (for linearization testing).
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

struct S_OPERATORS {
	unsigned int NvnS, NvnI, NfnI, NvnC,
	             *nOrdOutIn, *nOrdInOut;
	double       *w_fI, *ChiS_vI, **D_Weak, **ChiS_fI, ***GradChiS_fI, **I_vC_fI;
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

		OPS->ChiS_vI   = ELEMENT->ChiS_vIs[P][P][0];
		OPS->D_Weak    = ELEMENT->Ds_Weak_VV[P][P][0];
	} else {
		OPS->NvnI = ELEMENT->NvnIc[P];

		OPS->ChiS_vI   = ELEMENT->ChiS_vIc[P][P][0];
		OPS->D_Weak    = ELEMENT->Dc_Weak_VV[P][P][0];
	}
}

static void init_opsF(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const struct S_FACET *FACET,
                      const unsigned int IndFType)
{
	// Standard datatypes
	unsigned int PV, PF, Vtype, Vcurved, FtypeInt, IndOrdOutIn, IndOrdInOut;

	struct S_ELEMENT *ELEMENT, *ELEMENT_FACET;

	PV      = VOLUME->P;
	PF      = FACET->P;
	Vtype   = VOLUME->type;
	Vcurved = VOLUME->curved;

	FtypeInt = FACET->typeInt;
	IndOrdOutIn = FACET->IndOrdOutIn;
	IndOrdInOut = FACET->IndOrdInOut;

	ELEMENT       = get_ELEMENT_type(Vtype);
	ELEMENT_FACET = get_ELEMENT_FACET(Vtype,IndFType);

	OPS->NvnS = ELEMENT->NvnS[PV];
	if (FtypeInt == 's') {
		// Straight FACET Integration
		OPS->NfnI = ELEMENT->NfnIs[PF][IndFType];

		OPS->w_fI = ELEMENT->w_fIs[PF][IndFType];

		OPS->ChiS_fI     = ELEMENT->ChiS_fIs[PV][PF];
		OPS->GradChiS_fI = ELEMENT->GradChiS_fIs[PV][PF];

		OPS->nOrdInOut = ELEMENT_FACET->nOrd_fIs[PF][IndOrdInOut];
		OPS->nOrdOutIn = ELEMENT_FACET->nOrd_fIs[PF][IndOrdOutIn];

		if (!Vcurved) {
			OPS->NvnC = ELEMENT->NvnCs[PV];
			OPS->I_vC_fI = ELEMENT->I_vCs_fIs[PV][PF];
		} else {
			OPS->NvnC = ELEMENT->NvnCc[PV];
			OPS->I_vC_fI = ELEMENT->I_vCc_fIs[PV][PF];
		}
	} else {
		// Curved FACET Integration
		OPS->NfnI = ELEMENT->NfnIc[PF][IndFType];

		OPS->w_fI = ELEMENT->w_fIc[PF][IndFType];

		OPS->ChiS_fI     = ELEMENT->ChiS_fIc[PV][PF];
		OPS->GradChiS_fI = ELEMENT->GradChiS_fIc[PV][PF];

		OPS->nOrdInOut = ELEMENT_FACET->nOrd_fIc[PF][IndOrdInOut];
		OPS->nOrdOutIn = ELEMENT_FACET->nOrd_fIc[PF][IndOrdOutIn];

		if (!Vcurved) {
			OPS->NvnC = ELEMENT->NvnCs[PV];
			OPS->I_vC_fI = ELEMENT->I_vCs_fIc[PV][PF];
		} else {
			OPS->NvnC = ELEMENT->NvnCc[PV];
			OPS->I_vC_fI = ELEMENT->I_vCc_fIc[PV][PF];
		}
	}
}

void compute_qhat_VOLUME_c(void)
{

	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int   dim1, NvnS;
	double         *MInv, **DxyzChiS, *Sxyz;
	double complex *qhat;

	struct S_OPERATORS *OPS;
	struct S_VOLUME    *VOLUME;

	OPS = malloc(sizeof *OPS); // free

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		init_ops(OPS,VOLUME);

		NvnS = OPS->NvnS;

		// Compute RHS term
		MInv     = VOLUME->MInv;
		DxyzChiS = VOLUME->DxyzChiS;

		for (dim1 = 0; dim1 < d; dim1++) {
			Sxyz = mm_Alloc_d(CBRM,CBNT,CBT,NvnS,NvnS,NvnS,1.0,MInv,DxyzChiS[dim1]); // free

			// RHS
			qhat = malloc(NvnS*1 * sizeof *qhat); // keep
			mm_dcc(CBCM,CBT,CBNT,NvnS,1,NvnS,1.0,0.0,Sxyz,VOLUME->uhat_c,qhat);
			free(Sxyz);

			if (VOLUME->qhat_c[dim1])
				free(VOLUME->qhat_c[dim1]);
			VOLUME->qhat_c[dim1] = qhat;
		}
	}
	free(OPS);
}

static void boundary_Poisson_c(const unsigned int Nn, const unsigned int Nel, double *XYZ, double complex *uL,
                               double complex *uR, double complex *graduL, double complex *graduR,
                               const unsigned int BC, const unsigned int curved)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int n;
	double       *XYZproj, *uB_d, *graduB_d;

	if (Nel != 1)
		printf("Error: Vectorization unsupported.\n"), EXIT_MSG;

	// Modify XYZ coordinates for boundary condition evaluation
	XYZproj = malloc(Nn*d * sizeof *XYZproj); // free
	if (d > 1)
		project_to_sphere(Nn,XYZ,XYZproj,curved);
	else
		printf("Error: Unsupported.\n"), EXIT_MSG;

	if (BC == BC_DIRICHLET) {
		uB_d = malloc(Nn * sizeof *uB_d); // free
		compute_exact_solution(Nn*Nel,XYZ,uB_d,0);
		for (n = 0; n < Nn; n++) {
			uR[n]  = 2.0*uB_d[n];
			uR[n] -= uL[n];
		}
		free(uB_d);

		if (graduL) {
			for (n = 0; n < d*Nn; n++)
				graduR[n] = graduL[n];
		}
	} else if (BC == BC_NEUMANN) {
		for (n = 0; n < Nn; n++)
			uR[n] = uL[n];

		if (graduL) {
			graduB_d = malloc(d*Nn * sizeof *graduB_d); // free
			compute_exact_gradient(Nn*Nel,XYZ,graduB_d);
			for (n = 0; n < Nn*d; n++) {
				graduR[n]  = 2.0*graduB_d[n];
				graduR[n] -= graduL[n];
			}
			free(graduB_d);
		}
	} else {
		printf("Error: Unsupported BC_type.\n"), EXIT_MSG;
	}
	free(XYZproj);
}

void compute_qhat_FACET_c(void)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int   n, dim,
	               NfnI, NvnSIn, NvnSOut,
	               VfIn, VfOut, fIn, EclassIn, IndFType, BC, Boundary,
	               *nOrdOutIn, *nOrdInOut;
	double         *w_fI, *wnJ_fI, *n_fI, *detJF_fI;
	double complex *uIn_fI, *uOut_fIIn, *uOut_fI, *uNum_fI, *qhatIn, *qhatOut;

	struct S_OPERATORS *OPSIn, *OPSOut;
	struct S_VOLUME    *VIn, *VOut;
	struct S_FACET     *FACET;

	// silence
	uOut_fI = NULL;

	OPSIn  = malloc(sizeof *OPSIn);  // free
	OPSOut = malloc(sizeof *OPSOut); // free

	for (FACET = DB.FACET; FACET; FACET = FACET->next) {
		VIn  = FACET->VIn;
		VfIn = FACET->VfIn;
		fIn  = VfIn/NFREFMAX;

		EclassIn = VIn->Eclass;
		IndFType = get_IndFType(EclassIn,fIn);
		init_opsF(OPSIn,VIn,FACET,IndFType);

		VOut  = FACET->VOut;
		VfOut = FACET->VfOut;

		init_opsF(OPSOut,VOut,FACET,IndFType);

		BC       = FACET->BC;
		Boundary = FACET->Boundary;

		NfnI    = OPSIn->NfnI;
		NvnSIn  = OPSIn->NvnS;
		NvnSOut = OPSOut->NvnS;

		w_fI = OPSIn->w_fI;

		nOrdOutIn = OPSIn->nOrdOutIn;
		nOrdInOut = OPSIn->nOrdInOut;

		detJF_fI = FACET->detJF_fI;
		n_fI     = FACET->n_fI;

		// Compute uIn_fI
		uIn_fI = malloc(NfnI * sizeof *uIn_fI); // free
		mm_dcc(CBCM,CBT,CBNT,NfnI,1,NvnSIn,1.0,0.0,OPSIn->ChiS_fI[VfIn],VIn->uhat_c,uIn_fI);

		// Compute_uOut_fI (Taking BCs into account if applicable)
		uOut_fIIn = malloc(NfnI * sizeof *uOut_fIIn); // free
		if (!Boundary) {
			uOut_fI = malloc(NfnI * sizeof *uOut_fI); // free
			mm_dcc(CBCM,CBT,CBNT,NfnI,1,NvnSOut,1.0,0.0,OPSOut->ChiS_fI[VfOut],VOut->uhat_c,uOut_fI);

			// Reorder uOut_fI to correspond to uIn_fI
			for (n = 0; n < NfnI; n++)
				uOut_fIIn[n] = uOut_fI[nOrdOutIn[n]];
		} else {
			boundary_Poisson_c(NfnI,1,FACET->XYZ_fI,uIn_fI,uOut_fIIn,NULL,NULL,BC % BC_STEP_SC, BC / BC_STEP_SC);
		}

		// Compute numerical trace
		uNum_fI = malloc(NfnI * sizeof *uNum_fI); // free

	 	// Numerical trace: uNum = 0.5 * (uL+uR)
		for (n = 0; n < NfnI; n++) {
			uNum_fI[n] = 0.5*(uIn_fI[n]+uOut_fIIn[n]);
		}
		free(uOut_fIIn);

		wnJ_fI = malloc(NfnI*d * sizeof *wnJ_fI); // free
		for (dim = 0; dim < d; dim++) {
		for (n = 0; n < NfnI; n++) {
			wnJ_fI[dim*NfnI+n] = w_fI[n]*detJF_fI[n]*n_fI[n*d+dim];
		}}

		// Compute partial FACET RHS terms

		// Interior VOLUME
		for (dim = 0; dim < d; dim++) {
			// RHSIn (partial)
			qhatIn = malloc(NfnI*1 * sizeof *qhatIn); // keep
			for (n = 0; n < NfnI; n++) {
				qhatIn[n] = wnJ_fI[dim*NfnI+n]*(uNum_fI[n]-uIn_fI[n]);
			}

			if (FACET->qhatIn_c[dim])
				free(FACET->qhatIn_c[dim]);
			FACET->qhatIn_c[dim] = qhatIn;
		}
		free(uIn_fI);

		// Exterior VOLUME
		if (!Boundary) {
			// Use "-ve" normal for opposite VOLUME
			for (n = 0; n < d*NfnI; n++) {
				wnJ_fI[n] *= -1.0;
			}

			// Rearrange numerical trace to match node ordering from opposite VOLUME
			array_rearrange_d(NfnI,d,nOrdInOut,'C',wnJ_fI);
			array_rearrange_cmplx(NfnI,1,nOrdInOut,'R',uNum_fI);

			for (dim = 0; dim < d; dim++) {
				// RHSOut (partial)
				qhatOut = malloc(NfnI*1 * sizeof *qhatOut); // keep
				for (n = 0; n < NfnI; n++) {
					qhatOut[n] = wnJ_fI[dim*NfnI+n]*(uNum_fI[n]-uOut_fI[n]);
				}

				if (FACET->qhatOut_c[dim])
					free(FACET->qhatOut_c[dim]);
				FACET->qhatOut_c[dim] = qhatOut;
			}
			free(uOut_fI);
		}
		free(uNum_fI);
		free(wnJ_fI);
	}
	free(OPSIn);
	free(OPSOut);
}

void compute_uhat_VOLUME_c(void)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int   dim1, NvnS;
	double         **DxyzChiS;
	double complex *RHS;

	struct S_OPERATORS *OPS;
	struct S_VOLUME    *VOLUME;

	OPS = malloc(sizeof *OPS); // free

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		init_ops(OPS,VOLUME);

		NvnS = OPS->NvnS;

		// Compute RHS terms
		DxyzChiS = VOLUME->DxyzChiS;

		// RHS
		if (VOLUME->RHS_c)
			free(VOLUME->RHS_c);
		RHS = calloc(NvnS , sizeof *RHS); // keep
		VOLUME->RHS_c = RHS;

		for (dim1 = 0; dim1 < d; dim1++)
			mm_dcc(CBCM,CBT,CBNT,NvnS,1,NvnS,-1.0,1.0,DxyzChiS[dim1],VOLUME->qhat_c[dim1],RHS);
	}
	free(OPS);
}

void compute_uhat_FACET_c()
{
	// Initialize DB Parameters
	unsigned int d               = DB.d,
	             ViscousFluxType = DB.ViscousFluxType;

	// Standard datatypes
	unsigned int   i, j, n, dim, dim1,
	               NvnSIn, NvnSOut, NfnI,
	               BC, Boundary, VfIn, VfOut, fIn, EclassIn, IndFType,
	               *nOrdOutIn, *nOrdInOut;
	double         **GradxyzIn, **GradxyzOut, *ChiS_fI, *ChiS_fI_std, *ChiSMInv_fIIn, *ChiSMInv_fIOut, *I_jump,
	               *gradu_avg, *u_jump, *u_jump_part,
	               *detJVIn_fI, *detJVOut_fI, *h, *n_fI, *detJF_fI, *w_fI, VolL, VolR;
	double complex *uIn_fI, *grad_uIn_fI, *uOut_fIIn, *grad_uOut_fIIn, *uOut_fI,
	               *nqNum_fI, *RHSIn, *RHSOut;

	struct S_OPERATORS *OPSIn, *OPSOut;
	struct S_FACET     *FACET;
	struct S_VOLUME    *VIn, *VOut;

	OPSIn  = malloc(sizeof *OPSIn);  // free
	OPSOut = malloc(sizeof *OPSOut); // free

	for (FACET = DB.FACET; FACET; FACET = FACET->next) {
		VIn  = FACET->VIn;
		VfIn = FACET->VfIn;
		fIn  = VfIn/NFREFMAX;

		EclassIn = VIn->Eclass;
		IndFType = get_IndFType(EclassIn,fIn);
		init_opsF(OPSIn,VIn,FACET,IndFType);

		VOut  = FACET->VOut;
		VfOut = FACET->VfOut;

		init_opsF(OPSOut,VOut,FACET,IndFType);

		BC       = FACET->BC;
		Boundary = FACET->Boundary;

		NfnI    = OPSIn->NfnI;
		NvnSIn  = OPSIn->NvnS;
		NvnSOut = OPSOut->NvnS;

		w_fI = OPSIn->w_fI;

		nOrdOutIn = OPSIn->nOrdOutIn;
		nOrdInOut = OPSIn->nOrdInOut;

		n_fI = FACET->n_fI;

		detJVIn_fI = FACET->detJVIn_fI;
		if (!Boundary) {
			detJVOut_fI = malloc(NfnI * sizeof *detJVOut_fI); // free
			for (n = 0; n < NfnI; n++)
				detJVOut_fI[n] = FACET->detJVOut_fI[nOrdOutIn[n]];
		} else {
			detJVOut_fI = detJVIn_fI;
		}

		detJF_fI = FACET->detJF_fI;
		h = malloc(NfnI * sizeof *h); // free
		for (n = 0; n < NfnI; n++)
			h[n] = max(detJVIn_fI[n],detJVOut_fI[n])/detJF_fI[n];

		// Compute uIn_fI and gradu_In_fI
		uIn_fI      = malloc(NfnI   * sizeof *uIn_fI);      // free
		grad_uIn_fI = calloc(NfnI*d , sizeof *grad_uIn_fI); // free

		ChiSMInv_fIIn  = malloc(NfnI*NvnSIn  * sizeof *ChiSMInv_fIIn);  // free
		ChiSMInv_fIOut = malloc(NfnI*NvnSOut * sizeof *ChiSMInv_fIOut); // free

		// Note: GradxyzOut only needed if not on a boundary (ToBeDeleted)
		GradxyzIn  = malloc(d * sizeof *GradxyzIn);  // free
		GradxyzOut = malloc(d * sizeof *GradxyzOut); // free

		// Note: This is analogous to Hesthaven's BuildCurvedOPS2D (Chapter 9.1)
		// See solver_Poisson_weak under 'alternate/' for alternative computation of GradxyzIn/Out
		mm_d(CBRM,CBNT,CBNT,NfnI,NvnSIn, NvnSIn, 1.0,0.0,OPSIn->ChiS_fI[VfIn],  VIn->MInv, ChiSMInv_fIIn);
		mm_d(CBRM,CBNT,CBNT,NfnI,NvnSOut,NvnSOut,1.0,0.0,OPSOut->ChiS_fI[VfOut],VOut->MInv,ChiSMInv_fIOut);
		array_rearrange_d(NfnI,NvnSOut,nOrdOutIn,'R',ChiSMInv_fIOut);

		// Note: GradxyzOut is arranged according to the FACET node ordering of VIn.
		for (dim1 = 0; dim1 < d; dim1++) {
			GradxyzIn[dim1]  = malloc(NfnI*NvnSIn  * sizeof **GradxyzIn);  // free
			GradxyzOut[dim1] = malloc(NfnI*NvnSOut * sizeof **GradxyzOut); // free

			mm_d(CBRM,CBNT,CBT,NfnI,NvnSIn, NvnSIn, 1.0,0.0,ChiSMInv_fIIn, VIn->DxyzChiS[dim1], GradxyzIn[dim1]);
			mm_d(CBRM,CBNT,CBT,NfnI,NvnSOut,NvnSOut,1.0,0.0,ChiSMInv_fIOut,VOut->DxyzChiS[dim1],GradxyzOut[dim1]);
		}

		mm_dcc(CBCM,CBT,CBNT,NfnI,1,NvnSIn,1.0,0.0,OPSIn->ChiS_fI[VfIn],VIn->uhat_c,uIn_fI);

		for (dim = 0; dim < d; dim++)
			mm_dcc(CBCM,CBT,CBNT,NfnI,1,NvnSIn,1.0,0.0,GradxyzIn[dim],VIn->uhat_c,&grad_uIn_fI[NfnI*dim]);

		// Compute_uOut_fI (Taking BCs into account if applicable)
		uOut_fIIn      = malloc(NfnI   * sizeof *uOut_fIIn);      // free
		grad_uOut_fIIn = calloc(NfnI*d , sizeof *grad_uOut_fIIn); // free

		if (!Boundary) {
			uOut_fI = malloc(NfnI * sizeof *uOut_fI); // free

			mm_dcc(CBCM,CBT,CBNT,NfnI,1,NvnSOut,1.0,0.0,OPSOut->ChiS_fI[VfOut],VOut->uhat_c,uOut_fI);

			// Reorder uOut_fI to correspond to inner VOLUME ordering
			for (n = 0; n < NfnI; n++)
				uOut_fIIn[n] = uOut_fI[nOrdOutIn[n]];
			free(uOut_fI);

			for (dim = 0; dim < d; dim++)
				mm_dcc(CBCM,CBT,CBNT,NfnI,1,NvnSOut,1.0,0.0,GradxyzOut[dim],VOut->uhat_c,&grad_uOut_fIIn[NfnI*dim]);
		} else {
			boundary_Poisson_c(NfnI,1,FACET->XYZ_fI,uIn_fI,uOut_fIIn,grad_uIn_fI,grad_uOut_fIIn,BC % BC_STEP_SC,BC / BC_STEP_SC);
		}

		// Compute numerical flux
		nqNum_fI = calloc(NfnI , sizeof *nqNum_fI); // free

		gradu_avg = malloc(NfnI*d * sizeof *gradu_avg); // free
		u_jump    = malloc(NfnI*d * sizeof *u_jump);    // free

		jacobian_flux_coef(NfnI,1,n_fI,h,FACET->P,gradu_avg,u_jump,d,ViscousFluxType);
		free(h);

		switch(ViscousFluxType) {
		case FLUX_IP:
			// No modifications needed.
			break;
		case FLUX_BR2:
			// Assemble BR2 scaling term for u_jump
			I_jump = malloc(NfnI*NfnI * sizeof *I_jump); // free

			if (!Boundary) {
				ChiS_fI = OPSIn->ChiS_fI[VfIn];
				mm_d(CBRM,CBNT,CBT,NfnI,NfnI,NvnSIn,1.0,0.0,ChiSMInv_fIIn,ChiS_fI,I_jump);

				ChiS_fI_std = OPSOut->ChiS_fI[VfOut];
				ChiS_fI = malloc(NfnI*NvnSOut * sizeof *ChiS_fI); // free

				// Reordering
				for (i = 0; i < NfnI; i++) {
				for (j = 0; j < NvnSOut; j++) {
					ChiS_fI[i*NvnSOut+j] = ChiS_fI_std[nOrdInOut[i]*NvnSOut+j];
				}}
				mm_d(CBRM,CBNT,CBT,NfnI,NfnI,NvnSOut,1.0,1.0,ChiSMInv_fIOut,ChiS_fI,I_jump);
				free(ChiS_fI);
			} else {
				ChiS_fI = OPSIn->ChiS_fI[VfIn];
				mm_d(CBRM,CBNT,CBT,NfnI,NfnI,NvnSIn,2.0,0.0,ChiSMInv_fIIn,ChiS_fI,I_jump);
			}

			u_jump_part = malloc(NfnI*d * sizeof *u_jump_part); // free
			for (dim = 0; dim < d; dim++) {
			for (n = 0; n < NfnI; n++) {
				u_jump_part[NfnI*dim+n] = u_jump[NfnI*dim+n]*detJF_fI[n]*w_fI[n];
			}}
			mm_CTN_d(NfnI,d,NfnI,I_jump,u_jump_part,u_jump);

			free(I_jump);
			free(u_jump_part);
			break;
		case FLUX_CDG2:
			// Note: The simplification of the general term in Brdar(2012) Table 1 is from eq. (4.3) defining Beta based
			//       on the area switch.

			// Assemble 2*BR2 scaling term from appropriate side.
			I_jump = malloc(NfnI*NfnI * sizeof *I_jump); // free

			if (Boundary) {
				ChiS_fI = OPSIn->ChiS_fI[VfIn];
				mm_d(CBRM,CBNT,CBT,NfnI,NfnI,NvnSIn,2.0,0.0,ChiSMInv_fIIn,ChiS_fI,I_jump);
			} else {
				// Evaluate from which side scaling should be computed based on area switch (Brdar(2012))
				VolL = 0.0;
				VolR = 0.0;
				for (n = 0; n < NfnI; n++) {
					VolL = max(VolL,detJVIn_fI[n]);
					VolR = max(VolR,detJVOut_fI[n]);
				}

				if (VolL <= VolR) {
					ChiS_fI = OPSIn->ChiS_fI[VfIn];
					mm_d(CBRM,CBNT,CBT,NfnI,NfnI,NvnSIn,2.0,0.0,ChiSMInv_fIIn,ChiS_fI,I_jump);
				} else {
					ChiS_fI_std = OPSOut->ChiS_fI[VfOut];
					ChiS_fI = malloc(NfnI*NvnSOut * sizeof *ChiS_fI); // free

					// Reordering
					for (i = 0; i < NfnI; i++) {
					for (j = 0; j < NvnSOut; j++) {
						ChiS_fI[i*NvnSOut+j] = ChiS_fI_std[nOrdInOut[i]*NvnSOut+j];
					}}
					mm_d(CBRM,CBNT,CBT,NfnI,NfnI,NvnSOut,2.0,0.0,ChiSMInv_fIOut,ChiS_fI,I_jump);
					free(ChiS_fI);
				}
			}

			u_jump_part = malloc(NfnI*d * sizeof *u_jump_part); // free
			for (dim = 0; dim < d; dim++) {
			for (n = 0; n < NfnI; n++) {
				u_jump_part[NfnI*dim+n] = u_jump[NfnI*dim+n]*detJF_fI[n]*w_fI[n];
			}}
			mm_CTN_d(NfnI,d,NfnI,I_jump,u_jump_part,u_jump);

			free(I_jump);
			free(u_jump_part);
			break;
		default:
			printf("Error: Unsupported ViscousFluxType.\n"), EXIT_MSG;
			break;
		}
		if (!Boundary)
			free(detJVOut_fI);

		free(ChiSMInv_fIIn);
		free(ChiSMInv_fIOut);

		for (dim = 0; dim < d; dim++) {
		for (n = 0; n < NfnI; n++) {
			nqNum_fI[n] += n_fI[n*d+dim]*
				(  gradu_avg[NfnI*dim+n] * (grad_uIn_fI[NfnI*dim+n] + grad_uOut_fIIn[NfnI*dim+n])
				 + u_jump[NfnI*dim+n]    * (uIn_fI[n] - uOut_fIIn[n]));
		}}
		free(gradu_avg);
		free(u_jump);

		free(grad_uIn_fI);
		free(grad_uOut_fIIn);
		free(uIn_fI);
		free(uOut_fIIn);

		// Multiply by area element and cubature weights
		for (n = 0; n < NfnI; n++) {
			nqNum_fI[n] *= detJF_fI[n]*w_fI[n];
		}

		// Finalize FACET RHS terms

		// Add VOLUME contributions to RHS
		RHSIn  = calloc(NvnSIn  , sizeof *RHSIn);  // keep (requires external free)
		RHSOut = calloc(NvnSOut , sizeof *RHSOut); // keep (requires external free)

		if (FACET->RHSIn_c)
			free(FACET->RHSIn_c);
		FACET->RHSIn_c  = RHSIn;
		if (FACET->RHSOut_c)
			free(FACET->RHSOut_c);
		FACET->RHSOut_c = RHSOut;

		for (dim = 0; dim < d; dim++) {
			// RHS
			mm_dcc(CBCM,CBNT,CBNT,NvnSIn, 1,NfnI,-1.0,1.0,GradxyzIn[dim], FACET->qhatIn_c[dim], RHSIn);
			if (!Boundary) {
				array_rearrange_d(NfnI,NvnSOut,nOrdInOut,'R',GradxyzOut[dim]);
				mm_dcc(CBCM,CBNT,CBNT,NvnSOut,1,NfnI,-1.0,1.0,GradxyzOut[dim],FACET->qhatOut_c[dim],RHSOut);
			}
		}
		array_free2_d(d,GradxyzIn);
		array_free2_d(d,GradxyzOut);

		// Add FACET contributions to RHS

		// Interior FACET
		mm_dcc(CBCM,CBNT,CBNT,NvnSIn,1,NfnI,1.0,1.0,OPSIn->ChiS_fI[VfIn],nqNum_fI,RHSIn);

		// Exterior FACET
		if (!Boundary) {
			// RHS

			// Use -ve normal for opposite VOLUME
			for (n = 0; n < NfnI; n++)
				nqNum_fI[n] *= -1.0;

			// Rearrange nqNum to match node ordering from VOut
			array_rearrange_cmplx(NfnI,1,nOrdInOut,'R',nqNum_fI);

			mm_dcc(CBCM,CBNT,CBNT,NvnSOut,1,NfnI,1.0,1.0,OPSOut->ChiS_fI[VfOut],nqNum_fI,RHSOut);
		}
		free(nqNum_fI);
	}
	free(OPSIn);
	free(OPSOut);
}
