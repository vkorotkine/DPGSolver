// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "solver_poisson_c.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "S_FACET.h"

#include "solver_poisson.h"
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
	double       *ChiS_vI, **D_Weak, **ChiS_fI, **I_Weak_FF, ***GradChiS_fI, **I_vC_fI;
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

		OPS->ChiS_fI     = ELEMENT->ChiS_fIs[PV][PF];
		OPS->I_Weak_FF   = ELEMENT->Is_Weak_FF[PV][PF];
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

		OPS->ChiS_fI     = ELEMENT->ChiS_fIc[PV][PF];
		OPS->I_Weak_FF   = ELEMENT->Ic_Weak_FF[PV][PF];
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
			Sxyz = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS,NvnS,NvnS,-1.0,MInv,DxyzChiS[dim1]); // free

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

static void boundary_Dirichlet_c(const unsigned int Nn, const unsigned int Nel, double *XYZ, double complex *uL,
                                 double complex *uR)
{
	unsigned int n;
	double       *uB_d;

	if (Nel != 1)
		printf("Error: Vectorization unsupported.\n"), EXIT_MSG;

	uB_d = malloc(Nn * sizeof *uB_d); // free
	compute_exact_solution(Nn*Nel,XYZ,uB_d,0);

	for (n = 0; n < Nn; n++) {
		uR[n] = 2.0*uB_d[n];
		uR[n] -= uL[n];
	}
	free(uB_d);
}

void compute_qhat_FACET_c(void)
{
	// Initialize DB Parameters
	unsigned int d               = DB.d,
	             ViscousFluxType = DB.ViscousFluxType;

	// Standard datatypes
	unsigned int   n, dim, i, iMax,
	               NfnI, NvnSIn, NvnSOut,
	               VfIn, VfOut, fIn, EclassIn, IndFType, BC, Boundary,
	               *nOrdOutIn, *nOrdInOut;
	double         nJ, *n_fI, *detJF_fI, *u_avg, *u_jump,
	               *I_FF, *MInvI_FF;
	double complex *uIn_fI, *uOut_fIIn, *uOut_fI, *uNum_fI, *qhatIn, *qhatOut, *nuNum_fI;

	struct S_OPERATORS *OPSIn, *OPSOut;
	struct S_VOLUME    *VIn, *VOut;
	struct S_FACET     *FACET;

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
			free(uOut_fI);
		} else {
			if (BC % BC_STEP_SC == BC_DIRICHLET) {
				boundary_Dirichlet_c(NfnI,1,FACET->XYZ_fI,uIn_fI,uOut_fIIn);
			} else if (BC % BC_STEP_SC == BC_NEUMANN) {
				printf("Add support.\n"), EXIT_MSG;
			} else {
				printf("Error: Unsupported BC.\n"), EXIT_MSG;
			}
		}

		// Compute numerical trace
		u_avg   = malloc(NfnI * sizeof *u_avg);   // free
		u_jump  = malloc(NfnI * sizeof *u_jump);  // free
		uNum_fI = malloc(NfnI * sizeof *uNum_fI); // free
		switch (ViscousFluxType) {
		case FLUX_IP:
			trace_coef(NfnI,1,n_fI,u_avg,u_jump,d,"IP");
			break;
		default:
			printf("Error: Unsupported ViscousFluxType.\n"), EXIT_MSG;
			break;
		}

		for (n = 0; n < NfnI; n++)
			uNum_fI[n] = u_avg[n]*(uIn_fI[n]+uOut_fIIn[n]) + u_jump[n]*(uIn_fI[n]-uOut_fIIn[n]);

		free(u_avg);
		free(u_jump);

		free(uIn_fI);
		free(uOut_fIIn);

		// Multiply uNum by the normal vector and area element
		nuNum_fI = malloc(NfnI*d * sizeof *nuNum_fI); // free
		for (dim = 0; dim < d; dim++) {
		for (n = 0; n < NfnI; n++) {
			nJ = n_fI[n*d+dim]*detJF_fI[n];
			nuNum_fI[dim*NfnI+n] = uNum_fI[n]*nJ;
		}}
		free(uNum_fI);

		// Compute FACET RHS terms

		// Interior VOLUME

		I_FF     = OPSIn->I_Weak_FF[VfIn];
		MInvI_FF = mm_Alloc_d(CBRM,CBNT,CBNT,NvnSIn,NfnI,NvnSIn,-1.0,VIn->MInv,I_FF); // free

		for (dim = 0; dim < d; dim++) {
			// RHSIn
			qhatIn = malloc(NvnSIn*1 * sizeof *qhatIn); // keep
			mm_dcc(CBCM,CBT,CBNT,NvnSIn,1,NfnI,1.0,0.0,MInvI_FF,&nuNum_fI[NfnI*dim],qhatIn);

			if (FACET->qhatIn_c[dim])
				free(FACET->qhatIn_c[dim]);
			FACET->qhatIn_c[dim] = qhatIn;
		}

		// Exterior VOLUME
		if (!Boundary) {
			free(MInvI_FF);

			// Use "-ve" normal for opposite VOLUME
			for (i = 0, iMax = NfnI*d; i < iMax; i++) {
				nuNum_fI[i] *= -1.0;
			}

			// Rearrange numerical trace to match node ordering from opposite VOLUME
			array_rearrange_cmplx(NfnI,d,nOrdInOut,'C',nuNum_fI);

			I_FF     = OPSOut->I_Weak_FF[VfOut];
			MInvI_FF = mm_Alloc_d(CBRM,CBNT,CBNT,NvnSOut,NfnI,NvnSOut,-1.0,VOut->MInv,I_FF); // free

			for (dim = 0; dim < d; dim++) {
				// RHSOut
				qhatOut = malloc(NvnSOut*1 * sizeof *qhatOut); // keep
				mm_dcc(CBCM,CBT,CBNT,NvnSOut,1,NfnI,1.0,0.0,MInvI_FF,&nuNum_fI[NfnI*dim],qhatOut);

				if (FACET->qhatOut_c[dim])
					free(FACET->qhatOut_c[dim]);
				FACET->qhatOut_c[dim] = qhatOut;
			}
		}
		free(MInvI_FF);

		free(nuNum_fI);
	}
	free(OPSIn);
	free(OPSOut);
}

void finalize_qhat_c(void)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int   dim, iMax, NvnSIn, NvnSOut;
	double complex *VqhatIn_ptr, *FqhatIn_ptr, *VqhatOut_ptr, *FqhatOut_ptr;

	struct S_FACET  *FACET;
	struct S_VOLUME *VIn, *VOut;

	for (FACET = DB.FACET; FACET; FACET = FACET->next) {
		VIn    = FACET->VIn;
		NvnSIn = VIn->NvnS;

		for (dim = 0; dim < d; dim++) {
			VqhatIn_ptr  = VIn->qhat_c[dim];
			FqhatIn_ptr  = FACET->qhatIn_c[dim];

			for (iMax = NvnSIn; iMax--; )
				*VqhatIn_ptr++ += *FqhatIn_ptr;
		}

		if(!(FACET->Boundary)) {
			VOut    = FACET->VOut;
			NvnSOut = VOut->NvnS;

			for (dim = 0; dim < d; dim++) {
				VqhatOut_ptr = VOut->qhat_c[dim];
				FqhatOut_ptr = FACET->qhatOut_c[dim];

				for (iMax = NvnSOut; iMax--; )
					*VqhatOut_ptr++ += *FqhatOut_ptr;
			}
		}
	}
}

void compute_uhat_VOLUME_c(void)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int   dim1, NvnI, NvnS;
	double         *ChiS_vI, **DxyzChiS;
	double complex *q_vI, *RHS;

	struct S_OPERATORS *OPS;
	struct S_VOLUME    *VOLUME;

	OPS = malloc(sizeof *OPS); // free

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		init_ops(OPS,VOLUME);

		NvnI = OPS->NvnI;
		NvnS = OPS->NvnS;

		ChiS_vI = OPS->ChiS_vI;

		// Obtain q_vI
		q_vI = malloc(NvnI*d * sizeof *q_vI); // free
		for (dim1 = 0; dim1 < d; dim1++)
			mm_dcc(CBCM,CBT,CBNT,NvnI,1,NvnS,1.0,0.0,ChiS_vI,VOLUME->qhat_c[dim1],&q_vI[dim1*NvnI]);

		// Compute RHS terms
		DxyzChiS = VOLUME->DxyzChiS;

		// RHS
		if (VOLUME->RHS_c)
			free(VOLUME->RHS_c);
		RHS = calloc(NvnS , sizeof *RHS); // keep
		VOLUME->RHS_c = RHS;

		for (dim1 = 0; dim1 < d; dim1++)
			mm_dcc(CBCM,CBT,CBNT,NvnS,1,NvnS,-1.0,1.0,DxyzChiS[dim1],VOLUME->qhat_c[dim1],RHS);

		free(q_vI);
	}
	free(OPS);
}

void compute_uhat_FACET_c()
{
	// Initialize DB Parameters
	unsigned int d               = DB.d,
	             ViscousFluxType = DB.ViscousFluxType;

	// Standard datatypes
	unsigned int   j, n, dim, dim1, dim2, IndC,
	               NvnSIn, NvnSOut, NfnI, NvnCIn,
	               BC, Boundary, VfIn, VfOut, fIn, EclassIn, IndFType,
	               *nOrdOutIn, *nOrdInOut;
	double         ***GradChiS_fI, **GradxyzIn, **GradxyzOut,
	               *gradu_avg, *u_jump,
	               *detJVIn_fI, *detJVOut_fI, *C_fI, *h, *n_fI, *detJF_fI, *C_vC;
	double complex *uIn_fI, *grad_uIn_fI, *uOut_fIIn, *grad_uOut_fIIn, *uOut_fI, **qhatIn_fI, **qhatOut_fIIn,
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

		NfnI   = OPSIn->NfnI;
		NvnSIn  = OPSIn->NvnS;
		NvnSOut = OPSOut->NvnS;
		NvnCIn  = OPSIn->NvnC;

		nOrdOutIn = OPSIn->nOrdOutIn;
		nOrdInOut = OPSIn->nOrdInOut;

		C_vC = VIn->C_vC;
		C_fI = malloc(NfnI*d*d * sizeof *C_fI); // free
		mm_CTN_d(NfnI,d*d,NvnCIn,OPSIn->I_vC_fI[VfIn],C_vC,C_fI);
		n_fI = FACET->n_fI;

		detJVIn_fI = FACET->detJVIn_fI;
		if (!Boundary) {
			detJVOut_fI = FACET->detJVOut_fI;

			// Reorder detJVOut_fI
			array_rearrange_d(NfnI,1,nOrdOutIn,'R',detJVOut_fI);
		} else {
			detJVOut_fI = detJVIn_fI;
		}
		
		detJF_fI = FACET->detJF_fI;
		h = malloc(NfnI * sizeof *h); // free
		for (n = 0; n < NfnI; n++)
			h[n] = max(detJVIn_fI[n],detJVOut_fI[n])/detJF_fI[n];

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
			mm_dcc(CBCM,CBT,CBNT,NvnSIn, 1,NvnSIn, -1.0,1.0,VIn->DxyzChiS[dim], FACET->qhatIn_c[dim], RHSIn);
			if (!Boundary)
				mm_dcc(CBCM,CBT,CBNT,NvnSOut,1,NvnSOut,-1.0,1.0,VOut->DxyzChiS[dim],FACET->qhatOut_c[dim],RHSOut);
		}

		// Compute uIn_fI and gradu_In_fI

		uIn_fI      = malloc(NfnI   * sizeof *uIn_fI);      // free
		grad_uIn_fI = calloc(NfnI*d , sizeof *grad_uIn_fI); // free

		mm_dcc(CBCM,CBT,CBNT,NfnI,1,NvnSIn,1.0,0.0,OPSIn->ChiS_fI[VfIn],VIn->uhat_c,uIn_fI);

		GradChiS_fI = OPSIn->GradChiS_fI;

		GradxyzIn = malloc(d * sizeof *GradxyzIn); // free
		for (dim1 = 0; dim1 < d; dim1++) {
			GradxyzIn[dim1] = calloc(NfnI*NvnSIn , sizeof **GradxyzIn); // free

			for (dim2 = 0; dim2 < d; dim2++) {
				IndC = (dim1+dim2*d)*NfnI;
				for (n = 0; n < NfnI; n++) {
					for (j = 0; j < NvnSIn; j++)
						GradxyzIn[dim1][n*NvnSIn+j] += GradChiS_fI[VfIn][dim2][n*NvnSIn+j]*C_fI[IndC+n];
				}
			}
			for (n = 0; n < NfnI; n++) {
				for (j = 0; j < NvnSIn; j++)
					GradxyzIn[dim1][n*NvnSIn+j] /= detJVIn_fI[n];
			}
		}

		for (dim = 0; dim < d; dim++)
			mm_dcc(CBCM,CBT,CBNT,NfnI,1,NvnSIn,1.0,0.0,GradxyzIn[dim],VIn->uhat_c,&grad_uIn_fI[NfnI*dim]);
		array_free2_d(d,GradxyzIn);

		mm_CTN_d(NfnI,d*d,OPSOut->NvnC,OPSOut->I_vC_fI[VfOut],VOut->C_vC,C_fI);
		GradChiS_fI = OPSOut->GradChiS_fI;

		// Note: Rearrangement is embedded in the operator
		GradxyzOut = malloc(d * sizeof *GradxyzOut); // free
		for (dim1 = 0; dim1 < d; dim1++) {
			GradxyzOut[dim1] = calloc(NfnI*NvnSIn , sizeof **GradxyzOut); // free

			for (dim2 = 0; dim2 < d; dim2++) {
				IndC = (dim1+dim2*d)*NfnI;
				for (n = 0; n < NfnI; n++) {
					for (j = 0; j < NvnSOut; j++)
						GradxyzOut[dim1][n*NvnSOut+j] += GradChiS_fI[VfOut][dim2][nOrdOutIn[n]*NvnSOut+j]*C_fI[IndC+n];
				}
			}
			for (n = 0; n < NfnI; n++) {
				for (j = 0; j < NvnSOut; j++)
					GradxyzOut[dim1][n*NvnSOut+j] /= detJVOut_fI[n];
			}
		}
		free(C_fI);

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
			if (BC % BC_STEP_SC == BC_DIRICHLET) {
				boundary_Dirichlet_c(NfnI,1,FACET->XYZ_fI,uIn_fI,uOut_fIIn);
				for (dim = 0; dim < d; dim++) {
				for (n = 0; n < NfnI; n++) {
					grad_uOut_fIIn[dim*NfnI+n] = grad_uIn_fI[dim*NfnI+n];
				}}
			} else if (BC % BC_STEP_SC == BC_NEUMANN) {
				printf("Add support.\n"), EXIT_MSG;
			} else {
				printf("Error: Unsupported BC.\n"), EXIT_MSG;
			}
		}
		array_free2_d(d,GradxyzOut);

		// Compute numerical flux
		nqNum_fI = calloc(NfnI , sizeof *nqNum_fI); // free

		gradu_avg = malloc(NfnI*d * sizeof *gradu_avg); // free
		u_jump    = malloc(NfnI*d * sizeof *u_jump);    // free

		switch (ViscousFluxType) {
		case FLUX_IP:
			jacobian_flux_coef(NfnI,1,n_fI,h,FACET->P,gradu_avg,u_jump,d,"IP",'L');
			break;
		default:
			printf("Error: Unsupported PoissonFluxType.\n"), EXIT_MSG;
			break;
		}
		free(h);

		qhatIn_fI    = malloc(d * sizeof *qhatIn_fI);    // free
		qhatOut_fIIn = malloc(d * sizeof *qhatOut_fIIn); // free

		for (dim = 0; dim < d; dim++) {
			qhatIn_fI[dim]    = malloc(NfnI * sizeof *qhatIn_fI[dim]);    // free
			qhatOut_fIIn[dim] = malloc(NfnI * sizeof *qhatOut_fIIn[dim]); // free

			mm_dcc(CBCM,CBT,CBNT,NfnI,1,NvnSIn, 1.0,0.0,OPSIn->ChiS_fI[VfIn],  VIn->qhat_c[dim], qhatIn_fI[dim]);
			mm_dcc(CBCM,CBT,CBNT,NfnI,1,NvnSOut,1.0,0.0,OPSOut->ChiS_fI[VfOut],VOut->qhat_c[dim],qhatOut_fIIn[dim]);
			array_rearrange_cmplx(NfnI,1,nOrdOutIn,'C',qhatOut_fIIn[dim]);
		}

		for (dim = 0; dim < d; dim++) {
		for (n = 0; n < NfnI; n++) {
			nqNum_fI[n] += gradu_avg[NfnI*dim+n] * (grad_uIn_fI[NfnI*dim+n] + grad_uOut_fIIn[NfnI*dim+n])
			            +  u_jump[NfnI*dim+n]    * (uIn_fI[n] - uOut_fIIn[n]);
		}}
		free(gradu_avg);
		free(u_jump);

		free(grad_uIn_fI);
		free(grad_uOut_fIIn);
		free(uIn_fI);
		free(uOut_fIIn);

		array_free2_cmplx(d,qhatIn_fI);
		array_free2_cmplx(d,qhatOut_fIIn);

		// Multiply by area element
		for (n = 0; n < NfnI; n++) {
			nqNum_fI[n] *= detJF_fI[n];
		}

		// Finalize FACET RHS terms

		// Interior FACET
		mm_dcc(CBCM,CBT,CBNT,NvnSIn,1,NfnI,-1.0,1.0,OPSIn->I_Weak_FF[VfIn],nqNum_fI,RHSIn);

		// Exterior FACET
		if (!Boundary) {
			// RHS

			// Use -ve normal for opposite VOLUME
			for (n = 0; n < NfnI; n++)
				nqNum_fI[n] *= -1.0;

			// Rearrange nqNum to match node ordering from VOut
			array_rearrange_cmplx(NfnI,1,nOrdInOut,'C',nqNum_fI);

			mm_dcc(CBCM,CBT,CBNT,NvnSOut,1,NfnI,-1.0,1.0,OPSOut->I_Weak_FF[VfOut],nqNum_fI,RHSOut);
		}
		free(nqNum_fI);
	}
	free(OPSIn);
	free(OPSOut);
}
