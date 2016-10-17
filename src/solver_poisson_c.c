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

#include "array_print.h"

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
	double       *w_fI, *ChiS_vI, **D_Weak, **ChiS_fI, **I_Weak_FF, ***GradChiS_fI, **I_vC_fI;
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

		OPS->w_fI = ELEMENT->w_fIc[PF][IndFType];

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
			if (BC % BC_STEP_SC == BC_DIRICHLET) {
				boundary_Dirichlet_c(NfnI,1,FACET->XYZ_fI,uIn_fI,uOut_fIIn);
			} else if (BC % BC_STEP_SC == BC_NEUMANN) {
				printf("Add support.\n"), EXIT_MSG;
			} else {
				printf("Error: Unsupported BC.\n"), EXIT_MSG;
			}
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
	unsigned int   n, dim, dim1,
	               NvnSIn, NvnSOut, NfnI, NvnCIn,
	               BC, Boundary, VfIn, VfOut, fIn, EclassIn, IndFType,
	               *nOrdOutIn, *nOrdInOut;
	double         **GradxyzIn, **GradxyzOut, *SxyzIn, *SxyzOut,
	               *gradu_avg, *u_jump,
	               *detJVIn_fI, *detJVOut_fI, *C_fI, *h, *n_fI, *detJF_fI, *C_vC;
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

		SxyzIn  = malloc(NvnSIn*NvnSIn   * sizeof *SxyzIn);  // free
		SxyzOut = malloc(NvnSOut*NvnSOut * sizeof *SxyzOut); // free

		// Note: GradxyzOut only needed if not on a boundary (ToBeDeleted)
		GradxyzIn  = malloc(d * sizeof *GradxyzIn);  // free
		GradxyzOut = malloc(d * sizeof *GradxyzOut); // free
		for (dim1 = 0; dim1 < d; dim1++) {
			GradxyzIn[dim1] = malloc(NfnI*NvnSIn * sizeof **GradxyzIn); // free

			mm_d(CBRM,CBNT,CBT,NvnSIn,NvnSIn,NvnSIn,1.0,0.0,VIn->MInv,VIn->DxyzChiS[dim1],SxyzIn);
			mm_d(CBRM,CBNT,CBNT,NfnI,NvnSIn,NvnSIn,1.0,0.0,OPSIn->ChiS_fI[VfIn],SxyzIn,GradxyzIn[dim1]);

			GradxyzOut[dim1] = malloc(NfnI*NvnSOut * sizeof **GradxyzOut); // free

			mm_d(CBRM,CBNT,CBT,NvnSOut,NvnSOut,NvnSOut,1.0,0.0,VOut->MInv,VOut->DxyzChiS[dim1],SxyzOut);
			mm_d(CBRM,CBNT,CBNT,NfnI,NvnSOut,NvnSOut,1.0,0.0,OPSOut->ChiS_fI[VfOut],SxyzOut,GradxyzOut[dim1]);
			array_rearrange_d(NfnI,NvnSOut,nOrdOutIn,'R',GradxyzOut[dim1]);
		}
		free(SxyzIn);
		free(SxyzOut);

		mm_dcc(CBCM,CBT,CBNT,NfnI,1,NvnSIn,1.0,0.0,OPSIn->ChiS_fI[VfIn],VIn->uhat_c,uIn_fI);

		for (dim = 0; dim < d; dim++)
			mm_dcc(CBCM,CBT,CBNT,NfnI,1,NvnSIn,1.0,0.0,GradxyzIn[dim],VIn->uhat_c,&grad_uIn_fI[NfnI*dim]);

/*
unsigned int j, IndC, dim2;
double ***GradChiS_fI;
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
*/

		if (!Boundary)
			free(detJVOut_fI);

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

		// Multiply by area element
		for (n = 0; n < NfnI; n++) {
			nqNum_fI[n] *= detJF_fI[n];
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
		mm_dcc(CBCM,CBT,CBNT,NvnSIn,1,NfnI,-1.0,1.0,OPSIn->I_Weak_FF[VfIn],nqNum_fI,RHSIn);

		// Exterior FACET
		if (!Boundary) {
			// RHS

			// Use -ve normal for opposite VOLUME
			for (n = 0; n < NfnI; n++)
				nqNum_fI[n] *= -1.0;

			// Rearrange nqNum to match node ordering from VOut
			array_rearrange_cmplx(NfnI,1,nOrdInOut,'R',nqNum_fI);

			mm_dcc(CBCM,CBT,CBNT,NvnSOut,1,NfnI,-1.0,1.0,OPSOut->I_Weak_FF[VfOut],nqNum_fI,RHSOut);
		}
		free(nqNum_fI);
	}
	free(OPSIn);
	free(OPSOut);
}
