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

#include "array_print.h" // ToBeDeleted

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
	double       *ChiS_vI, **D_Weak, **ChiS_fI, **I_Weak_FF, *I_Weak_VV, *I_vG_vI, ***GradChiS_fI, **I_vC_fI;
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
		OPS->I_Weak_VV = ELEMENT->Is_Weak_VV[P][P][0];
		OPS->D_Weak    = ELEMENT->Ds_Weak_VV[P][P][0];
		OPS->I_vG_vI   = ELEMENT->I_vGs_vIs[1][P][0];
	} else {
		OPS->NvnI = ELEMENT->NvnIc[P];

		OPS->ChiS_vI   = ELEMENT->ChiS_vIc[P][P][0];
		OPS->I_Weak_VV = ELEMENT->Ic_Weak_VV[P][P][0];
		OPS->D_Weak    = ELEMENT->Dc_Weak_VV[P][P][0];
		OPS->I_vG_vI   = ELEMENT->I_vGc_vIc[P][P][0];
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
	unsigned int   dim1, NvnI, NvnS;
	double         *ChiS_vI, *MInv, **DxyzChiS, *Sxyz;
	double complex *u_vI, *qhat;

	struct S_OPERATORS *OPS;
	struct S_VOLUME    *VOLUME;

	OPS = malloc(sizeof *OPS); // free

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		init_ops(OPS,VOLUME);

		NvnI = OPS->NvnI;
		NvnS = OPS->NvnS;

		ChiS_vI = OPS->ChiS_vI;

		// Obtain u_vI
		u_vI = malloc(NvnI * sizeof *u_vI); // free
		mm_dcc(CBCM,CBT,CBNT,NvnI,1,NvnS,1.0,0.0,ChiS_vI,VOLUME->uhat_c,u_vI);

		MInv     = VOLUME->MInv;
		DxyzChiS = VOLUME->DxyzChiS;

		// Compute RHS term
		for (dim1 = 0; dim1 < d; dim1++) {
			Sxyz = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS,NvnS,NvnS,-1.0,MInv,DxyzChiS[dim1]); // keep

			// RHS
			qhat = malloc(NvnS*1 * sizeof *qhat); // keep
			mm_dcc(CBCM,CBT,CBNT,NvnS,1,NvnS,1.0,0.0,Sxyz,VOLUME->uhat_c,qhat);
			VOLUME->qhat_c[dim1] = qhat;
		}
		free(u_vI);
	}
	free(OPS);
}

static void boundary_Dirichlet_c(const unsigned int Nn, const unsigned int Nel, double *XYZ, double complex *uL,
                                 double complex *uR)
{
	unsigned int n;
	double       *uL_d, *uR_d;

	if (Nel != 1)
		printf("Error: Vectorization unsupported.\n"), EXIT_MSG;

	uL_d = malloc(Nn * sizeof *uL_d); // free
	uR_d = malloc(Nn * sizeof *uR_d); // free
	for (n = 0; n < Nn; n++) {
		uL_d[n] = creal(uL[n]);
		uR_d[n] = creal(uR[n]);
	}

	boundary_Dirichlet(Nn,Nel,XYZ,uL_d,uR_d);
	free(uL_d);

	for (n = 0; n < Nn; n++)
		uR[n] = uR_d[n];

	free(uR_d);
}

static void trace_IP_c(const unsigned int Nn, const unsigned int Nel, double complex *uL, double complex *uR,
                       double complex *uNum)
{
	unsigned int n, NnTotal;
	double       *uL_d, *uR_d, *uNum_d;

	NnTotal = Nn*Nel;

	uL_d   = malloc(NnTotal * sizeof *uL_d);   // free
	uR_d   = malloc(NnTotal * sizeof *uR_d);   // free
	uNum_d = malloc(NnTotal * sizeof *uNum_d); // free

	for (n = 0; n < NnTotal; n++) {
		uL_d[n] = creal(uL[n]);
		uR_d[n] = creal(uR[n]);
	}
	trace_IP(Nn,Nel,uL_d,uR_d,uNum_d);
	free(uL_d);
	free(uR_d);

	for (n = 0; n < NnTotal; n++)
		uNum[n] = uNum_d[n];
	free(uNum_d);
}

void compute_qhat_FACET_c(void)
{
	// Initialize DB Parameters
	unsigned int d               = DB.d,
	             ViscousFluxType = DB.ViscousFluxType;

	// Standard datatypes
	unsigned int   n, dim, i, j, iMax,
	               NfnI, NvnSIn, NvnSOut,
	               VfIn, VfOut, fIn, EclassIn, IndFType, BC, Boundary,
	               *nOrdOutIn, *nOrdInOut;
	double         nJ, *n_fI, *detJF_fI,
	               *I_FF, *MInvI_FF, *ChiS_fI, *ChiS_fIInOut;
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
		uNum_fI = malloc(NfnI * sizeof *uNum_fI); // free
		switch (ViscousFluxType) {
		case FLUX_IP:
			trace_IP_c(NfnI,1,uIn_fI,uOut_fIIn,uNum_fI);
			break;
		default:
			printf("Error: Unsupported ViscousFluxType.\n"), EXIT_MSG;
			break;
		}
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
		MInvI_FF = mm_Alloc_d(CBRM,CBNT,CBNT,NvnSIn,NfnI,NvnSIn,1.0,VIn->MInv,I_FF); // free

		for (dim = 0; dim < d; dim++) {
			// RHSIn
			qhatIn = malloc(NvnSIn*1 * sizeof *qhatIn); // keep
			mm_dcc(CBCM,CBT,CBNT,NvnSIn,1,NfnI,1.0,0.0,MInvI_FF,&nuNum_fI[NfnI*dim],qhatIn);
			FACET->qhatIn_c[dim] = qhatIn;
		}

		// Exterior VOLUME
		if (!Boundary) {
			ChiS_fI = OPSOut->ChiS_fI[VfOut];

			free(MInvI_FF);

			// Use "-ve" normal for opposite VOLUME
			for (i = 0, iMax = NfnI*d; i < iMax; i++) {
				nuNum_fI[i] *= -1.0;
			}

			// Rearrange numerical trace to match node ordering from opposite VOLUME
			array_rearrange_cmplx(NfnI,d,nOrdInOut,nuNum_fI);

			I_FF     = OPSOut->I_Weak_FF[VfOut];
			MInvI_FF = mm_Alloc_d(CBRM,CBNT,CBNT,NvnSOut,NfnI,NvnSOut,1.0,VOut->MInv,I_FF); // free

			ChiS_fIInOut = malloc(NvnSOut*NfnI * sizeof *ChiS_fIInOut); // free

			ChiS_fI = OPSIn->ChiS_fI[VfIn];
			for (i = 0; i < NfnI; i++) {
			for (j = 0; j < NvnSIn; j++) {
				ChiS_fIInOut[i*NvnSIn+j] = ChiS_fI[nOrdInOut[i]*NvnSIn+j];
			}}

			for (dim = 0; dim < d; dim++) {
				// RHSOut
				qhatOut = malloc(NvnSOut*1 * sizeof *qhatOut); // keep
				mm_dcc(CBCM,CBT,CBNT,NvnSOut,1,NfnI,1.0,0.0,MInvI_FF,&nuNum_fI[NfnI*dim],qhatOut);
				FACET->qhatOut_c[dim] = qhatOut;
			}
			free(ChiS_fIInOut);
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

		for (dim = 0; dim < d; dim++) {
			free(FACET->qhatIn_c[dim]);
			free(FACET->qhatOut_c[dim]);
		}
	}
}

static void compute_source_c(const unsigned int Nn, double *XYZ, double complex *source)
{
	unsigned int   n;
	double         *source_d;
//	double complex *tmp;

	// silence
//	tmp = source; source = tmp;

	source_d = malloc(Nn * sizeof *source_d); // free
	compute_source(Nn,XYZ,source_d);

	for (n = 0; n < Nn; n++)
		source[n] = source_d[n];
	free(source_d);
}

void compute_uhat_VOLUME_c(void)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int   i, j, dim1,
	               NvnI, NvnS;
	double         *ChiS_vI, **DxyzChiS, *detJV_vI, *I_Weak, *I_Weak_detJV_vI, *XYZ_vI;
	double complex *q_vI, *f_vI, *RHS;

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

		// RHS (Source)
		I_Weak = OPS->I_Weak_VV;
		detJV_vI = VOLUME->detJV_vI;

		I_Weak_detJV_vI = malloc(NvnS*NvnI * sizeof *I_Weak_detJV_vI); // free
		for (i = 0; i < NvnS; i++) {
		for (j = 0; j < NvnI; j++) {
			I_Weak_detJV_vI[i*NvnI+j] = I_Weak[i*NvnI+j]*detJV_vI[j];
		}}

		XYZ_vI = malloc(NvnI*d * sizeof *XYZ_vI); // free
		mm_CTN_d(NvnI,d,VOLUME->NvnG,OPS->I_vG_vI,VOLUME->XYZ,XYZ_vI);

		f_vI = malloc(NvnI * sizeof *f_vI); // free
		compute_source_c(NvnI,XYZ_vI,f_vI);
		free(XYZ_vI);

		mm_dcc(CBCM,CBT,CBNT,NvnS,1,NvnI,-1.0,1.0,I_Weak_detJV_vI,f_vI,RHS);
		free(I_Weak_detJV_vI);
		free(f_vI);
	}
	free(OPS);
}

static void flux_IP_c(const unsigned int Nn, const unsigned int Nel, double complex *uL, double complex *uR,
                      double complex *grad_uL, double complex *grad_uR, double complex *qL, double complex *qR,
                      double *h, const unsigned int P, double complex *nqNum, double *nL, const unsigned int d)
{
	unsigned int n;
	double       *uL_d, *uR_d, *grad_uL_d, *grad_uR_d, *qL_d, *qR_d, *nqNum_d;

	// silence
//	double complex *tmp;
//	tmp = qL; qL = tmp;
//	tmp = qR; qR = tmp;

	if (Nel != 1) // Vectorization not supported
		printf("Error: Unsupported Nel.\n"), EXIT_MSG;

	uL_d      = malloc(Nn * sizeof *uL_d);      // free
	uR_d      = malloc(Nn * sizeof *uR_d);      // free
	grad_uL_d = malloc(Nn * sizeof *grad_uL_d); // free
	grad_uR_d = malloc(Nn * sizeof *grad_uR_d); // free
	qL_d      = malloc(Nn * sizeof *qL_d);      // free
	qR_d      = malloc(Nn * sizeof *qR_d);      // free
	nqNum_d   = malloc(Nn * sizeof *nqNum_d);   // free

	for (n = 0; n < Nn; n++) {
		uL_d[n]      = creal(uL[n]);
		uR_d[n]      = creal(uR[n]);
		grad_uL_d[n] = creal(grad_uL[n]);
		grad_uR_d[n] = creal(grad_uR[n]);
	}

	// Remove if NULL is no longer passed for qL/qR
	if (qL) {
		for (n = 0; n < Nn; n++) {
			qL_d[n] = creal(qL[n]);
			qR_d[n] = creal(qR[n]);
		}
	}

	flux_IP(Nn,Nel,uL_d,uR_d,grad_uL_d,grad_uR_d,qL_d,qR_d,h,P,nqNum_d,nL,d);
	free(uL_d);
	free(uR_d);
	free(grad_uL_d);
	free(grad_uR_d);
	free(qL_d);
	free(qR_d);

	for (n = 0; n < Nn; n++)
		nqNum[n] = nqNum_d[n];
	free(nqNum_d);
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
	               *detJV_fI, *C_fI, *h, *n_fI, *detJF_fI, *C_vC;
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

		detJV_fI = FACET->detJV_fI;
		detJF_fI = FACET->detJF_fI;
		h = malloc(NfnI * sizeof *h); // free
		for (n = 0; n < NfnI; n++)
			h[n] = detJV_fI[n]/detJF_fI[n];

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
					GradxyzIn[dim1][n*NvnSIn+j] /= detJV_fI[n];
			}
		}

		for (dim = 0; dim < d; dim++)
			mm_dcc(CBCM,CBT,CBNT,NfnI,1,NvnSIn,1.0,0.0,GradxyzIn[dim],VIn->uhat_c,&grad_uIn_fI[NfnI*dim]);

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
					GradxyzOut[dim1][n*NvnSOut+j] /= detJV_fI[n];
			}
		}

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
		nqNum_fI = malloc(NfnI * sizeof *nqNum_fI); // free

		switch (ViscousFluxType) {
		case FLUX_IP:
			flux_IP_c(NfnI,1,uIn_fI,uOut_fIIn,grad_uIn_fI,grad_uOut_fIIn,NULL,NULL,h,FACET->P,nqNum_fI,n_fI,d);
			break;
		default:
			printf("Error: Unsupported PoissonFluxType.\n"), EXIT_MSG;
			break;
		}

		// Multiply by area element
		for (n = 0; n < NfnI; n++) {
			nqNum_fI[n] *= detJF_fI[n];
		}

		// Finalize FACET RHS terms
		RHSIn  = calloc(NvnSIn  , sizeof *RHSIn);  // keep (requires external free)
		RHSOut = calloc(NvnSOut , sizeof *RHSOut); // keep (requires external free)
		FACET->RHSIn_c  = RHSIn;
		FACET->RHSOut_c = RHSOut;

		// Interior FACET
		mm_dcc(CBCM,CBT,CBNT,NvnSIn,1,NfnI,1.0,0.0,OPSIn->I_Weak_FF[VfIn],nqNum_fI,RHSIn);

		// Exterior FACET
		if (!Boundary) {
			// RHS

			// Use -ve normal for opposite VOLUME
			for (n = 0; n < NfnI; n++)
				nqNum_fI[n] *= -1.0;

			// Rearrange nqNum to match node ordering from VOut
			array_rearrange_cmplx(NfnI,1,nOrdInOut,nqNum_fI);

			mm_dcc(CBCM,CBT,CBNT,NvnSOut,1,NfnI,1.0,0.0,OPSOut->I_Weak_FF[VfOut],nqNum_fI,RHSOut);
		}
	}
}
