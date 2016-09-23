// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "solver_poisson.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "S_FACET.h"

// When errors are fixed, comment all headers below and include ALL that are necessary (ToBeDeleted)
#include "element_functions.h"
#include "matrix_functions.h"
#include "exact_solutions.h"



#include "array_print.h" // ToBeDeleted

/*
 *	Purpose:
 *		Perform the implicit solve for the Poisson equation.
 *
 *	Comments:
 *		CHECK FOR MEMORY LEAKS. ToBeDeleted
 *
 *		Many of the RHS terms computed are 0. They are included as they are used to check the linearization. Further,
 *		the computational cost is dominated by the global system solve making this additional cost negligible.
 *		When finished with the implicit functions, include these contributions, include these redundant terms only in
 *		the explicit functions. (ToBeDeleted)
 *
 *	Notation:
 *
 *	References:
 */

struct S_OPERATORS {
	unsigned int NvnS, NvnI, NfnI,
	             *nOrdOutIn, *nOrdInOut;
	double       *ChiS_vI, **D_Weak, **ChiS_fI, **I_Weak_FF;
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

		OPS->ChiS_vI = ELEMENT->ChiS_vIs[P][P][0];
		OPS->D_Weak  = ELEMENT->Ds_Weak_VV[P][P][0];
	} else {
		OPS->NvnI = ELEMENT->NvnIc[P];

		OPS->ChiS_vI = ELEMENT->ChiS_vIc[P][P][0];
		OPS->D_Weak  = ELEMENT->Dc_Weak_VV[P][P][0];
	}
}

static void init_opsF(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const struct S_FACET *FACET,
                      const unsigned int IndFType)
{
	// Standard datatypes
	unsigned int PV, PF, Vtype, FtypeInt, IndOrdOutIn, IndOrdInOut;

	struct S_ELEMENT *ELEMENT, *ELEMENT_FACET;

	PV     = VOLUME->P;
	PF     = FACET->P;
	Vtype  = VOLUME->type;

	FtypeInt = FACET->typeInt;
	IndOrdOutIn = FACET->IndOrdOutIn;
	IndOrdInOut = FACET->IndOrdInOut;

	ELEMENT       = get_ELEMENT_type(Vtype);
	ELEMENT_FACET = get_ELEMENT_FACET(Vtype,IndFType);

	if (FtypeInt == 's') {
		// Straight FACET Integration
		OPS->NfnI = ELEMENT->NfnIs[PF][IndFType];

		OPS->ChiS_fI   = ELEMENT->ChiS_fIs[PV][PF];
		OPS->I_Weak_FF = ELEMENT->Is_Weak_FF[PV][PF];

		OPS->nOrdInOut = ELEMENT_FACET->nOrd_fIs[PF][IndOrdInOut];
		OPS->nOrdOutIn = ELEMENT_FACET->nOrd_fIs[PF][IndOrdOutIn];
	} else {
		// Curved FACET Integration
		OPS->NfnI = ELEMENT->NfnIc[PF][IndFType];

		OPS->ChiS_fI   = ELEMENT->ChiS_fIc[PV][PF];
		OPS->I_Weak_FF = ELEMENT->Ic_Weak_FF[PV][PF];

		OPS->nOrdInOut = ELEMENT_FACET->nOrd_fIc[PF][IndOrdInOut];
		OPS->nOrdOutIn = ELEMENT_FACET->nOrd_fIc[PF][IndOrdOutIn];
	}
}

static void compute_qhat_VOLUME(void)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int i, j, dim1, dim2, 
	             IndD, IndC,
	             NvnI, NvnS;
	double       *u_vI;
	double       *ChiS_vI, *MInv, **D, *C_vI, *Dxyz, **DxyzChiS, *Sxyz;

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
		mm_CTN_d(NvnI,1,NvnS,ChiS_vI,VOLUME->uhat,u_vI);

		MInv = VOLUME->MInv;
		C_vI = VOLUME->C_vI;

		// Construct physical derivative operator matrices
		D = OPS->D_Weak;

		Dxyz     = malloc(NvnS*NvnI * sizeof *Dxyz); // free
		DxyzChiS = malloc(d * sizeof *Dxyz); // keep
		for (dim1 = 0; dim1 < d; dim1++)
			DxyzChiS[dim1] = malloc(NvnS*NvnS * sizeof *DxyzChiS[dim1]); // free

		for (dim1 = 0; dim1 < d; dim1++) {
			IndD = 0;
			for (dim2 = 0; dim2 < d; dim2++) {
				IndC = (dim1+dim2*d)*NvnI;
				for (i = 0; i < NvnS; i++) {
					for (j = 0; j < NvnI; j++)
						Dxyz[IndD+j] = D[dim2][IndD+j]*C_vI[IndC+j];
					IndD += NvnI;
				}
			}
			DxyzChiS[dim1] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS,NvnS,NvnI,1.0,Dxyz,ChiS_vI); // keep
		}

// Note: Storage of DxyzChiS is not necessary if cofactor terms are included in uhat^ref and qhat^ref as in the paper.
// ToBeDeleted
		VOLUME->DxyzChiS = DxyzChiS;

		// Compute RHS and LHS terms
		for (dim1 = 0; dim1 < d; dim1++) {
			Sxyz = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS,NvnI,NvnS,-1.0,MInv,DxyzChiS[dim1]); // keep

			// RHS
			VOLUME->qhat[dim1] = mm_Alloc_d(CBCM,CBT,CBNT,NvnS,1,NvnI,1.0,Sxyz,VOLUME->uhat); // keep

			// LHS
			VOLUME->qhat_uhat[dim1] = Sxyz;
		}
		free(u_vI);
	}
	free(OPS);
}

static void boundary_Dirichlet(const unsigned int Nn, const unsigned int Nel, double *XYZ, double *uB)
{
	// Is this right? Or should uB == -uIn? ToBeDeleted
	compute_exact_solution(Nn*Nel,XYZ,uB,NULL,0);
}

static void jacobian_boundary_Dirichlet(const unsigned int Nn, const unsigned int Nel, double *duOutduIn)
{
	unsigned int n;

	if (Nel != 1)
		printf("Error: Vectorization unsupported.\n"), EXIT_MSG;

	for (n = 0; n < Nn; n++)
		duOutduIn[n] = 0.0;
}

static void trace_IP(const unsigned int Nn, const unsigned int Nel, double *uL, double *uR, double *uNum)
{
	/*
	 *	Comments:
	 *		(T)race (flux) for the (I)nternal (P)enalty method:
	 *			uNum = 1/2*(uL+uR)
	 */

	unsigned int n, NnTotal;

	NnTotal = Nn*Nel;

	for (n = 0; n < NnTotal; n++)
		uNum[n] = 0.5*(uL[n]+uR[n]);
}

static void jacobian_trace_IP(const unsigned int Nn, const unsigned int Nel, double *duNumdu, const char side)
{
	unsigned int n, NnTotal;

	if (side != 'L' || side != 'R')
		printf("Error: Invalid side.\n"), EXIT_MSG;

	NnTotal = Nn*Nel;

	for (n = 0; n < NnTotal; n++)
		duNumdu[n] = 0.5;
}

static void compute_qhat_FACET(void)
{
	// Initialize DB Parameters
	unsigned int d               = DB.d,
	             ViscousFluxType = DB.ViscousFluxType;

	// Standard datatypes
	unsigned int n, dim, i, j;
	unsigned int NfnI, NvnSIn, NvnSOut;
	unsigned int VfIn, VfOut, fIn, EclassIn, IndFType, BC, Boundary;
	unsigned int *nOrdOutIn, *nOrdInOut;
	double       nJ, *n_fI, *detJF_fI;
	double       *uIn_fI, *uOut_fIIn, *uOut_fI, *uNum_fI, *duNumduIn_fI, *duNumduOut_fI, *duOutduIn,
	             *nuNum_fI, *dnuNumduIn_fI, *dnuNumduOut_fI;
	double       *I_FF, *MInvI_FF;

	struct S_OPERATORS *OPSIn, *OPSOut;
	struct S_VOLUME    *VIn, *VOut;
	struct S_FACET     *FACET;

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
		NvnSIn = OPSIn->NvnS;

		nOrdOutIn = OPSIn->nOrdOutIn;
		nOrdInOut = OPSIn->nOrdInOut;

		// Compute uIn_fI
		uIn_fI = malloc(NfnI * sizeof *uIn_fI); // free
		mm_CTN_d(NfnI,1,NvnSIn,OPSIn->ChiS_fI[VfIn],VIn->uhat,uIn_fI);

		// Compute_uOut_fI (Taking BCs into account if applicable)
		uOut_fIIn = malloc(NfnI * sizeof *uOut_fIIn); // free
		if (!Boundary) {
			uOut_fI = malloc(NfnI * sizeof *uOut_fI); // free
			mm_CTN_d(NfnI,1,NvnSOut,OPSOut->ChiS_fI[VfOut],VOut->uhat,uOut_fI);

			// Reorder uOut_fI to correspond to uIn_fI
			for (n = 0; n < NfnI; n++)
				uOut_fIIn[n] = uOut_fI[nOrdOutIn[n]];
			free(uOut_fI);
		} else {
			if (BC % BC_STEP_SC == BC_DIRICHLET) {
				boundary_Dirichlet(NfnI,1,FACET->XYZ_fI,uOut_fIIn);
			} else if (BC % BC_STEP_SC == BC_NEUMANN) {
				printf("Add support.\n"), EXIT_MSG;
			} else {
				printf("Error: Unsupported BC.\n"), EXIT_MSG;
			}
		}

		// Compute numerical trace and its Jacobians
		switch (ViscousFluxType) {
		case FLUX_IP:
			trace_IP(NfnI,1,uIn_fI,uOut_fIIn,uNum_fI);
			jacobian_trace_IP(NfnI,1,duNumduIn_fI,'L');
			jacobian_trace_IP(NfnI,1,duNumduOut_fI,'R');
			break;
		default:
			printf("Error: Unsupported ViscousFluxType.\n"), EXIT_MSG;
			break;
		}
		free(uIn_fI);
		free(uOut_fIIn);

		// Include BC information in duNumduIn_fI if on a boundary
		if (Boundary) {
			duOutduIn = malloc(NfnI * sizeof *duOutduIn); // free

			if (BC % BC_STEP_SC == BC_DIRICHLET)
				jacobian_boundary_Dirichlet(NfnI,1,duOutduIn);
			else if (BC % BC_STEP_SC == BC_NEUMANN)
				printf("Error: Add support.\n"), EXIT_MSG;
			else
				printf("Error: Unsupported BC.\n"), EXIT_MSG;

			for (n = 0; n < NfnI; n++)
				duNumduIn_fI[n] += duNumduOut_fI[n]*duOutduIn[n];

if (BC % BC_STEP_SC == BC_DIRICHLET) { // ToBeDeleted
printf("Should be zero.\n");
array_print_d(NfnI,1,duNumduIn_fI,'R');
EXIT_MSG;
}

			free(duOutduIn);
		}

		// Multiply uNum and its Jacobian by the normal vector and area element
		nuNum_fI       = malloc(NfnI*d * sizeof *nuNum_fI); // tbd
		dnuNumduIn_fI  = malloc(NfnI*d * sizeof *dnuNumduIn_fI); // tbd
		dnuNumduOut_fI = malloc(NfnI*d * sizeof *dnuNumduOut_fI); // tbd
		for (dim = 0; dim < d; dim++) {
		for (n = 0; n < NfnI; n++) {
			nJ = n_fI[n*d+dim]*detJF_fI[n];
			nuNum_fI[dim*NfnI+n]       = uNum_fI[n]*nJ;
			dnuNumduIn_fI[dim*NfnI+n]  = duNumduIn_fI[n]*nJ;
			dnuNumduOut_fI[dim*NfnI+n] = duNumduOut_fI[n]*nJ;
		}}


		// Compute FACET RHS and LHS terms

		// Interior VOLUME

		I_FF     = OPSIn->I_Weak_FF[VfIn];
		MInvI_FF = mm_Alloc_d(CBRM,CBNT,CBNT,NvnSIn,NfnI,NvnSIn,1.0,VIn->MInv,I_FF); // free

		MInvIdnuNumdu = malloc(NvnSIn*NfnI * sizeof *MInvIdnuNumdu); // free
		for (dim = 0; dim < d; dim++) {
			// RHSIn
			qhatIn = mm_Alloc_d(CBCM,CBT,CBNT,NvnSIn,1,NfnI,1.0,MInvI_FF,&nuNum_fI[NfnI*dim]); // keep
			FACET->qhatIn[dim] = qhatIn;

			// LHS (InIn)
			for (i = 0; i < NvnSIn; i++) {
			for (j = 0; j < NvnI; j++) {
				MInvIdnuNumdu[i*NvnI+j] = MInvI_FF[i*NvnI+j]*dnuNumduIn_fI[NfnI*dim+j];
			}}

			FACET->qhat_uhatInIn[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnSIn,NvnSIn,NfnI,1.0,MInvIdnuNumdu,OPSIn->ChiS_fI[VfIn]);
		}

		// Exterior VOLUME
		if (!Boundary) {
			ChiS_fIOutIn = malloc(NvnSOut*NfnI * sizeof *ChiS_fIOutIn); // free

			ChiS_fI = OPSOut->ChiS_fI[VfOut];
			for (i = 0; i < NfnI; i++) {
			for (j = 0; j < NvnSOut; j++) {
				ChiS_fIOutIn[i*NvnSOut+j] = ChiS_fI[nOrdOutIn[i]*NvnSOut+j];
			}}

			for (dim = 0; dim < d; dim++) {
				// LHS (OutIn)
				for (i = 0; i < NvnSIn; i++) {
				for (j = 0; j < NvnI; j++) {
					MInvIdnuNumdu[i*NvnI+j] = MInvI_FF[i*NvnI+j]*dnuNumduOut_fI[NfnI*dim+j];
				}}

				FACET->qhat_uhatOutIn[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnSIn,NvnSIn,NfnI,1.0,MInvIdnuNumdu,ChiS_fIOutIn);
			}
			free(ChiS_fIOutIn);
			free(MInvIdnuNumdu);

			free(MInvI_FF);

			// Use "-ve" normal for opposite VOLUME
			for (i = 0, iMax = NfnI*d; i < iMax; i++) {
				nuNum_fI[i]      *= -1.0;
				duNumduIn_fI[i]  *= -1.0;
				duNumduOut_fI[i] *= -1.0;
			}

			// Rearrange numerical trace and its Jacobians to match node ordering from opposite VOLUME
			array_rearrange(NfnI,d,nOrdInOut,nuNum);
			array_rearrange(NfnI,d,nOrdInOut,dnuNumduIn);
			array_rearrange(NfnI,d,nOrdInOut,dnuNumduOut);

			I_FF     = OPSOut->I_Weak_FF[VfOut];
			MInvI_FF = mm_Alloc_d(CBRM,CBNT,CBNT,NvnSOut,NfnI,NvnSOut,1.0,VOut->MInv,I_FF); // free

			MInvIdnuNumdu = malloc(NvnSOut*NfnI * sizeof *MInvIdnuNumdu); // free

			ChiS_fIOutIn = malloc(NvnSOut*NfnI * sizeof *ChiS_fIOutIn); // free

			ChiS_fI = OPSIn->ChiS_fI[VfIn];
			for (i = 0; i < NfnI; i++) {
			for (j = 0; j < NvnSIn; j++) {
				ChiS_fIInOut[i*NvnSIn+j] = ChiS_fI[nOrdInOut[i]*NvnSIn+j];
			}}

			for (dim = 0; dim < d; dim++) {
				// RHSOut
				qhatOut = mm_Alloc_d(CBCM,CBT,CBNT,NvnSOut,1,NfnI,1.0,MInvI_FF,&nuNum_fI[NfnI*dim]); // keep
				FACET->qhatOut[dim] = qhatIn;

				// LHS (InOut)
				for (i = 0; i < NvnSOut; i++) {
				for (j = 0; j < NvnI; j++) {
					MInvIdnuNumdu[i*NvnI+j] = MInvI_FF[i*NvnI+j]*dnuNumduIn_fI[NfnI*dim+j];
				}}

				FACET->qhat_uhatInOut[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnSOut,NvnSIn,NfnI,1.0,MInvIdnuNumdu,ChiS_fIInOut);

				// LHS (OutOut)
				for (i = 0; i < NvnSOut; i++) {
				for (j = 0; j < NvnI; j++) {
					MInvIdnuNumdu[i*NvnI+j] = MInvI_FF[i*NvnI+j]*dnuNumduOut_fI[NfnI*dim+j];
				}}

				FACET->qhat_uhatOutOut[dim] =
					mm_Alloc_d(CBRM,CBNT,CBNT,NvnSOut,NvnSOut,NfnI,1.0,MInvIdnuNumdu,OPSOut->ChiS_fI[VfOut]);
			}
			free(ChiS_fIInOut);
		}
		free(MInvIdnuNumdu);
		free(MInvI_FF);

		free(nuNum_fI);
		free(dnuNumduIn_fI);
		free(dnuNumduOut_fI);
		free(RowTracker);
	}

	free(OPSIn);
	free(OPSOut);
}

static void finalize_qhat(void)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int dim, iMax, NvnSIn, NvnSOut;
	double       *VqhatIn_ptr, *FqhatIn_ptr, *VqhatOut_ptr, *FqhatOut_ptr;

	struct S_FACET  *FACET;
	struct S_VOLUME *VIn, *VOut;

	for (FACET = DB.FACET; FACET; FACET = FACET->next) {
		VIn    = FACET->VIn;
		NvnSIn = VIn->NvnS;

		for (dim = 0; dim < d; dim++) {
			VqhatIn_ptr  = VIn->qhat[dim];
			FqhatIn_ptr  = FACET->qhatIn[dim];

			for (iMax = NvnSIn; iMax--; )
				*VqhatIn_ptr++ += *FqhatIn_ptr;
		}

		if(!(FACET->Boundary)) {
			VOut    = FACET->VOut;
			NvnSOut = VOut->NvnS;

			for (dim = 0; dim < d; dim++) {
				VqhatOut_ptr = VOut->qhat[dim];
				FqhatOut_ptr = FACET->qhatOut[dim];

				for (iMax = NvnSOut; iMax--; )
					*VqhatOut_ptr++ += *FqhatOut_ptr;
			}
		}

		free(FACET->qhatIn);
		free(FACET->qhatOut);
	}
}

static void compute_uhat_VOLUME(void)
{
	// Standard datatypes

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
			mm_CTN_d(NvnI,1,NvnS,ChiS_vI,VOLUME->qhat[dim1],&q_vI[dim1*NvnI]);

		MInv = VOLUME->MInv;
		C_vI = VOLUME->C_vI;

		// Compute RHS and LHS terms
		DxyzChiS = VOLUME->DxyzChiS;

		// RHS
		if (VOLUME->RHS)
			free(VOLUME->RHS);
		RHS = calloc(NvnS , sizeof *RHS); // keep
		VOLUME->RHS = RHS;

		Dq = malloc(NvnS * sizeof *Dq); // free
		for (dim1 = 0; dim1 < d; dim1++) {
			mm_CTN_d(NvnS,1,NvnI,DxyzChiS[dim1],VOLUME->qhat[dim1],Dq);

			for (n = 0; n < NvnS; n++)
				RHS[n] -= Dq[n];
		}
		free(Dq);
		free(q_vI);

		// RHS (Source)
		w_vI = OPS->w_vI;
		detJV_vI = VOLUME->detJV_vI;

		for (n = 0; n < NvnI; n++)
			wdetJV_vI[n] = w_vI[n]*detJV_vI[n];

		ChiSTwdetJV_vI = malloc(NvnS*NvnI * sizeof *ChiSwdetJV_vI); // free
		for (i = 0; i < NvnS; i++) {
		for (j = 0; j < NvnI; j++) {
			ChiSTwdetJV_vI[i*NvnI+j] = ChiS_vI[i+j*NvnS]*wdetJV_vI[j];
		}}

		XYZ_vI = malloc(NvnI*d * sizeof *XYZ_vI); // free
		mm_CTN_d(NvnI,d,VOLUME->NvnG,OPS->I_vG_vI,VOLUME->XYZ,XYZ_vI);

		compute_source(NfnI,XYZ_vI,f_vI);
		free(XYZ_vI);

		RHSs = malloc(NvnS * sizeof *RHSs); // free
		mm_CTN_d(NvnS,1,NvnI,ChiSTwdetJV_vI,f_vI,RHSs);
		free(ChiSwdetJV_vI);

		for (n = 0; n < NvnS; n++)
			RHS[n] -= RHSs[n];

		free(RHSs);

		// LHS
		if (VOLUME->LHS)
			free(VOLUME->LHS);
		LHS = calloc(NvnS*NvnS , sizeof *LHS); // keep
		VOLUME->LHS = LHS;

		for (dim1 = 0; dim1 < d; dim1++) {
			LHSd = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS,NvnS,NvnS,-1.0,DxyzChiS,VOLUME->qhat_uhat[dim1]); // free

			for (i = 0, iMax = NvnS*NvnS; i < iMax; i++)
				LHS[i] += LHSd[i];
		
			free(LHSd);
		}

	}
	
	free(OPS);
}

static void flux_IP(const unsigned int Nn, const unsigned int Nel, double *uIn, double *uOut, double *grad_uIn,
                    double *grad_uOut, double *qIn, double *qOut, double *h, const unsigned int P, double *nqNum,
                    double *nIn, const unsigned int d)
{
	unsigned int n, dim;
	double       tau;

	if (Nel != 1) // Vectorization not supported
		printf("Error: Unsupported Nel.\n"), EXIT_MSG;

	for (n = 0; n < Nn; n++) {
		tau = CONST_IP*(P+1)*(P+1)/h[n];

		for (dim = 0; dim < d; dim++)
			nqNum[n] = nIn[n*d+dim]*(0.5*(grad_uIn[Nn*dim+n]+grad_uOut[Nn*dim+n]) - tau*(nIn[n*d+dim]*(uIn[n]-uOut[n])));
	}
}

static void jacobian_flux_coef(const unsigned int Nn, const unsigned int Nel, const double *nIn, const double *h,
                               const unsigned int P, const double *gradu_avg, double *u_jump, double *q_avg, double *q_jump,
                               const unsigned int d, char *flux_type, char side)
{
	/*
	 *	Comments:
	 *		Likely a good idea to define flux_IP similarly (i.e. simply returning coefficients for various terms. Then
	 *		only one of the two functions would be required.
	 *
	 *	References:
	 *		Add references for definition of flux coefficients (e.g. Hesthaven(2008) Table 7.3) ToBeModified
	 */

	unsigned int dim, n;
	double       tau;

	if (Nel != 1)
		printf("Error: Unsupported Nel.\n"), EXIT_MSG;

	if (strstr(flux_type,"IP")) {
		for (dim = 0; dim < d; dim++) {
		for (n = 0; n < Nn; n++) {
			tau = CONST_IP*(P+1)*(P+1)/h[n];

			gradu_avg[Nn*dim+n] = 0.5*nIn[n*d+dim];
			q_avg[Nn*dim+n]     = 0.0;

			u_jump[Nn*dim+n] = tau*nIn[n*d+dim];
			q_jump[Nn*dim+n] = 0.0;

			if (side == 'L') {
				u_jump[Nn*dim+n] *= -1.0;
				q_jump[Nn*dim+n] *= -1.0;
			}
		}}
	} else {
		printf("Error: Unsupported flux_type.\n"), EXIT_MSG;
	}
}

static void compute_uhat_FACET()
{
	// Standard datatypes

	struct S_OPERATORS *OPSIn, *OPSOut;
	struct S_FACET     *FACET;
	struct S_VOLUME    *VIn, *VOut;

	for (FACET = DB.FACET; FACET; FACET = FACET->next) {
		VIn = FACET->VIn;

		VOut = FACET->VOut;

		Boundary = FACET->Boundary;



		DxyzChiS = VIn->DxyzChiS;

		LHSInIn   = calloc(NvnSIn*NvnSIn   , sizeof *LHSInIn);   // keep
		LHSOutIn  = calloc(NvnSIn*NvnSOut  , sizeof *LHSOutIn);  // keep
		LHSInOut  = calloc(NvnSOut*NvnSIn  , sizeof *LHSInOut);  // keep
		LHSOutOut = calloc(NvnSOut*NvnSOut , sizeof *LHSOutOut); // keep

		FACET->LHSInIn   = LHSInIn;
		FACET->LHSOutIn  = LHSOutIn;
		FACET->LHSInOut  = LHSInOut;
		FACET->LHSOutOut = LHSOutOut;

		for (dim = 0; dim < d; dim++) {
			// LHS - VOLUME terms
			mm_d(CBRM,CBNT,CBNT,NvnSIn, NvnSIn, NvnSIn, -1.0,1.0,DxyzChiS[dim],FACET->qhat_uhatInIn[dim],  LHSInIn);
// Something relating to whether this is a boundary VOLUME here? (ToBeDeleted)
			mm_d(CBRM,CBNT,CBNT,NvnSIn, NvnSOut,NvnSIn, -1.0,1.0,DxyzChiS[dim],FACET->qhat_uhatOutIn[dim], LHSOutIn);
			mm_d(CBRM,CBNT,CBNT,NvnSOut,NvnSIn, NvnSOut,-1.0,1.0,DxyzChiS[dim],FACET->qhat_uhatInOut[dim], LHSInOut);
			mm_d(CBRM,CBNT,CBNT,NvnSOut,NvnSOut,NvnSOut,-1.0,1.0,DxyzChiS[dim],FACET->qhat_uhatOutOut[dim],LHSOutOut);
		}

		// Compute uIn_fI and gradu_In_fI
		detJV_fI = FACET->detJV_fI;

		uIn_fI         = malloc(NfnI   * sizeof *uIn_fI);         // tbd
		gradref_uIn_fI = malloc(NfnI*d * sizeof *gradref_uIn_fI); // free
		grad_uIn_fI    = calloc(NfnI*d , sizeof *graduIn_fI);     // tbd

		mm_CTN_d(NfnI,1,NvnSIn,OPSIn->ChiS_fI[VfIn],VIn->uhat,uIn_fI);

		for (dim = 0; dim < d; dim++)
			mm_CTN_d(NfnI,1,NvnSIn,OPSIn->GradChiS_fI[VfIn][dim],VIn->uhat,&gradref_uIn_fI[NfnI*dim]);

		for (dim1 = 0; dim1 < d; dim1++) {
			for (dim2 = 0; dim2 < d; dim2++) {
				IndC = (dim1+dim2*d)*NfnI;
				for (n = 0; n < NfnI; n++)
					grad_uIn_fI[dim1*NfnI+n] += gradref_uIn_fI[dim2*NfnI+n]*C_fI[IndC+n];
			}
			for (n = 0; n < NfnI; n++)
				grad_uIn_fI[dim1*NfnI+n] /= detJV_fI[n];
		}
		free(gradref_uIn_fI);

array_print_d(NfnI,d,grad_uIn_fI,'C');

// Alternate computation (first assembling operator)
// Don't need gradref_uIn_fI in this case (ToBeDeleted)
		GradChiS_fI = OPSIn->GradChiS_fI;

		GradxyzIn = malloc(d * sizeof *GradxyzIn); // tbd
		for (dim1 = 0; dim1 < d; dim1++) {
			GradxyzIn[dim1] = calloc(NfnI*NvnSIn , sizeof **GradxyzIn); // tbd

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
			mm_CTN_d(NfnI,1,NvnSIn,GradxyzIn[dim],VIn->uhat,&grad_uIn_fI[NfnI*dim]);

array_print_d(NfnI,d,graduIn_fI,'C');



		// Compute_uOut_fI (Taking BCs into account if applicable)
		uOut_fIIn         = malloc(NfnI   * sizeof *uOut_fIIn);       // tbd
		gradref_uOut_fIIn = malloc(NfnI*d * sizeof *gradref_uOut_fI); // free
		grad_uOut_fIIn    = calloc(NfnI*d , sizeof *graduOut_fI);     // tbd

		if (!Boundary) {
			uOut_fI         = malloc(NfnI   * sizeof *uOut_fI);         // free
			gradref_uOut_fI = malloc(NfnI*d * sizeof *gradref_uOut_fI); // free

			mm_CTN_d(NfnI,1,NvnSOut,OPSOut->ChiS_fI[VfOut],VOut->uhat,uOut_fI);

			for (dim = 0; dim < d; dim++)
				mm_CTN_d(NfnI,1,NvnSOut,OPSOut->GradChiS_fI[VfOut][dim],VOut->uhat,&gradref_uOut_fI[NfnI*dim]);

			// Reorder uOut_fI and gradref_uOut_fI to correspond to inner VOLUME ordering
			for (n = 0; n < NfnI; n++) {
				uOut_fIIn[n]         = uOut_fI[nOrdOutIn[n]];
				gradref_uOut_fIIn[n] = gradref_uOut_fI[nOrdOutIn[n]];
			}
			free(uOut_fI);
			free(gradref_uOut_fI);

			for (dim1 = 0; dim1 < d; dim1++) {
				for (dim2 = 0; dim2 < d; dim2++) {
					IndC = (dim1+dim2*d)*NfnI;
					for (n = 0; n < NfnI; n++)
						grad_uOut_fIIn[dim1*NfnI+n] += gradref_uOut_fIIn[dim2*NfnI+n]*C_fI[IndC+n];
				}
				for (n = 0; n < NfnI; n++)
					grad_uOut_fIIn[dim1*NfnI+n] /= detJV_fI[n];
			}
			free(gradref_uOut_fIIn);

array_print_d(NfnI,d,grad_uOut_fIIn,'C');
// Alternate (ToBeDeleted (comment))
// note: already re-arranged
		GradChiS_fI = OPSOut->GradChiS_fI;

		GradxyzOut = malloc(d * sizeof *GradxyzOut); // tbd
		for (dim1 = 0; dim1 < d; dim1++) {
			GradxyzOut[dim1] = calloc(NfnI*NvnSIn , sizeof **GradxyzOut); // tbd

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

		for (dim = 0; dim < d; dim++)
			mm_CTN_d(NfnI,1,NvnSOut,GradxyzOut[dim],VOut->uhat,&grad_uOut_fIIn[NfnI*dim]);

array_print_d(NfnI,d,graduOut_fIIn,'C');
		} else {
			if (BC % BC_STEP_SC == BC_DIRICHLET) {
				boundary_Dirichlet(NfnI,1,FACET->XYZ_fI,uOut_fIIn);
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

		// Compute numerical flux and its Jacobians
		nqNum_fI = malloc(NfnI * sizeof *nqNum_fI); // tbd

		switch (ViscousFluxType) {
		case FLUX_IP:
			flux_IP(NfnI,1,uIn_fI,uOut_fIIn,graduIn_fI,graduOut_fIIn,NULL,NULL,h,FACET->P,nqNum_fI,n_fI,d);
			break;
		default:
			printf("Error: Unsupported ViscousFluxType.\n"), EXIT_MSG;
			break;
		}

// If using ViscousFluxType == "IP", what is done below is quite inefficient as there is no dependence on q. It is very
// general however and supports many different FluxTypes.
		dnqNumduhatIn_fI    = calloc(NfnI*NvnSIn  , sizeof *dnqNumduhatIn_fI);    // tbd
		dnqNumduhatOut_fIIn = calloc(NfnI*NvnSOut , sizeof *dnqNumduhatOut_fIIn); // tbd


		// InIn contributions
		jacobian_flux_coef(NfnI,1,n_fI,h,FACET->P,gradu_avg,u_jump,q_avg,q_jump,d,"IP",'L');
		ChiS_fI = OPSIn->ChiS_fI[VfIn];

		for (dim = 0; dim < d; dim++) {
			q_uhatV = mm_Alloc_d(CBRM,CBNT,CBNT,NfnI,NvnSIn,NvnSIn,1.0,ChiS_fI,VIn->qhat_uhat[dim]);       // free
			q_uhatF = mm_Alloc_d(CBRM,CBNT,CBNT,NfnI,NvnSIn,NvnSIn,1.0,ChiS_fI,FACET->qhat_uhatInIn[dim]); // free
			for (n = 0; n < NfnI; n++) {
			for (j = 0; j < NvnSIn; j++) {
				dnqNumduhatIn_fI[n*NvnSIn+j] // row-major
					+= u_avg[NfnI*dim+n]*GradxyzIn[dim][n*NvnSIn+j]
					+  u_jump[NfnI*dim+n]*ChiS_fI[n*NvnSIn+j]
					+  (q_avg[NfnI*dim+n]+q_jump[NfnI*dim+n])*(q_uhatV[n*NvnSIn+j]+q_uhatF[n*NvnSIn+j]);
			}}
			free(q_uhatV);
			free(q_uhatF);
		}

		// If on a boundary, there is no qhat contribution other than qhat_uhatInIn. However, there is a still a
		// contribution from terms depending on the solution, which are added to dnqNumduhatIn_fI using the chain
		// rule.
		// Note: Both dnqNumduhatIn and dnqNumduhatOut are ordered corresponding to VIn. 

		// OutOut contribution (u)
		jacobian_flux_coef(NfnI,1,n_fI,h,FACET->P,gradu_avg,u_jump,q_avg,q_jump,d,"IP",'R');
		ChiS_fI_std = OPSOut->ChiS_fI[VfOut];
		ChiS_fI = malloc(NfnI*NvnSOut * sizeof *ChiS_fI); // free

		// Reordering
		for (i = 0; i < NfnI; i++) {
		for (j = 0; j < NvnSOut; j++) {
			ChiS_fI[i*NvnSOut+j] = ChiS_fI_std[nOrdOutIn[i]*NvnSOut+j];
		}}

		for (dim = 0; dim < d; dim++) {
			for (n = 0; n < NfnI; n++) {
			for (j = 0; j < NvnSOut; j++) {
				dnqNumduhatOut_fI[n*NvnSOut+j] += u_avg[NfnI*dim+n]*GradxyzOut[dim][n*NvnSOut+j]
				                               +  u_jump[NfnI*dim+n]*ChiS_fI[n*NvnSOut+j];
			}}
		}

		if (!Boundary) {
			// InOut contribution
			jacobian_flux_coef(NfnI,1,n_fI,h,FACET->P,gradu_avg,u_jump,q_avg,q_jump,d,"IP",'R');

			// ChiS_fI reordered above.
			for (dim = 0; dim < d; dim++) {
				q_uhatF = mm_Alloc_d(CBRM,CBNT,CBNT,NfnI,NvnSIn,NvnSOut,1.0,ChiS_fI,FACET->qhat_uhatInOut[dim]); // free
				for (n = 0; n < NfnI; n++) {
				for (j = 0; j < NvnSIn; j++) {
					dnqNumduhatIn_fI[n*NvnSIn+j] += (q_avg[NfnI*dim+n]+q_jump[NfnI*dim+n])*(q_uhatF[n*NvnSIn+j]);
				}}
				free(q_uhatF);
			}
			free(ChiS_fI);

			// OutIn contribution
			jacobian_flux_coef(NfnI,1,n_fI,h,FACET->P,gradu_avg,u_jump,q_avg,q_jump,d,"IP",'L');
			ChiS_fI = OPSIn->ChiS_fI[VfIn];

			for (dim = 0; dim < d; dim++) {
				q_uhatF = mm_Alloc_d(CBRM,CBNT,CBNT,NfnI,NvnSOut,NvnSIn,1.0,ChiS_fI,FACET->qhat_uhatOutIn[dim]); // free
				for (n = 0; n < NfnI; n++) {
				for (j = 0; j < NvnSOut; j++) {
					dnqNumduhatOut_fI[n*NvnSOut+j] += (q_avg[NfnI*dim+n]+q_jump[NfnI*dim+n])*(q_uhatF[n*NvnSOut+j]);
				}}
				free(q_uhatF);
			}

			// OutOut contribution
			jacobian_flux_coef(NfnI,1,n_fI,h,FACET->P,gradu_avg,u_jump,q_avg,q_jump,d,"IP",'R');

			ChiS_fI_std = OPSOut->ChiS_fI[VfOut];
			ChiS_fI = malloc(NfnI*NvnSOut * sizeof *ChiS_fI); // free

			// Reordering
			for (i = 0; i < NfnI; i++) {
			for (j = 0; j < NvnSOut; j++) {
				ChiS_fI[i*NvnSOut+j] = ChiS_fI_std[nOrdOutIn[i]*NvnSOut+j];
			}}

			for (dim = 0; dim < d; dim++) {
				q_uhatV = mm_Alloc_d(CBRM,CBNT,CBNT,NfnI,NvnSOut,NvnSOut,1.0,ChiS_fI,VOut->qhat_uhat[dim]);        // free
				q_uhatF = mm_Alloc_d(CBRM,CBNT,CBNT,NfnI,NvnSOut,NvnSOut,1.0,ChiS_fI,FACET->qhat_uhatOutOut[dim]); // free
				for (n = 0; n < NfnI; n++) {
				for (j = 0; j < NvnSOut; j++) {
					dnqNumduhatOut_fI[n*NvnSOut+j]
						+= (q_avg[NfnI*dim+n]+q_jump[NfnI*dim+n])*(q_uhatV[n*NvnSOut+j]+q_uhatF[n*NvnSOut+j]);
				}}
				free(q_uhatV);
				free(q_uhatF);
			}
			free(ChiS_fI);
		} else {
			free(ChiS_fI);

			// Include BC information in dnqNumduhatIn_fI

			// Note: This approach is only possible because duOutduIn is a constant here. Otherwise, it would be
			//       better to linearize the flux wrt the solution itself (not the coefficients), add the duOutduIn
			//       term, and then add the final contribution for the linearization (du/duhat).
			duOutduIn = malloc(NfnI * sizeof *duOutduIn); // free

			if (BC % BC_STEP_SC == BC_DIRICHLET)
				jacobian_boundary_Dirichlet(NfnI,1,duOutduIn);
			else if (BC % BC_STEP_SC == BC_NEUMANN)
				printf("Error: Add support.\n"), EXIT_MSG;
			else
				printf("Error: Unsupported BC.\n"), EXIT_MSG;

			for (n = 0; n < NfnI; n++) {
			for (j = 0; j < NvnSIn; j++) {
				dnqNumduhatIn_fI[n] += dnqNumduhatOut[n]*duOutduIn[n];
			}}

			free(duOutduIn);
		}

		// Multiply by area element
		for (n = 0; n < NfnI; n++) {
			nqNum_fI[n] *= detJF_fI[n];
			for (j = 0; j < NvnSIn; j++)
				dnqNumduhatIn_fI[n*NvnSIn+j] *= detJF_fI[n];
			for (j = 0; j < NvnSOut; j++)
				dnqNumduhatOut_fI[n*NvnSOut+j] *= detJF_fI[n];
		}

		// Finalize FACET RHS and LHS terms

		// Interior FACET

		mm_CTN_d(NvnSIn,1,NfnI,OPSIn->I_Weak_FF[VfIn],nqNum_fI,RHSIn);
		mm_d(CBRM,CBNT,CBNT,NvnSIn,NvnSIn,NfnI,1.0,1.0,OPSIn->I_Weak_FF[VfIn],dnqNumduhatIn_fI,&LHSInIn);

		// Exterior FACET
		if (!Boundary) {
			// RHS

			// Use -ve normal for opposite VOLUME
			for (n = 0; n < NfnI; n++)
				nqNum_fI[n] *= -1.0;

			// Rearrange nqNum to match node ordering from VOut
			array_rearrange(NfnI,1,nOrdInOut,nqNum_fI);

			mm_CTN_d(NvnSOut,1,NfnI,OPSOut->I_Weak_FF[VfOut],nqNum_fI,RHSOut);

			// LHS

			// OutIn
			mm_d(CBRM,CBNT,CBNT,NvnSIn,NvnSOut,NfnI,1.0,1.0,OPSIn->I_Weak_FF[VfIn],dnqNumduhatOut_fI,&LHSOutIn);

			// Use -ve normal for opposite VOLUME
			for (n = 0; n < NfnI; n++) {
				for (j = 0; j < NvnSIn; j++)
					dnqNumduhatIn_fI[n*NvnSIn+j] *= -1.0;
				for (j = 0; j < NvnSOut; j++)
					dnqNumduhatOut_fI[n*NvnSOut+j] *= -1.0;
			}

			// Rearrange to match ordering from VOut
			array_rearrange(NfnI,NvnSIn, nOrdInOut,dnqNumduhatIn_fI);
			array_rearrange(NfnI,NvnSOut,nOrdInOut,dnqNumduhatOut_fI);

			// InOut
			mm_d(CBRM,CBNT,CBNT,NvnSOut,NvnSIn,NfnI,1.0,1.0,OPSOut->I_Weak_FF[VfOut],dnqNumduhatIn_fI,&LHSInOut);

			// OutOut
			mm_d(CBRM,CBNT,CBNT,NvnSOut,NvnSOut,NfnI,1.0,1.0,OPSOut->I_Weak_FF[VfOut],dnqNumduhatOut_fI,&LHSOutOut);
		}
	}
}

void implicit_info_Poisson(void)
{
	compute_qhat_VOLUME();
	compute_qhat_FACET();
	finalize_qhat();

	compute_uhat_VOLUME();
	compute_uhat_FACET();
}

void solver_Poisson(void)
{
	// Initialize DB Parameters
	unsigned int OutputInterval = DB.OutputInterval,
	             Nvar           = DB.Nvar;

	// Standard datatypes
	double maxRHS;

	implicit_info_Poisson();
	maxRHS = finalize_LHS();

	// add source contribution to VOLUME->RHS

}
