// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "solver_poisson.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"

/*
#include "Parameters.h"
#include "S_DB.h"
#include "S_VOLUME.h"

#include "adaptation.h"
#include "update_VOLUMEs.h"
#include "implicit_VOLUME_info.h"
#include "implicit_FACET_info.h"
#include "finalize_LHS.h"
#include "output_to_paraview.h"

#include "Macros.h"
#include "array_print.h"
*/

/*
 *	Purpose:
 *		Perform the implicit solve for the Poisson equation.
 *
 *	Comments:
 *		Many of the RHS terms computed are 0. They are included as they are used to check the linearization. Further,
 *		the computational cost is dominated by the global system solve making this additional cost negligible.
 *		When finished with the implicit functions, include these contributions, include these redundant terms only in
 *		the explicit functions. (ToBeDeleted)
 *
 *	Notation:
 *
 *	References:
 */

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

static void compute_qhat_VOLUME(void)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	double *MInv, **D, *C_vI, **Dxyz, *Sxyz;

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
		mm_CTN_d(NvnI,1,NvnS,ChiS_vI,VOLUME->What,u_vI);

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
			mm_Alloc_d(CBCM,CBT,CBNT,NvnS,1,NvnI,1.0,Sxyz,VOLUME->What,qhat); // keep
			VOLUME->qhat[dim1] = qhat;

			// LHS
			VOLUME->qhat_uhat[dim1] = Sxyz;
		}
		free(u_vI);
	}
	free(OPS);
}

static void boundary_Dirichlet(const unsigned int Nn, const unsigned int Nel, double *XYZ, double *uB)
{
	compute_exact_solution(Nn*Nel,XYZ,*uB,NULL,0);
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
	// Standard datatypes

	struct S_OPERATORS *OPSIn, *OPSOut;
	struct S_VOLUME    *VIn, *VOut;
	struct S_FACET     *FACET;

	for (FACET = DB.FACET; FACET; FACET = FACET->next) {
		VIn  = FACET->VIn;
		VfIn = FACET->VfIn;
		fIn  = FACET->fIn;

		EclassIn = VIn->Eclass;
		IndFType = get_IndFType(EclassIn,fIn);
		init_opsF(OPSIn,VIn,FACET,0);

		VOut  = FACET->VOut;
		VfOut = FACET->VfOut;

		init_opsF(OPSOut,VOut,FACET,0);

		BC       = FACET->BC;
		Boundary = FACET->Boundary;

		NfnI   = OPSIn[IndFType]->NfnI;
		NvnSIn = OPSIn[IndFType]->NvnS;

		// Compute uIn_fI
		uIn_fI = calloc(nfnI * sizeof *uIn_fI); // free

		// Compute_uOut_fI (Taking BCs into account if applicable)
		uOut_fIIn = calloc(nfnI * sizeof *uOut_fIIn); // free
		if (Boundary) {
			if (BC % BC_STEP_SC == BC_DIRICHLET) {
				boundary_Dirichlet(NfnI,1,FACET->XYZ_fI,uOut_fIIn);
			} else if (BC % BC_STEP_SC == BC_NEUMANN) {
				printf("Add support.\n"), EXIT_MSG;
			} else {
				printf("Error: Unsupported BC.\n"), EXIT_MSG;
			}
		}

		// Compute numerical trace and its Jacobian
		switch (PoissonTraceFluxType) {
		case IP:
			trace_IP(NfnI,1,uIn_fI,uOut_fIIn,uNum_fI);
			jacobian_trace_IP(NfnI,1,duNumduIn_fI,'L');
			jacobian_trace_IP(NfnI,1,duNumduOut_fI,'R');
			break;
		default:
			printf("Error: Unsupported PoissonTraceFluxType.\n"), EXIT_MSG;
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
				duNumduIn_fI[n] += duNumduOut[n]*duOutduIn[n];

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
			nuNum_fI[dim*NfnI+n]      = uNum_fI[n]*nJ;
			duNumduIn_fI[dim*NfnI+n]  = duNumduIn_fI[n]*nJ;
			duNumduOut_fI[dim*NfnI+n] = duNumduOut_fI[n]*nJ;
		}}

		// Compute FACET RHS and LHS terms
		RowTracker = malloc(NfnI * sizeof *RowTracker); // free

		// Interior VOLUME

		I_FF     = OPSIn->I_Weak_FF[VfIn];
		MInvI_FF = mm_Alloc_d(CBRM,CBNT,CBNT,NvnSIn,NfnI,NvnSIn,1.0,VIn->MInv,I_FF); // free

		nuNum         = malloc(NfnI * sizeof *nuNum); // free
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
			for (i = 0; i < NfnI; i++)
				RowTracker[i] = i;

			for (RowInd = 0; RowInd < NfnI; RowInd++) {
				ReOrder = nOrdInOut[RowInd];
				for (RowSub = ReOrder; RowTracker[RowSub] != ReOrder; RowSub = RowTracker[RowSub])
					;

				if (RowInd != RowSub) {
					array_swap_d(&nuNum[RowInd],&nuNum[RowSub],d,NfnI);
					array_swap_d(&dnuNumduIn[RowInd],&dnuNumduIn[RowSub],d,NfnI);
					array_swap_d(&dnuNumduOut[RowInd],&dnuNumduOut[RowSub],d,NfnI);
					array_swap_ui(&RowTracker[RowInd],&RowTracker[RowSub],1,1);
				}
			}

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

static void compute_source(const unsigned int Nn, double *XYZ, double *source)
{
	// Initialize DB Parameters
	char *TestCase = DB.TestCase;

	// Standard datatypes
	double *X, *Y;

	X = &XYZ[0*Nn];
	Y = &XYZ[1*Nn];

	if (strstr(TestCase,"Poisson")) {
		for (n = 0; n < Nn; n++)
			source[n] = -2.0*PI*PI*sin(PI*X[n])*sin(PI*Y[n]);
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

		// LHS (InIn) - VOLUME term
		DxyzChiS = VIn->DxyzChiS;

		LHSInIn   = calloc(NvnSIn*NvnSIn   , sizeof *LHSInIn); // keep

		FACET->LHSInIn   = LHSInIn;

		LHSd = malloc(NvnSIn*NvnSIn * sizeof *LHSd); // free

		for (dim = 0; dim < d; dim++) {
			LHSd = mm_Alloc_d(CBRM,CBNT,CBNT,NvnSIn,NvnSIn,NvnSIn,-1.0,DxyzChiS[dim],FACET->qhat_uhatInIn[dim]); // free

			for (i = 0, iMax = NvnSIn*NvnSIn; i < iMax; i++)
				LHSInIn[i] += LHSd;
		}

		// LHS (OutIn) - VOLUME term

		free(LHSd);

	}
}

void solver_poisson(void)
{
	// Initialize DB Parameters
	unsigned int OutputInterval = DB.OutputInterval,
	             Nvar           = DB.Nvar;

	// Standard datatypes

	compute_qhat_VOLUME();
	compute_qhat_FACET();
	finalize_qhat();

	compute_uhat_VOLUME();
	compute_uhat_FACET();

	// add source contribution to VOLUME->RHS

}
