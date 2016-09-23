// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "solver_poisson_c.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Parameters.h"
/*
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
 *		Compute RHS for Poisson solver using complex variables (for linearization testing).
 *
 *	Comments:
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

void compute_qhat_VOLUME_c(void)
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
		mm_dcc(CBCM,CBT,CBNT,NvnI,1,NvnS,1.0,0.0,ChiS_vI,VOLUME->uhat_c,u_vI);

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

		VOLUME->DxyzChiS = DxyzChiS;

		// Compute RHS term
		for (dim1 = 0; dim1 < d; dim1++) {
			Sxyz = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS,NvnI,NvnS,-1.0,MInv,DxyzChiS[dim1]); // keep

			// RHS
			qhat = malloc(NvnS*1 * sizeof *qhat); // keep
			mm_dcc(CBCM,CBT,CBNT,NvnI,1,NvnS,1.0,0.0,Sxyz,VOLUME->uhat_c,qhat);
			VOLUME->qhat_c[dim1] = qhat;
		}
		free(u_vI);
	}
	free(OPS);
}

static void boundary_Dirichlet(const unsigned int Nn, const unsigned int Nel, double *XYZ, double complex *uB,
                               double complex *grad_uB)
{
	unsigned int n, nMax;
	double       *uB_d;

	uB_d = malloc(Nn*Nel * sizeof *uB_d); // free
	compute_exact_solution(Nn*Nel,XYZ,*uB_d,NULL,0);

	for (n = 0, nMax = Nn*Nel; n < nMax; n++)
		uB[n] = uB_d[n];
}

static void trace_IP(const unsigned int Nn, const unsigned int Nel, double complex *uL, double complex *uR,
                     double complex *uNum)
{
	unsigned int n, NnTotal;

	NnTotal = Nn*Nel;

	for (n = 0; n < NnTotal; n++)
		uNum[n] = 0.5*(uL[n]+uR[n]);
}

static void jacobian_trace_IP(const unsigned int Nn, const unsigned int Nel, double complex *duNumdu, const char side)
{
	unsigned int n, NnTotal;

	if (side != 'L' || side != 'R')
		printf("Error: Invalid side.\n"), EXIT_MSG;

	NnTotal = Nn*Nel;

	for (n = 0; n < NnTotal; n++)
		duNumdu[n] = 0.5;
}

void compute_qhat_FACET_c(void)
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
		uIn_fI = malloc(nfnI * sizeof *uIn_fI); // free
		mm_dcc(CBCM,CBT,CBNT,NfnI,1,NvnSIn,1.0,0.0,OPSIn->ChiS_fI[VfIn],VIn->uhat_c,uIn_fI);

		// Compute_uOut_fI (Taking BCs into account if applicable)
		uOut_fIIn = malloc(nfnI * sizeof *uOut_fIIn); // free
		if (!Boundary) {
			uOut_fI = malloc(nfnI * sizeof *uOut_fI); // free
			mm_dcc(CBCM,CBT,CBNT,NfnI,1,NvnSOut,1.0,0.0,OPSOut->ChiS_fI[VfOut],VOut->uhat_c,uOut_fI);

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

		// Compute numerical trace
		switch (PoissonTraceType) {
		case IP:
			trace_IP(NfnI,1,uIn_fI,uOut_fIIn,uNum_fI);
			break;
		default:
			printf("Error: Unsupported PoissonTraceType.\n"), EXIT_MSG;
			break;
		}
		free(uIn_fI);
		free(uOut_fIIn);

		// Multiply uNum and its Jacobian by the normal vector and area element
		nuNum_fI       = malloc(NfnI*d * sizeof *nuNum_fI); // tbd
		for (dim = 0; dim < d; dim++) {
		for (n = 0; n < NfnI; n++) {
			nJ = n_fI[n*d+dim]*detJF_fI[n];
			nuNum_fI[dim*NfnI+n] = uNum_fI[n]*nJ;
		}}


		// Compute FACET RHS terms
		RowTracker = malloc(NfnI * sizeof *RowTracker); // free

		// Interior VOLUME

		I_FF     = OPSIn->I_Weak_FF[VfIn];
		MInvI_FF = mm_Alloc_d(CBRM,CBNT,CBNT,NvnSIn,NfnI,NvnSIn,1.0,VIn->MInv,I_FF); // free

		nuNum         = malloc(NfnI * sizeof *nuNum); // free
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
				nuNum_fI[i]      *= -1.0;
			}

			// Rearrange numerical trace to match node ordering from opposite VOLUME
			for (i = 0; i < NfnI; i++)
				RowTracker[i] = i;

			for (RowInd = 0; RowInd < NfnI; RowInd++) {
				ReOrder = nOrdInOut[RowInd];
				for (RowSub = ReOrder; RowTracker[RowSub] != ReOrder; RowSub = RowTracker[RowSub])
					;

				if (RowInd != RowSub) {
					array_swap_d(&nuNum[RowInd],&nuNum[RowSub],d,NfnI);
					array_swap_ui(&RowTracker[RowInd],&RowTracker[RowSub],1,1);
				}
			}

			I_FF     = OPSOut->I_Weak_FF[VfOut];
			MInvI_FF = mm_Alloc_d(CBRM,CBNT,CBNT,NvnSOut,NfnI,NvnSOut,1.0,VOut->MInv,I_FF); // free

			ChiS_fIOutIn = malloc(NvnSOut*NfnI * sizeof *ChiS_fIOutIn); // free

			ChiS_fI = OPSIn->ChiS_fI[VfIn];
			for (i = 0; i < NfnI; i++) {
			for (j = 0; j < NvnSIn; j++) {
				ChiS_fIInOut[i*NvnSIn+j] = ChiS_fI[nOrdInOut[i]*NvnSIn+j];
			}}

			for (dim = 0; dim < d; dim++) {
				// RHSOut
				qhatOut = malloc(NvnSOut*1 * sizeof *qhatOut); // keep
				mm_dcc(CBCM,CBT,CBNT,NvnSOut,1,NfnI,1.0,0.0,MInvI_FF,&nuNum_fI[NfnI*dim],qhatOut);
				FACET->qhatOut_c[dim] = qhatIn;
			}
			free(ChiS_fIInOut);
		}
		free(MInvI_FF);

		free(nuNum_fI);
		free(RowTracker);
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

		free(FACET->qhatIn_c);
		free(FACET->qhatOut_c);
	}
}

static void compute_source_c(const unsigned int Nn, const double *XYZ, double complex *source)
{
	double *source_d;

	source_d = malloc(Nn * sizeof *source_d); // free
	compute_source(Nn,XYZ,source_d);

	for (n = 0; n < Nn; n++)
		source[n] = source_d[n];
	free(source_d);
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
			mm_dcc(CBCM,CBT,CBNT,NvnI,1,NvnS,1.0,0.0,ChiS_vI,VOLUME->qhat_c[dim1],&q_vI[dim1*NvnI]);

		MInv = VOLUME->MInv;
		C_vI = VOLUME->C_vI;

		// Compute RHS terms
		DxyzChiS = VOLUME->DxyzChiS;

		// RHS
		if (VOLUME->RHS_c)
			free(VOLUME->RHS_c);
		RHS = calloc(NvnS , sizeof *RHS); // keep
		VOLUME->RHS_c = RHS;

		Dq = malloc(NvnS * sizeof *Dq); // free
		for (dim1 = 0; dim1 < d; dim1++) {
			mm_dcc(CBCM,CBT,CBNT,NvnS,1,NvnI,1.0,0.0,DxyzChiS[dim1],VOLUME->qhat_c[dim1],Dq);

			for (n = 0; n < NvnS; n++)
				RHS[n] -= Dq[n];
		}
		free(Dq);
		free(q_vI);

		// RHS (Source)
// Can replace this with I_Weak_VV (ToBeDeleted)
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

// Can modify this to use mm_d and not need new memory allocation (ToBeDeleted)
		RHSs = malloc(NvnS * sizeof *RHSs); // free
		mm_dcc(CBCM,CBT,CBNT,NvnS,1,NvnI,1.0,0.0,ChiSTwdetJV_vI,f_vI,RHSs);
		free(ChiSwdetJV_vI);

		for (n = 0; n < NvnS; n++)
			RHS[n] -= RHSs[n];

		free(RHSs);
	}
	free(OPS);
}

static void flux_IP(const unsigned int Nn, const unsigned int Nel, double complex *uIn, double complex *uOut, double complex *grad_uIn,
                    double complex *grad_uOut, double complex *qIn, double complex *qOut, double *h, const unsigned int P, double complex *nqNum,
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

		// Compute uIn_fI and gradu_In_fI
		detJV_fI = FACET->detJV_fI;

		uIn_fI         = malloc(NfnI   * sizeof *uIn_fI);         // tbd
		grad_uIn_fI    = calloc(NfnI*d , sizeof *graduIn_fI);     // tbd

		mm_dcc(CBCM,CBT,CBNT,NfnI,1,NvnSIn,1.0,0.0,OPSIn->ChiS_fI[VfIn],VIn->uhat_c,uIn_fI);

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
			mm_dcc(CBCM,CBT,CBNT,NfnI,1,NvnSIn,1.0,0.0,GradxyzIn[dim],VIn->uhat_c,&grad_uIn_fI[NfnI*dim]);


		// Compute_uOut_fI (Taking BCs into account if applicable)
		uOut_fIIn         = malloc(NfnI   * sizeof *uOut_fIIn);       // tbd
		grad_uOut_fIIn    = calloc(NfnI*d , sizeof *graduOut_fI);     // tbd

		if (!Boundary) {
			uOut_fI         = malloc(NfnI   * sizeof *uOut_fI);         // free

			mm_CTN_d(NfnI,1,NvnSOut,OPSOut->ChiS_fI[VfOut],VOut->uhat,uOut_fI);

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
				mm_ddc(CBCM,CBT,CBNT,NfnI,1,NvnSOut,1.0,0.0,GradxyzOut[dim],VOut->uhat_c,&grad_uOut_fIIn[NfnI*dim]);
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

		// Compute numerical flux
		nqNum_fI = malloc(NfnI * sizeof *nqNum_fI); // tbd

		switch (PoissonFluxType) {
		case IP:
			flux_IP(NfnI,1,uIn_fI,uOut_fIIn,graduIn_fI,graduOut_fIIn,NULL,NULL,h,FACET->P,nqNum_fI,n_fI,d);
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

		// Interior FACET
		mm_ddc(CBCM,CBT,CBNT,NvnSIn,1,NfnI,1.0,0.0,OPSIn->I_Weak_FF[VfIn],nqNum_fI,RHSIn);

		// Exterior FACET
		if (!Boundary) {
			// RHS

			// Use -ve normal for opposite VOLUME
			for (n = 0; n < NfnI; n++)
				nqNum_fI[n] *= -1.0;

			// Rearrange nqNum to match node ordering from VOut
			array_rearrange(NfnI,1,nOrdInOut,nqNum_fI);

			mm_dcc(CBCM,CBT,CBNT,NvnSOut,1,NfnI,1.0,0.0,OPSOut->I_Weak_FF[VfOut],nqNum_fI,RHSOut);
		}
	}
}
