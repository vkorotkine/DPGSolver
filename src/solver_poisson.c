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

#include "element_functions.h"
#include "update_VOLUMEs.h"
#include "matrix_functions.h"
#include "exact_solutions.h"
#include "setup_geom_factors.h"
#include "array_swap.h"
#include "array_free.h"
#include "finalize_LHS.h"
#include "solver_implicit.h"
#include "output_to_paraview.h"

#include "array_print.h" // ToBeDeleted

/*
 *	Purpose:
 *		Perform the implicit solve for the Poisson equation.
 *
 *	Comments:
 *		Many of the RHS terms computed are 0. They are included as they are used to check the linearization. Further,
 *		the computational cost is dominated by the global system solve making this additional cost negligible.
 *		When finished with the implicit functions, include these redundant terms only in the explicit functions.
 *		(ToBeDeleted)
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

static void compute_qhat_VOLUME(void)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int i, j, dim1, dim2, IndD, IndC, NvnI, NvnS;
	double       *ChiS_vI, *MInv, **D, *C_vI, *Dxyz, **DxyzChiS, *Sxyz;

	struct S_OPERATORS *OPS;
	struct S_VOLUME    *VOLUME;

	OPS = malloc(sizeof *OPS); // free

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		compute_inverse_mass(VOLUME);

		init_ops(OPS,VOLUME);

		NvnI = OPS->NvnI;
		NvnS = OPS->NvnS;

		ChiS_vI = OPS->ChiS_vI;

		// Construct physical derivative operator matrices
		MInv = VOLUME->MInv;
		C_vI = VOLUME->C_vI;

		D = OPS->D_Weak;

// Note: Storage of DxyzChiS is not necessary if cofactor terms are included in uhat^ref and qhat^ref as in the paper.
// ToBeDeleted
//		DxyzChiS = malloc(d * sizeof *Dxyz); // keep
		DxyzChiS = VOLUME->DxyzChiS;
		for (dim1 = 0; dim1 < d; dim1++) {
			Dxyz = calloc(NvnS*NvnI , sizeof *Dxyz); // free
			for (dim2 = 0; dim2 < d; dim2++) {
				IndD = 0;
				IndC = (dim1+dim2*d)*NvnI;
				for (i = 0; i < NvnS; i++) {
					for (j = 0; j < NvnI; j++) {
						Dxyz[IndD+j] += D[dim2][IndD+j]*C_vI[IndC+j];
					}
					IndD += NvnI;
				}
			}
			DxyzChiS[dim1] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS,NvnS,NvnI,1.0,Dxyz,ChiS_vI); // keep
			free(Dxyz);
		}

		// Compute RHS and LHS terms
		for (dim1 = 0; dim1 < d; dim1++) {
			Sxyz = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS,NvnS,NvnS,-1.0,MInv,DxyzChiS[dim1]); // keep

			// RHS
			VOLUME->qhat[dim1] = mm_Alloc_d(CBCM,CBT,CBNT,NvnS,1,NvnS,1.0,Sxyz,VOLUME->uhat); // keep

			// LHS
			VOLUME->qhat_uhat[dim1] = Sxyz;
		}
	}
	free(OPS);
}

void boundary_Dirichlet(const unsigned int Nn, const unsigned int Nel, double *XYZ, double *uL, double *uR)
{
	/*
	 *	Comments:
	 *		uR = -uL + 2*uB.
	 */

	unsigned int n;

	if (Nel != 1)
		printf("Error: Vectorization unsupported.\n"), EXIT_MSG;

	compute_exact_solution(Nn*Nel,XYZ,uR,0);

	for (n = 0; n < Nn; n++) {
		uR[n] *= 2.0;
		uR[n] -= uL[n];
	}
}

static void jacobian_boundary_Dirichlet(const unsigned int Nn, const unsigned int Nel, double *duRduL)
{
	unsigned int n;

	if (Nel != 1)
		printf("Error: Vectorization unsupported.\n"), EXIT_MSG;

	for (n = 0; n < Nn; n++)
		duRduL[n] = -1.0;
}

void trace_coef(const unsigned int Nn, const unsigned int Nel, const double *nL, double *u_avg, double *u_jump,
                const unsigned int d, const char *trace_type)
{
	/*
	 *	Comments:
	 *		(T)race coefs for the various options are given as follows
	 *			uNum = u_avg * (uL+uR) + u_jump * (uL-uR)
	 *
	 *		Likely include normal here in future (ToBeDeleted).
	 *		Mistake in u_jump here. Should have components in 'd' directions. If normal is included, then u_avg would
	 *		have components in the 'd' directions instead (ToBeDeleted).
	 */

	unsigned int n, dim;

	if (Nel != 1)
		printf("Error: Unsupported Nel.\n"), EXIT_MSG;

	if (strstr(trace_type,"IP")) {
		for (n = 0; n < Nn; n++) {
			u_avg[n]  = 0.5;
			u_jump[n] = 0.0;
			// ToBeDeleted
			for (dim = 0; dim < d; dim++)
				u_jump[n] += 0.0*nL[n*d+dim];
		}
	} else {
		printf("Error: Unsupported trace_type.\n"), EXIT_MSG;
	}
}

static void compute_qhat_FACET(void)
{
	// Initialize DB Parameters
	unsigned int d               = DB.d,
	             ViscousFluxType = DB.ViscousFluxType;

	// Standard datatypes
	unsigned int n, dim, i, j, iMax,
	             NfnI, NvnSIn, NvnSOut,
	             VfIn, VfOut, fIn, EclassIn, IndFType, BC, Boundary,
	             *nOrdOutIn, *nOrdInOut;
	double       nJ, *n_fI, *detJF_fI,
	             *uIn_fI, *uOut_fIIn, *uOut_fI, *uNum_fI, *duNumduIn_fI, *duNumduOut_fI, *duOutduIn,
	             *u_avg, *u_jump, *nuNum_fI, *dnuNumduIn_fI, *dnuNumduOut_fI,
	             *I_FF, *MInvI_FF, *ChiS_fI, *ChiS_fIOutIn, *ChiS_fIInOut,
	             *MInvIdnuNumdu;

	struct S_OPERATORS *OPSIn, *OPSOut;
	struct S_VOLUME    *VIn, *VOut;
	struct S_FACET     *FACET;

	OPSIn  = malloc(sizeof *OPSIn);  // free
	OPSOut = malloc(sizeof *OPSOut); // free

	for (FACET = DB.FACET; FACET; FACET = FACET->next) {
		setup_geom_factors_highorder(FACET);

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
				boundary_Dirichlet(NfnI,1,FACET->XYZ_fI,uIn_fI,uOut_fIIn);
			} else if (BC % BC_STEP_SC == BC_NEUMANN) {
				printf("Add support.\n"), EXIT_MSG;
			} else {
				printf("Error: Unsupported BC.\n"), EXIT_MSG;
			}
		}

		// Compute numerical trace and its Jacobians
		u_avg         = malloc(NfnI * sizeof *u_avg);         // free
		u_jump        = malloc(NfnI * sizeof *u_jump);        // free
		uNum_fI       = malloc(NfnI * sizeof *uNum_fI);       // free
		duNumduIn_fI  = malloc(NfnI * sizeof *duNumduIn_fI);  // free
		duNumduOut_fI = malloc(NfnI * sizeof *duNumduOut_fI); // free
		switch (ViscousFluxType) {
		case FLUX_IP:
			trace_coef(NfnI,1,n_fI,u_avg,u_jump,d,"IP");
			break;
		default:
			printf("Error: Unsupported ViscousFluxType.\n"), EXIT_MSG;
			break;
		}

		for (n = 0; n < NfnI; n++) {
			uNum_fI[n]       = u_avg[n]*(uIn_fI[n]+uOut_fIIn[n]) + u_jump[n]*(uIn_fI[n]-uOut_fIIn[n]);
			duNumduIn_fI[n]  = u_avg[n] + u_jump[n];
			duNumduOut_fI[n] = u_avg[n] - u_jump[n];
		}
		free(u_avg);
		free(u_jump);

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

			free(duOutduIn);
		}

		// Multiply uNum and its Jacobian by the normal vector and area element
		nuNum_fI       = malloc(NfnI*d * sizeof *nuNum_fI);       // free
		dnuNumduIn_fI  = malloc(NfnI*d * sizeof *dnuNumduIn_fI);  // free
		dnuNumduOut_fI = malloc(NfnI*d * sizeof *dnuNumduOut_fI); // free
		for (dim = 0; dim < d; dim++) {
		for (n = 0; n < NfnI; n++) {
			nJ = n_fI[n*d+dim]*detJF_fI[n];
			nuNum_fI[dim*NfnI+n]       = uNum_fI[n]*nJ;
			dnuNumduIn_fI[dim*NfnI+n]  = duNumduIn_fI[n]*nJ;
			dnuNumduOut_fI[dim*NfnI+n] = duNumduOut_fI[n]*nJ;
		}}
		free(uNum_fI);
		free(duNumduIn_fI);
		free(duNumduOut_fI);

		// Compute FACET RHS and LHS terms

		// Interior VOLUME

		// Note: There is a "-ve" sign embedded in the I_Weak_FF operator!
		I_FF     = OPSIn->I_Weak_FF[VfIn];
		MInvI_FF = mm_Alloc_d(CBRM,CBNT,CBNT,NvnSIn,NfnI,NvnSIn,-1.0,VIn->MInv,I_FF); // free

		MInvIdnuNumdu = malloc(NvnSIn*NfnI * sizeof *MInvIdnuNumdu); // free
		for (dim = 0; dim < d; dim++) {
			// RHSIn
			FACET->qhatIn[dim] = mm_Alloc_d(CBCM,CBT,CBNT,NvnSIn,1,NfnI,1.0,MInvI_FF,&nuNum_fI[NfnI*dim]); // keep

			// LHS (InIn)
			mm_diag_d(NvnSIn,NfnI,&dnuNumduIn_fI[NfnI*dim],MInvI_FF,MInvIdnuNumdu,'R','R');

			FACET->qhat_uhatInIn[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnSIn,NvnSIn,NfnI,1.0,MInvIdnuNumdu,OPSIn->ChiS_fI[VfIn]);
		}

		// Exterior VOLUME
		if (!Boundary) {
			ChiS_fIOutIn = malloc(NfnI*NvnSOut * sizeof *ChiS_fIOutIn); // free

			ChiS_fI = OPSOut->ChiS_fI[VfOut];
			for (i = 0; i < NfnI; i++) {
			for (j = 0; j < NvnSOut; j++) {
				ChiS_fIOutIn[i*NvnSOut+j] = ChiS_fI[nOrdOutIn[i]*NvnSOut+j];
			}}

			for (dim = 0; dim < d; dim++) {
				// LHS (OutIn)
				mm_diag_d(NvnSIn,NfnI,&dnuNumduOut_fI[NfnI*dim],MInvI_FF,MInvIdnuNumdu,'R','R');

				FACET->qhat_uhatOutIn[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnSIn,NvnSOut,NfnI,1.0,MInvIdnuNumdu,ChiS_fIOutIn);
			}
			free(ChiS_fIOutIn);
			free(MInvIdnuNumdu);

			free(MInvI_FF);

			// Use "-ve" normal for opposite VOLUME
			for (i = 0, iMax = NfnI*d; i < iMax; i++) {
				nuNum_fI[i]       *= -1.0;
				dnuNumduIn_fI[i]  *= -1.0;
				dnuNumduOut_fI[i] *= -1.0;
			}

			// Rearrange numerical trace and its Jacobians to match node ordering from opposite VOLUME
			array_rearrange_d(NfnI,d,nOrdInOut,'C',nuNum_fI);
			array_rearrange_d(NfnI,d,nOrdInOut,'C',dnuNumduIn_fI);
			array_rearrange_d(NfnI,d,nOrdInOut,'C',dnuNumduOut_fI);

			I_FF     = OPSOut->I_Weak_FF[VfOut];
			MInvI_FF = mm_Alloc_d(CBRM,CBNT,CBNT,NvnSOut,NfnI,NvnSOut,-1.0,VOut->MInv,I_FF); // free

			MInvIdnuNumdu = malloc(NvnSOut*NfnI * sizeof *MInvIdnuNumdu); // free

			ChiS_fIInOut = malloc(NfnI*NvnSOut * sizeof *ChiS_fIInOut); // free

			ChiS_fI = OPSIn->ChiS_fI[VfIn];
			for (i = 0; i < NfnI; i++) {
			for (j = 0; j < NvnSIn; j++) {
				ChiS_fIInOut[i*NvnSIn+j] = ChiS_fI[nOrdInOut[i]*NvnSIn+j];
			}}

			for (dim = 0; dim < d; dim++) {
				// RHSOut
				FACET->qhatOut[dim] = mm_Alloc_d(CBCM,CBT,CBNT,NvnSOut,1,NfnI,1.0,MInvI_FF,&nuNum_fI[NfnI*dim]); // keep

				// LHS (InOut)
				mm_diag_d(NvnSOut,NfnI,&dnuNumduIn_fI[NfnI*dim],MInvI_FF,MInvIdnuNumdu,'R','R');

				FACET->qhat_uhatInOut[dim] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnSOut,NvnSIn,NfnI,1.0,MInvIdnuNumdu,ChiS_fIInOut);

				// LHS (OutOut)
				mm_diag_d(NvnSOut,NfnI,&dnuNumduOut_fI[NfnI*dim],MInvI_FF,MInvIdnuNumdu,'R','R');

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
				*VqhatIn_ptr++ += *FqhatIn_ptr++;
		}

		if(!(FACET->Boundary)) {
			VOut    = FACET->VOut;
			NvnSOut = VOut->NvnS;

			for (dim = 0; dim < d; dim++) {
				VqhatOut_ptr = VOut->qhat[dim];
				FqhatOut_ptr = FACET->qhatOut[dim];

				for (iMax = NvnSOut; iMax--; )
					*VqhatOut_ptr++ += *FqhatOut_ptr++;
			}
		}

		for (dim = 0; dim < d; dim++) {
			free(FACET->qhatIn[dim]);
			FACET->qhatIn[dim] = NULL;
			free(FACET->qhatOut[dim]);
			FACET->qhatOut[dim] = NULL;
		}
	}
}

static void compute_uhat_VOLUME(void)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int dim1, NvnS;
	double       **DxyzChiS, *RHS, *LHS;

	struct S_OPERATORS *OPS;
	struct S_VOLUME    *VOLUME;

	OPS = malloc(sizeof *OPS); // free

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		init_ops(OPS,VOLUME);

		NvnS = OPS->NvnS;

		// Compute RHS and LHS terms
		DxyzChiS = VOLUME->DxyzChiS;

		// RHS
		if (VOLUME->RHS)
			free(VOLUME->RHS);
		RHS = calloc(NvnS , sizeof *RHS); // keep
		VOLUME->RHS = RHS;

		for (dim1 = 0; dim1 < d; dim1++)
			mm_d(CBCM,CBT,CBNT,NvnS,1,NvnS,-1.0,1.0,DxyzChiS[dim1],VOLUME->qhat[dim1],RHS);

		// LHS
		if (VOLUME->LHS)
			free(VOLUME->LHS);
		LHS = calloc(NvnS*NvnS , sizeof *LHS); // keep
		VOLUME->LHS = LHS;

		for (dim1 = 0; dim1 < d; dim1++)
			mm_d(CBRM,CBNT,CBNT,NvnS,NvnS,NvnS,-1.0,1.0,DxyzChiS[dim1],VOLUME->qhat_uhat[dim1],LHS);
	}
	free(OPS);
}

void jacobian_flux_coef(const unsigned int Nn, const unsigned int Nel, const double *nIn, const double *h,
                        const unsigned int P, double *gradu_avg, double *u_jump, const unsigned int d, char *flux_type,
                        char side)
{
	/*
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
			tau = 1e2*(P+1)*(P+1)/h[n];

			gradu_avg[Nn*dim+n] =  0.5*nIn[n*d+dim];
			u_jump[Nn*dim+n]    = -tau*nIn[n*d+dim]*nIn[n*d+dim];

			if (side == 'R')
				u_jump[Nn*dim+n] *= -1.0;
		}}
	} else {
		printf("Error: Unsupported flux_type.\n"), EXIT_MSG;
	}
}

static void compute_uhat_FACET()
{
	// Initialize DB Parameters
	unsigned int d               = DB.d,
	             ViscousFluxType = DB.ViscousFluxType;

	// Standard datatypes
	unsigned int i, j, n, dim, dim1, dim2, IndC,
	             NvnSIn, NvnSOut, NfnI, NvnCIn,
	             BC, Boundary, VfIn, VfOut, fIn, EclassIn, IndFType,
	             *nOrdOutIn, *nOrdInOut;
	double       **GradChiS_fI, **GradxyzIn, **GradxyzOut, *ChiS_fI, *ChiS_fI_std,
	             *LHSInIn, *LHSOutIn, *LHSInOut, *LHSOutOut,
	             *detJVIn_fI, *detJVOut_fI, *C_fI, *h, *n_fI, *detJF_fI, *C_vC,
	             *uIn_fI, *grad_uIn_fI, *uOut_fIIn, *grad_uOut_fIIn, *uOut_fI,
	             *nqNum_fI, *dnqNumduhatIn_fI, *dnqNumduhatOut_fI, *duOutduIn,
	             *gradu_avg, *u_jump,
	             *RHSIn, *RHSOut;

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

		// Add VOLUME contributions to RHS and LHS
		RHSIn  = calloc(NvnSIn  , sizeof *RHSIn);  // keep (requires external free)
		RHSOut = calloc(NvnSOut , sizeof *RHSOut); // keep (requires external free)

		FACET->RHSIn  = RHSIn;
		FACET->RHSOut = RHSOut;

		LHSInIn   = calloc(NvnSIn*NvnSIn   , sizeof *LHSInIn);   // keep
		LHSOutIn  = calloc(NvnSIn*NvnSOut  , sizeof *LHSOutIn);  // keep
		LHSInOut  = calloc(NvnSOut*NvnSIn  , sizeof *LHSInOut);  // keep
		LHSOutOut = calloc(NvnSOut*NvnSOut , sizeof *LHSOutOut); // keep

		FACET->LHSInIn   = LHSInIn;
		FACET->LHSOutIn  = LHSOutIn;
		FACET->LHSInOut  = LHSInOut;
		FACET->LHSOutOut = LHSOutOut;

		for (dim = 0; dim < d; dim++) {
			// RHS
			mm_d(CBCM,CBT,CBNT,NvnSIn, 1,NvnSIn, -1.0,1.0,VIn->DxyzChiS[dim], FACET->qhatIn[dim], RHSIn);
			if (!Boundary)
				mm_d(CBCM,CBT,CBNT,NvnSOut,1,NvnSOut,-1.0,1.0,VOut->DxyzChiS[dim],FACET->qhatOut[dim],RHSOut);

			// LHS
			mm_d(CBRM,CBNT,CBNT,NvnSIn, NvnSIn, NvnSIn, -1.0,1.0,VIn->DxyzChiS[dim], FACET->qhat_uhatInIn[dim],  LHSInIn);
			if (!Boundary) {
				mm_d(CBRM,CBNT,CBNT,NvnSIn, NvnSOut,NvnSIn, -1.0,1.0,VIn->DxyzChiS[dim], FACET->qhat_uhatOutIn[dim], LHSOutIn);
				mm_d(CBRM,CBNT,CBNT,NvnSOut,NvnSIn, NvnSOut,-1.0,1.0,VOut->DxyzChiS[dim],FACET->qhat_uhatInOut[dim], LHSInOut);
				mm_d(CBRM,CBNT,CBNT,NvnSOut,NvnSOut,NvnSOut,-1.0,1.0,VOut->DxyzChiS[dim],FACET->qhat_uhatOutOut[dim],LHSOutOut);
			}
		}

		// Compute uIn_fI and gradu_In_fI

		uIn_fI      = malloc(NfnI   * sizeof *uIn_fI);      // free
		grad_uIn_fI = calloc(NfnI*d , sizeof *grad_uIn_fI); // free

		mm_CTN_d(NfnI,1,NvnSIn,OPSIn->ChiS_fI[VfIn],VIn->uhat,uIn_fI);

		GradChiS_fI = OPSIn->GradChiS_fI[VfIn];

		GradxyzIn = malloc(d * sizeof *GradxyzIn); // free
		for (dim1 = 0; dim1 < d; dim1++) {
			GradxyzIn[dim1] = calloc(NfnI*NvnSIn , sizeof **GradxyzIn); // free

			for (dim2 = 0; dim2 < d; dim2++) {
				IndC = (dim1+dim2*d)*NfnI;
				for (n = 0; n < NfnI; n++) {
					for (j = 0; j < NvnSIn; j++)
						GradxyzIn[dim1][n*NvnSIn+j] += GradChiS_fI[dim2][n*NvnSIn+j]*C_fI[IndC+n];
				}
			}
			for (n = 0; n < NfnI; n++) {
				for (j = 0; j < NvnSIn; j++)
					GradxyzIn[dim1][n*NvnSIn+j] /= detJVIn_fI[n];
			}
		}

		for (dim = 0; dim < d; dim++)
			mm_CTN_d(NfnI,1,NvnSIn,GradxyzIn[dim],VIn->uhat,&grad_uIn_fI[NfnI*dim]);

		mm_CTN_d(NfnI,d*d,OPSOut->NvnC,OPSOut->I_vC_fI[VfOut],VOut->C_vC,C_fI);
		GradChiS_fI = OPSOut->GradChiS_fI[VfOut];

		// Note: Rearrangement is embedded in the operator
		GradxyzOut = malloc(d * sizeof *GradxyzOut); // free
		for (dim1 = 0; dim1 < d; dim1++) {
			GradxyzOut[dim1] = calloc(NfnI*NvnSOut , sizeof **GradxyzOut); // free

			for (dim2 = 0; dim2 < d; dim2++) {
				IndC = (dim1+dim2*d)*NfnI;
				for (n = 0; n < NfnI; n++) {
					for (j = 0; j < NvnSOut; j++)
						GradxyzOut[dim1][n*NvnSOut+j] += GradChiS_fI[dim2][nOrdOutIn[n]*NvnSOut+j]*C_fI[IndC+nOrdOutIn[n]];
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

			mm_CTN_d(NfnI,1,NvnSOut,OPSOut->ChiS_fI[VfOut],VOut->uhat,uOut_fI);

			// Reorder uOut_fI to correspond to inner VOLUME ordering
			for (n = 0; n < NfnI; n++)
				uOut_fIIn[n] = uOut_fI[nOrdOutIn[n]];
			free(uOut_fI);

			for (dim = 0; dim < d; dim++)
				mm_CTN_d(NfnI,1,NvnSOut,GradxyzOut[dim],VOut->uhat,&grad_uOut_fIIn[NfnI*dim]);
		} else {
			if (BC % BC_STEP_SC == BC_DIRICHLET) {
				boundary_Dirichlet(NfnI,1,FACET->XYZ_fI,uIn_fI,uOut_fIIn);
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
		nqNum_fI = calloc(NfnI , sizeof *nqNum_fI); // free

		gradu_avg = malloc(NfnI*d * sizeof *gradu_avg); // free
		u_jump    = malloc(NfnI*d * sizeof *u_jump);    // free

		switch (ViscousFluxType) {
		case FLUX_IP:
			jacobian_flux_coef(NfnI,1,n_fI,h,FACET->P,gradu_avg,u_jump,d,"IP",'L');
			break;
		default:
			printf("Error: Unsupported ViscousFluxType.\n"), EXIT_MSG;
			break;
		}

		for (n = 0; n < NfnI; n++) {
		for (dim = 0; dim < d; dim++) {
			nqNum_fI[n] += gradu_avg[NfnI*dim+n] * (grad_uIn_fI[NfnI*dim+n] + grad_uOut_fIIn[NfnI*dim+n])
			            +  u_jump[NfnI*dim+n]    * (uIn_fI[n] - uOut_fIIn[n]);
		}}
		free(grad_uIn_fI);
		free(grad_uOut_fIIn);
		free(uIn_fI);
		free(uOut_fIIn);

// If using ViscousFluxType == "IP", what is done below is quite inefficient as there is no dependence on q. It is very
// general however and supports many different FluxTypes.
		dnqNumduhatIn_fI    = calloc(NfnI*NvnSIn  , sizeof *dnqNumduhatIn_fI);    // free
		dnqNumduhatOut_fI   = calloc(NfnI*NvnSOut , sizeof *dnqNumduhatOut_fI);   // free


		// InIn contributions
		jacobian_flux_coef(NfnI,1,n_fI,h,FACET->P,gradu_avg,u_jump,d,"IP",'L');
		ChiS_fI = OPSIn->ChiS_fI[VfIn];

		for (dim = 0; dim < d; dim++) {
			for (n = 0; n < NfnI; n++) {
			for (j = 0; j < NvnSIn; j++) {
				dnqNumduhatIn_fI[n*NvnSIn+j] // row-major
					+= gradu_avg[NfnI*dim+n]*GradxyzIn[dim][n*NvnSIn+j] + u_jump[NfnI*dim+n]*ChiS_fI[n*NvnSIn+j];
			}}
		}
		array_free2_d(d,GradxyzIn);

		// Note: Both dnqNumduhatIn and dnqNumduhatOut are ordered corresponding to VIn.
		if (!Boundary) {
			ChiS_fI_std = OPSOut->ChiS_fI[VfOut];
			ChiS_fI = malloc(NfnI*NvnSOut * sizeof *ChiS_fI); // free

			// Reordering
			for (i = 0; i < NfnI; i++) {
			for (j = 0; j < NvnSOut; j++) {
				ChiS_fI[i*NvnSOut+j] = ChiS_fI_std[nOrdOutIn[i]*NvnSOut+j];
			}}

			jacobian_flux_coef(NfnI,1,n_fI,h,FACET->P,gradu_avg,u_jump,d,"IP",'R');

			for (dim = 0; dim < d; dim++) {
				for (n = 0; n < NfnI; n++) {
				for (j = 0; j < NvnSOut; j++) {
					dnqNumduhatOut_fI[n*NvnSOut+j]
						+= gradu_avg[NfnI*dim+n]*GradxyzOut[dim][n*NvnSOut+j] + u_jump[NfnI*dim+n]*ChiS_fI[n*NvnSOut+j];
				}}
			}
			free(ChiS_fI);
		} else {
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

			// OutOut contribution (u)
			jacobian_flux_coef(NfnI,1,n_fI,h,FACET->P,gradu_avg,u_jump,d,"IP",'R');
			ChiS_fI = OPSIn->ChiS_fI[VfIn];

			// Note: dgraduOutduIn == -duOutduIn.
			for (dim = 0; dim < d; dim++) {
				for (n = 0; n < NfnI; n++) {
				for (j = 0; j < NvnSIn; j++) {
					dnqNumduhatIn_fI[n*NvnSIn+j] += (- duOutduIn[n]*gradu_avg[NfnI*dim+n]*GradxyzOut[dim][n*NvnSIn+j]
				                                     + duOutduIn[n]*u_jump[NfnI*dim+n]*ChiS_fI[n*NvnSIn+j]);
				}}
			}
			free(duOutduIn);
		}
		free(gradu_avg);
		free(u_jump);
		free(h);

		array_free2_d(d,GradxyzOut);

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

		mm_d(CBCM,CBT,CBNT,NvnSIn,1,NfnI,-1.0,1.0,OPSIn->I_Weak_FF[VfIn],nqNum_fI,RHSIn);
		mm_d(CBRM,CBNT,CBNT,NvnSIn,NvnSIn,NfnI,-1.0,1.0,OPSIn->I_Weak_FF[VfIn],dnqNumduhatIn_fI,LHSInIn);

		// Exterior FACET
		if (!Boundary) {
			// RHS

			// Use -ve normal for opposite VOLUME
			for (n = 0; n < NfnI; n++)
				nqNum_fI[n] *= -1.0;

			// Rearrange nqNum to match node ordering from VOut
			array_rearrange_d(NfnI,1,nOrdInOut,'C',nqNum_fI);

			mm_d(CBCM,CBT,CBNT,NvnSOut,1,NfnI,-1.0,1.0,OPSOut->I_Weak_FF[VfOut],nqNum_fI,RHSOut);

			// LHS

			// OutIn
			mm_d(CBRM,CBNT,CBNT,NvnSIn,NvnSOut,NfnI,-1.0,1.0,OPSIn->I_Weak_FF[VfIn],dnqNumduhatOut_fI,LHSOutIn);

			// Use -ve normal for opposite VOLUME
			for (n = 0; n < NfnI; n++) {
				for (j = 0; j < NvnSIn; j++)
					dnqNumduhatIn_fI[n*NvnSIn+j] *= -1.0;
				for (j = 0; j < NvnSOut; j++)
					dnqNumduhatOut_fI[n*NvnSOut+j] *= -1.0;
			}

			// Rearrange to match ordering from VOut
			array_rearrange_d(NfnI,NvnSIn,nOrdInOut,'R',dnqNumduhatIn_fI);
			array_rearrange_d(NfnI,NvnSOut,nOrdInOut,'R',dnqNumduhatOut_fI);

			// InOut
			mm_d(CBRM,CBNT,CBNT,NvnSOut,NvnSIn,NfnI,-1.0,1.0,OPSOut->I_Weak_FF[VfOut],dnqNumduhatIn_fI,LHSInOut);

			// OutOut
			mm_d(CBRM,CBNT,CBNT,NvnSOut,NvnSOut,NfnI,-1.0,1.0,OPSOut->I_Weak_FF[VfOut],dnqNumduhatOut_fI,LHSOutOut);
		}
		free(nqNum_fI);
		free(dnqNumduhatIn_fI);
		free(dnqNumduhatOut_fI);
	}
	free(OPSIn);
	free(OPSOut);
}

void implicit_info_Poisson(void)
{
	compute_qhat_VOLUME();
	compute_qhat_FACET();

	compute_uhat_VOLUME();
	compute_uhat_FACET();
}

void solver_Poisson(void)
{
	// Initialize DB Parameters
	unsigned int Nvar = DB.Nvar;

	// Standard datatypes
	char         *string, *fNameOut;
	unsigned int i, iMax, IndA, NvnS;
	int          iteration_ksp;
	double       maxRHS, *uhat, *duhat;

	struct S_VOLUME *VOLUME;

	fNameOut = malloc(STRLEN_MAX * sizeof *fNameOut); // free
	string   = malloc(STRLEN_MIN * sizeof *string);   // free

	// Build the RHS and LHS terms
	Mat                A = NULL;
	Vec                b = NULL, x = NULL;
	KSP                ksp;
	KSPConvergedReason reason;

	PetscInt  *ix;
	PetscReal emax, emin;

	printf("RL"); implicit_info_Poisson();
	printf("F "); maxRHS = finalize_LHS(&A,&b,&x,0);

//	MatView(A,PETSC_VIEWER_STDOUT_SELF);
//	VecView(b,PETSC_VIEWER_STDOUT_SELF);
//	EXIT_MSG;

	// Solve linear system
	printf("S");
	KSPCreate(MPI_COMM_WORLD,&ksp);
	setup_KSP(A,ksp);

	printf("S ");
	KSPSolve(ksp,b,x);
	KSPGetConvergedReason(ksp,&reason);
	KSPGetIterationNumber(ksp,&iteration_ksp);
	KSPComputeExtremeSingularValues(ksp,&emax,&emin);

//	KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);

	// Update uhat
	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		IndA = VOLUME->IndA;
		NvnS = VOLUME->NvnS;
		uhat = VOLUME->uhat;

		iMax = NvnS*Nvar;
		ix    = malloc(iMax * sizeof *ix);    // free
		duhat = malloc(iMax * sizeof *duhat); // free
		for (i = 0; i < iMax; i++)
			ix[i] = IndA+i;

		VecGetValues(x,iMax,ix,duhat);
		free(ix);

		for (i = 0; i < iMax; i++)
			(*uhat++) += duhat[i];
		free(duhat);
	}

	KSPDestroy(&ksp);
	finalize_ksp(&A,&b,&x,2);

	// Update qhat
	compute_qhat_VOLUME();
	compute_qhat_FACET();
	finalize_qhat();

	// Output to paraview
	strcpy(fNameOut,"SolFinal_");
	sprintf(string,"%dD_",DB.d);   strcat(fNameOut,string);
	                               strcat(fNameOut,DB.MeshType);
	sprintf(string,"_ML%d",DB.ML); strcat(fNameOut,string);
	if (DB.Adapt == ADAPT_0)
		sprintf(string,"P%d_",DB.PGlobal), strcat(fNameOut,string);
	output_to_paraview(fNameOut);

	printf("KSP iterations (cond, reason): %5d (% .3e, %d), maxRHS (no MInv): % .3e\n",
	       iteration_ksp,emax/emin,reason,maxRHS);

	free(fNameOut);
	free(string);

	// Potentially adaptation option (ToBeDeleted)
}
