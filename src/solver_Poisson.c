// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "solver_Poisson.h"

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

#include "array_print.h"

/*
 *	Purpose:
 *		Perform the implicit solve for the Poisson equation.
 *
 *	Comments:
 *		Many of the RHS terms computed are 0; they are included as they are used to check the linearization. Further,
 *		the computational cost is dominated by the global system solve making this additional cost negligible.
 *
 *		Originally, the Weak form was employed for both of the equations resulting from the transformation of the 2nd
 *		order equation to the system of 1st order equations. It was then determined that, for high enough order
 *		solutions (P > 5 for TRIs), the global system matrix failed to be symmetric (even when using straight elements).
 *		This was found to be a result of the impossibility of exactly integrating by parts the VOLUME term in the
 *		equation for qhat (even with very high cubature order (3P)). While it seems that the exact IBP should be
 *		possible in this case, the investigation was not pursued. Further, the strong form of the scheme very easily
 *		retains the symmetry required and is hence a more intuitive discretization.
 *
 *	Notation:
 *
 *	References:
 *		See alternate/solver_Poisson_weak.c for weak form implementation.
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
			if (DxyzChiS[dim1])
				free(DxyzChiS[dim1]);
			DxyzChiS[dim1] = mm_Alloc_d(CBRM,CBNT,CBNT,NvnS,NvnS,NvnI,1.0,Dxyz,ChiS_vI); // keep
			free(Dxyz);
		}

		// Compute RHS and LHS terms
		for (dim1 = 0; dim1 < d; dim1++) {
			Sxyz = mm_Alloc_d(CBRM,CBNT,CBT,NvnS,NvnS,NvnS,1.0,MInv,DxyzChiS[dim1]); // keep

			// RHS
			if (VOLUME->qhat[dim1])
				free(VOLUME->qhat[dim1]);
			VOLUME->qhat[dim1] = mm_Alloc_d(CBCM,CBT,CBNT,NvnS,1,NvnS,1.0,Sxyz,VOLUME->uhat); // keep

			// LHS
			if (VOLUME->qhat_uhat[dim1])
				free(VOLUME->qhat_uhat[dim1]);
			VOLUME->qhat_uhat[dim1] = Sxyz;
		}
	}
	free(OPS);
}

void project_to_sphere(const unsigned int Nn, double *XYZIn, double *XYZOut, const unsigned int curved)
{
	/*
	 *	Purpose:
	 *		Project coordinates to the sphere surface.
	 */

	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int n;
	double       r, r2, t, p, rIn, rOut, norm_rIn, norm_rOut, *XIn, *YIn, *ZIn, *XOut, *YOut, *ZOut;

	rIn  = 0.5;
	rOut = 1.0;

	XIn = &XYZIn[0*Nn];
	YIn = &XYZIn[1*Nn];
	ZIn = &XYZIn[(d-1)*Nn];

	XOut = &XYZOut[0*Nn];
	YOut = &XYZOut[1*Nn];
	ZOut = &XYZOut[(d-1)*Nn];

	norm_rIn  = 0.0;
	norm_rOut = 0.0;

	for (n = 0; n < Nn; n++) {
		r2 = XIn[n]*XIn[n]+YIn[n]*YIn[n];
		if (d == 3)
			r2 += ZIn[n]*ZIn[n];
		r          = sqrt(r2);
		norm_rIn  += sqrt((r-rIn)* (r-rIn));
		norm_rOut += sqrt((r-rOut)*(r-rOut));
	}
	norm_rIn  /= Nn;
	norm_rOut /= Nn;

	if (0)
		printf("%e %e\n",XOut[0],YOut[0]);

//	if (curved == 1) {
	if (1||curved == 1) {
		for (n = 0; n < Nn; n++) {
			XOut[n] = XIn[n];
			YOut[n] = YIn[n];
			if (d == 3)
				ZOut[n] = ZIn[n];
		}
	} else {
		// Project to correct radius
		if (norm_rIn < 2e-1) {
			r = rIn;
		} else if (norm_rOut < 2e-1) {
			r = rOut;
		} else {
			printf("% .3e % .3e\n",norm_rIn,norm_rOut);
			printf("Error: Found unsupported curved BC.\n"), EXIT_MSG;
		}

		for (n = 0; n < Nn; n++) {
			t = atan2(YIn[n],XIn[n]);
			XOut[n] = r*cos(t);
			YOut[n] = r*sin(t);

			if (d == 3) {
				p = acos(ZIn[n]/r);

				XOut[n] *= sin(p);
				YOut[n] *= sin(p);
				ZOut[n]  = r*cos(p);
			}
		}
	}
}

void boundary_Poisson(const unsigned int Nn, const unsigned int Nel, double *XYZ, double *uL, double *uR,
                      double *graduL, double *graduR, const unsigned int BC, const unsigned int curved)
{
	/*
	 *	Comments:
	 *		Available BC types:
	 *			(D)irichlet: uR = -uL + 2*uB; graduR =  graduL.
	 *			(N)eumann:   uR =  uL;        graduR = -graduL + 2*graduB.
	 *
	 *		The boundary solution must be computed on the exact geometry for the curving of the elements to become
	 *		necessary. Otherwise, the problem is simply being solved on the inexact domain.
	 */

	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int n;
	double       *XYZproj;

	if (Nel != 1)
		printf("Error: Vectorization unsupported.\n"), EXIT_MSG;

	// Modify XYZ coordinates for boundary condition evaluation
	XYZproj = malloc(Nn*d * sizeof *XYZproj); // free
	if (d > 1)
		project_to_sphere(Nn,XYZ,XYZproj,curved);
	else
		printf("Error: Unsupported.\n"), EXIT_MSG;

	if (BC == BC_DIRICHLET) {
		compute_exact_solution(Nn*Nel,XYZproj,uR,0);
		for (n = 0; n < Nn; n++) {
			uR[n] *= 2.0;
			uR[n] -= uL[n];
		}

		if (graduL) {
			for (n = 0; n < d*Nn; n++)
				graduR[n] = graduL[n];
		}
	} else if (BC == BC_NEUMANN) {
		for (n = 0; n < Nn; n++)
			uR[n] = uL[n];

		if (graduL) {
			compute_exact_gradient(Nn*Nel,XYZproj,graduR);
			for (n = 0; n < Nn*d; n++) {
				graduR[n] *= 2.0;
				graduR[n] -= graduL[n];
			}
		}
	} else {
		printf("Error: Unsupported BC_type.\n"), EXIT_MSG;
	}
	free(XYZproj);
}

static void jacobian_boundary_Poisson(const unsigned int Nn, const unsigned int Nel, double *duRduL,
                                      const unsigned int BC)
{
	unsigned int n;

	if (Nel != 1)
		printf("Error: Vectorization unsupported.\n"), EXIT_MSG;

	if (BC == BC_DIRICHLET) {
		for (n = 0; n < Nn; n++)
			duRduL[n] = -1.0;
	} else if (BC == BC_NEUMANN) {
		for (n = 0; n < Nn; n++)
			duRduL[n] = 1.0;
	} else {
		printf("Error: Unsupported BC_type.\n"), EXIT_MSG;
	}
}

static void compute_qhat_FACET(void)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int n, dim, i, j,
	             NfnI, NvnSIn, NvnSOut,
	             VfIn, VfOut, fIn, EclassIn, IndFType, BC, Boundary,
	             *nOrdOutIn, *nOrdInOut;
	double       *w_fI, *n_fI, *detJF_fI, *wnJ_fI, *mult_diag,
	             *uIn_fI, *uOut_fIIn, *uOut_fI, *uNum_fI, *duNumduIn_fI, *duNumduOut_fI, *duOutduIn,
	             *ChiS_fI, *ChiS_fIOutIn, *ChiS_fIInOut;

	struct S_OPERATORS *OPSIn, *OPSOut;
	struct S_VOLUME    *VIn, *VOut;
	struct S_FACET     *FACET;

	// silence
	uOut_fI = NULL;

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

		w_fI = OPSIn->w_fI;

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
		} else {
			boundary_Poisson(NfnI,1,FACET->XYZ_fI,uIn_fI,uOut_fIIn,NULL,NULL,BC % BC_STEP_SC,BC / BC_STEP_SC);
		}

		// Compute numerical trace and its Jacobians
		uNum_fI       = malloc(NfnI * sizeof *uNum_fI);       // free
		duNumduIn_fI  = malloc(NfnI * sizeof *duNumduIn_fI);  // free
		duNumduOut_fI = malloc(NfnI * sizeof *duNumduOut_fI); // free

	 	// Numerical trace: uNum = 0.5 * (uL+uR)
		for (n = 0; n < NfnI; n++) {
			uNum_fI[n]       = 0.5*(uIn_fI[n]+uOut_fIIn[n]);
			duNumduIn_fI[n]  = 0.5;
			duNumduOut_fI[n] = 0.5;
		}
		free(uOut_fIIn);

		// Include BC information in duNumduIn_fI if on a boundary
		if (Boundary) {
			duOutduIn = malloc(NfnI * sizeof *duOutduIn); // free

			jacobian_boundary_Poisson(NfnI,1,duOutduIn,BC % BC_STEP_SC);
			for (n = 0; n < NfnI; n++)
				duNumduIn_fI[n] += duNumduOut_fI[n]*duOutduIn[n];

			free(duOutduIn);
		}

		wnJ_fI = malloc(NfnI*d * sizeof *wnJ_fI); // free
		for (dim = 0; dim < d; dim++) {
		for (n = 0; n < NfnI; n++) {
			wnJ_fI[dim*NfnI+n] = w_fI[n]*detJF_fI[n]*n_fI[n*d+dim];
		}}

		// Compute partial FACET RHS and LHS terms

		// Interior VOLUME
		for (dim = 0; dim < d; dim++) {
			// RHSIn (partial)
			if (FACET->qhatIn[dim])
				free(FACET->qhatIn[dim]);
			FACET->qhatIn[dim] = malloc(NfnI * sizeof *(FACET->qhatIn[dim])); // keep
			for (n = 0; n < NfnI; n++) {
				FACET->qhatIn[dim][n] = wnJ_fI[dim*NfnI+n]*(uNum_fI[n]-uIn_fI[n]);
			}

			// LHSInIn (partial)
			if (FACET->qhat_uhatInIn[dim])
				free(FACET->qhat_uhatInIn[dim]);
			FACET->qhat_uhatInIn[dim] = calloc(NfnI*NvnSIn , sizeof *(FACET->qhat_uhatInIn[dim])); // keep

			mult_diag = malloc(NfnI * sizeof *mult_diag); // free
			for (n = 0; n < NfnI; n++)
				mult_diag[n] = wnJ_fI[dim*NfnI+n]*(duNumduIn_fI[n]-1.0);

			mm_diag_d(NfnI,NvnSIn,mult_diag,OPSIn->ChiS_fI[VfIn],FACET->qhat_uhatInIn[dim],1.0,'L','R');
			free(mult_diag);
		}
		free(uIn_fI);

		// Exterior VOLUME
		if (!Boundary) {
			ChiS_fIOutIn = malloc(NfnI*NvnSOut * sizeof *ChiS_fIOutIn); // free

			ChiS_fI = OPSOut->ChiS_fI[VfOut];
			for (i = 0; i < NfnI; i++) {
			for (j = 0; j < NvnSOut; j++) {
				ChiS_fIOutIn[i*NvnSOut+j] = ChiS_fI[nOrdOutIn[i]*NvnSOut+j];
			}}

			for (dim = 0; dim < d; dim++) {
				// LHSOutIn (partial)
				if (FACET->qhat_uhatOutIn[dim])
					free(FACET->qhat_uhatOutIn[dim]);
				FACET->qhat_uhatOutIn[dim] = calloc(NfnI*NvnSOut , sizeof *(FACET->qhat_uhatOutIn[dim])); // keep

				mult_diag = malloc(NfnI * sizeof *mult_diag); // free
				for (n = 0; n < NfnI; n++)
					mult_diag[n] = wnJ_fI[dim*NfnI+n]*(duNumduOut_fI[n]);

				mm_diag_d(NfnI,NvnSOut,mult_diag,ChiS_fIOutIn,FACET->qhat_uhatOutIn[dim],1.0,'L','R');
				free(mult_diag);
			}
			free(ChiS_fIOutIn);

			// Use "-ve" normal for opposite VOLUME
			for (n = 0; n < d*NfnI; n++) {
				wnJ_fI[n] *= -1.0;
			}

			array_rearrange_d(NfnI,d,nOrdInOut,'C',wnJ_fI);
			array_rearrange_d(NfnI,1,nOrdInOut,'R',uNum_fI);
			array_rearrange_d(NfnI,1,nOrdInOut,'R',duNumduIn_fI);
			array_rearrange_d(NfnI,1,nOrdInOut,'R',duNumduOut_fI);

			ChiS_fIInOut = malloc(NfnI*NvnSIn * sizeof *ChiS_fIInOut); // free

			ChiS_fI = OPSIn->ChiS_fI[VfIn];
			for (i = 0; i < NfnI; i++) {
			for (j = 0; j < NvnSIn; j++) {
				ChiS_fIInOut[i*NvnSIn+j] = ChiS_fI[nOrdInOut[i]*NvnSIn+j];
			}}

			for (dim = 0; dim < d; dim++) {
				if (FACET->qhatOut[dim])
					free(FACET->qhatOut[dim]);
				FACET->qhatOut[dim] = malloc(NfnI * sizeof *(FACET->qhatOut[dim])); // keep
				// RHSOut (partial)
				for (n = 0; n < NfnI; n++) {
					FACET->qhatOut[dim][n] = wnJ_fI[dim*NfnI+n]*(uNum_fI[n]-uOut_fI[n]);
				}

				// LHSInOut (partial)
				if (FACET->qhat_uhatInOut[dim])
					free(FACET->qhat_uhatInOut[dim]);
				FACET->qhat_uhatInOut[dim] = calloc(NfnI*NvnSIn , sizeof *(FACET->qhat_uhatInOut[dim])); // keep

				mult_diag = malloc(NfnI * sizeof *mult_diag); // free
				for (n = 0; n < NfnI; n++)
					mult_diag[n] = wnJ_fI[dim*NfnI+n]*(duNumduIn_fI[n]);

				mm_diag_d(NfnI,NvnSIn,mult_diag,ChiS_fIInOut,FACET->qhat_uhatInOut[dim],1.0,'L','R');
				free(mult_diag);

				// LHSOutOut (partial)
				if (FACET->qhat_uhatOutOut[dim])
					free(FACET->qhat_uhatOutOut[dim]);
				FACET->qhat_uhatOutOut[dim] = calloc(NfnI*NvnSOut , sizeof *(FACET->qhat_uhatOutOut[dim])); // keep

				mult_diag = malloc(NfnI * sizeof *mult_diag); // free
				for (n = 0; n < NfnI; n++)
					mult_diag[n] = wnJ_fI[dim*NfnI+n]*(duNumduOut_fI[n]-1.0);

				mm_diag_d(NfnI,NvnSOut,mult_diag,OPSOut->ChiS_fI[VfOut],FACET->qhat_uhatOutOut[dim],1.0,'L','R');
				free(mult_diag);
			}
			free(ChiS_fIInOut);

			free(uOut_fI);
		}
		free(uNum_fI);
		free(duNumduIn_fI);
		free(duNumduOut_fI);
		free(wnJ_fI);
	}
	free(OPSIn);
	free(OPSOut);
}

static void finalize_qhat(void)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int dim, iMax, NvnSIn, NvnSOut, NfnI, VfIn, VfOut;
	double       *MInvChiS_fI;
	double       *VqhatIn_ptr, *FqhatIn_ptr, *VqhatOut_ptr, *FqhatOut_ptr, *Fqhat;

	struct S_OPERATORS *OPSIn, *OPSOut;
	struct S_FACET     *FACET;
	struct S_VOLUME    *VIn, *VOut;

	OPSIn  = malloc(sizeof *OPSIn);  // free
	OPSOut = malloc(sizeof *OPSOut); // free

	for (FACET = DB.FACET; FACET; FACET = FACET->next) {
		VIn    = FACET->VIn;
		VfIn   = FACET->VfIn;
		NvnSIn = VIn->NvnS;

		init_opsF(OPSIn,VIn,FACET,get_IndFType(VIn->Eclass,VfIn/NFREFMAX));

		NfnI = OPSIn->NfnI;

		MInvChiS_fI = mm_Alloc_d(CBRM,CBNT,CBT,NvnSIn,NfnI,NvnSIn,1.0,VIn->MInv,OPSIn->ChiS_fI[VfIn]); // free

		Fqhat = malloc(NvnSIn * sizeof *Fqhat); // free
		for (dim = 0; dim < d; dim++) {
			mm_d(CBCM,CBT,CBNT,NvnSIn,1,NfnI,1.0,0.0,MInvChiS_fI,FACET->qhatIn[dim],Fqhat);

			VqhatIn_ptr  = VIn->qhat[dim];
			FqhatIn_ptr  = Fqhat;

			for (iMax = NvnSIn; iMax--; )
				*VqhatIn_ptr++ += *FqhatIn_ptr++;
		}
		free(MInvChiS_fI);
		free(Fqhat);

		if(!(FACET->Boundary)) {
			VOut    = FACET->VOut;
			VfOut   = FACET->VfOut;
			NvnSOut = VOut->NvnS;

			init_opsF(OPSOut,VOut,FACET,get_IndFType(VOut->Eclass,VfOut/NFREFMAX));

			MInvChiS_fI = mm_Alloc_d(CBRM,CBNT,CBT,NvnSOut,NfnI,NvnSOut,1.0,VOut->MInv,OPSOut->ChiS_fI[VfOut]); // free

			Fqhat = malloc(NvnSOut * sizeof *Fqhat); // free
			for (dim = 0; dim < d; dim++) {
				mm_d(CBCM,CBT,CBNT,NvnSOut,1,NfnI,1.0,0.0,MInvChiS_fI,FACET->qhatOut[dim],Fqhat);

				VqhatOut_ptr = VOut->qhat[dim];
				FqhatOut_ptr = Fqhat;

				for (iMax = NvnSOut; iMax--; )
					*VqhatOut_ptr++ += *FqhatOut_ptr++;
			}
			free(MInvChiS_fI);
			free(Fqhat);
		}

		for (dim = 0; dim < d; dim++) {
			free(FACET->qhatIn[dim]);  FACET->qhatIn[dim]  = NULL;
			free(FACET->qhatOut[dim]); FACET->qhatOut[dim] = NULL;
		}
	}
	free(OPSIn);
	free(OPSOut);
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
                        const unsigned int P, double *gradu_avg, double *u_jump, const unsigned int d,
                        const  unsigned int flux_type)
{
	/*
	 *	Comments:
	 *		Potentially change so that only the determination of tau depends on flux_type. (ToBeModified)
	 *
	 *	References:
	 *		Brdar(2012)-Compact_and_Stable_Discontinuous_Galerkin_Methods_for_Convection-Diffusion_Problems
	 *
	 *		tau:
	 *			Shahbazi(2005)-An_Explicit_Expression_for_the_Penalty_Parameter_of_the_Interior_Penalty_Method
	 *			Hesthave(2008)-Nodal_Discontinuous_Galerkin_Methods (section 7.2)
	 *			Brdar(2012): Theorem 2 (b)
	 */

	unsigned int dim, n;
	double       tau;

	if (Nel != 1)
		printf("Error: Unsupported Nel.\n"), EXIT_MSG;

	switch (flux_type) {
	case FLUX_IP:
		for (dim = 0; dim < d; dim++) {
		for (n = 0; n < Nn; n++) {
			tau = 1e2*(P+1)*(P+1)/h[n];

			gradu_avg[Nn*dim+n] = 0.5;
			u_jump[Nn*dim+n]    = -tau*nIn[n*d+dim];
		}}
		break;
	case FLUX_BR2:
		tau = 1.0*(DB.NfMax);
		for (dim = 0; dim < d; dim++) {
		for (n = 0; n < Nn; n++) {
			gradu_avg[Nn*dim+n] = 0.5;
			u_jump[Nn*dim+n]    = -0.5*tau*nIn[n*d+dim];
		}}
		break;
	case FLUX_CDG2:
		tau = 0.5*(DB.NfMax);
		for (dim = 0; dim < d; dim++) {
		for (n = 0; n < Nn; n++) {
			gradu_avg[Nn*dim+n] = 0.5;
			u_jump[Nn*dim+n]    = -0.5*tau*nIn[n*d+dim];
		}}
		break;
	default:
		printf("Error: Unsupported flux_type.\n"), EXIT_MSG;
		break;
	}
}

static void compute_uhat_FACET()
{
	// Initialize DB Parameters
	unsigned int d               = DB.d,
	             ViscousFluxType = DB.ViscousFluxType;

	// Standard datatypes
	unsigned int i, j, n, dim, dim1,
	             NvnSIn, NvnSOut, NfnI,
	             BC, Boundary, VfIn, VfOut, fIn, EclassIn, IndFType,
	             *nOrdOutIn, *nOrdInOut;
	double       **GradxyzIn, **GradxyzOut, *ChiS_fI, *ChiS_fI_std, *ChiSMInv_fIIn, *ChiSMInv_fIOut, *I_jump,
	             *LHSInIn, *LHSOutIn, *LHSInOut, *LHSOutOut,
	             *detJVIn_fI, *detJVOut_fI, *h, *n_fI, *detJF_fI, *w_fI, VolL, VolR,
	             *uIn_fI, *grad_uIn_fI, *uOut_fIIn, *grad_uOut_fIIn, *uOut_fI,
	             *nqNum_fI, *dnqNumduhatIn_fI, *dnqNumduhatOut_fI, *duOutduIn,
	             *gradu_avg, *u_jump, *u_jump_part,
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

		mm_CTN_d(NfnI,1,NvnSIn,OPSIn->ChiS_fI[VfIn],VIn->uhat,uIn_fI);

		for (dim = 0; dim < d; dim++)
			mm_CTN_d(NfnI,1,NvnSIn,GradxyzIn[dim],VIn->uhat,&grad_uIn_fI[NfnI*dim]);

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
			boundary_Poisson(NfnI,1,FACET->XYZ_fI,uIn_fI,uOut_fIIn,grad_uIn_fI,grad_uOut_fIIn,BC % BC_STEP_SC,BC / BC_STEP_SC);
		}

		// Compute numerical flux and its Jacobians
		nqNum_fI = calloc(NfnI   , sizeof *nqNum_fI); // free

		gradu_avg = malloc(NfnI*d * sizeof *gradu_avg); // free
		u_jump    = malloc(NfnI*d * sizeof *u_jump);    // free

		jacobian_flux_coef(NfnI,1,n_fI,h,FACET->P,gradu_avg,u_jump,d,ViscousFluxType);

		// Modify gradu_avg and u_jump scaling terms if required
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

		for (n = 0; n < NfnI; n++) {
		for (dim = 0; dim < d; dim++) {
			nqNum_fI[n] += n_fI[n*d+dim] *
				(  gradu_avg[NfnI*dim+n] * (grad_uIn_fI[NfnI*dim+n] + grad_uOut_fIIn[NfnI*dim+n])
				 + u_jump[NfnI*dim+n]    * (uIn_fI[n] - uOut_fIIn[n]));
		}}

		free(grad_uIn_fI);
		free(grad_uOut_fIIn);
		free(uIn_fI);
		free(uOut_fIIn);

		dnqNumduhatIn_fI  = calloc(NfnI*NvnSIn  , sizeof *dnqNumduhatIn_fI);  // free
		dnqNumduhatOut_fI = calloc(NfnI*NvnSOut , sizeof *dnqNumduhatOut_fI); // free

		// InIn contributions
		ChiS_fI = OPSIn->ChiS_fI[VfIn];

		for (dim = 0; dim < d; dim++) {
			for (n = 0; n < NfnI; n++) {
			for (j = 0; j < NvnSIn; j++) {
				// row-major
				dnqNumduhatIn_fI[n*NvnSIn+j] += n_fI[n*d+dim]*
					(gradu_avg[NfnI*dim+n]*GradxyzIn[dim][n*NvnSIn+j] + u_jump[NfnI*dim+n]*ChiS_fI[n*NvnSIn+j]);
			}}
		}

		for (n = 0; n < NfnI*d; n++)
			u_jump[n] *= -1.0;

		// Note: Both dnqNumduhatIn and dnqNumduhatOut are ordered corresponding to VIn.
		if (!Boundary) {
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
					dnqNumduhatOut_fI[n*NvnSOut+j] += n_fI[n*d+dim] *
						(gradu_avg[NfnI*dim+n]*GradxyzOut[dim][n*NvnSOut+j] + u_jump[NfnI*dim+n]*ChiS_fI[n*NvnSOut+j]);
				}}
			}
			free(ChiS_fI);
		} else {
			// Include BC information in dnqNumduhatIn_fI

			// Note: This approach is only possible because duOutduIn is a constant here. Otherwise, it would be
			//       better to linearize the flux wrt the solution itself (not the coefficients), add the duOutduIn
			//       term, and then add the final contribution for the linearization (du/duhat).
			duOutduIn = malloc(NfnI * sizeof *duOutduIn); // free

			jacobian_boundary_Poisson(NfnI,1,duOutduIn,BC % BC_STEP_SC);

			// OutOut contribution (u)
			ChiS_fI = OPSIn->ChiS_fI[VfIn];

			// Note: dgraduOutduIn == -duOutduIn.
			for (dim = 0; dim < d; dim++) {
				for (n = 0; n < NfnI; n++) {
				for (j = 0; j < NvnSIn; j++) {
					dnqNumduhatIn_fI[n*NvnSIn+j] += n_fI[n*d+dim] *
						(- duOutduIn[n]*gradu_avg[NfnI*dim+n]*GradxyzOut[dim][n*NvnSIn+j]
				    	 + duOutduIn[n]*u_jump[NfnI*dim+n]*ChiS_fI[n*NvnSIn+j]);
				}}
			}
			free(duOutduIn);
		}
		free(gradu_avg);
		free(u_jump);
		free(h);

		// Multiply by area element and cubature weights
		for (n = 0; n < NfnI; n++) {
			nqNum_fI[n] *= detJF_fI[n]*w_fI[n];
			for (j = 0; j < NvnSIn; j++)
				dnqNumduhatIn_fI[n*NvnSIn+j] *= detJF_fI[n]*w_fI[n];
			for (j = 0; j < NvnSOut; j++)
				dnqNumduhatOut_fI[n*NvnSOut+j] *= detJF_fI[n]*w_fI[n];
		}

		// Finalize FACET RHS and LHS terms

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
			mm_d(CBCM,CBNT,CBNT,NvnSIn, 1,NfnI,-1.0,1.0,GradxyzIn[dim], FACET->qhatIn[dim], RHSIn);
			mm_d(CBRM,CBT,CBNT,NvnSIn, NvnSIn, NfnI,-1.0,1.0,GradxyzIn[dim], FACET->qhat_uhatInIn[dim],  LHSInIn);
			if (!Boundary) {
				array_rearrange_d(NfnI,NvnSOut,nOrdInOut,'R',GradxyzOut[dim]);

				mm_d(CBCM,CBNT,CBNT,NvnSOut,1,NfnI,-1.0,1.0,GradxyzOut[dim],FACET->qhatOut[dim],RHSOut);
				mm_d(CBRM,CBT,CBNT,NvnSIn, NvnSOut,NfnI,-1.0,1.0,GradxyzIn[dim], FACET->qhat_uhatOutIn[dim], LHSOutIn);
				mm_d(CBRM,CBT,CBNT,NvnSOut,NvnSIn, NfnI,-1.0,1.0,GradxyzOut[dim],FACET->qhat_uhatInOut[dim], LHSInOut);
				mm_d(CBRM,CBT,CBNT,NvnSOut,NvnSOut,NfnI,-1.0,1.0,GradxyzOut[dim],FACET->qhat_uhatOutOut[dim],LHSOutOut);
			}
		}
		array_free2_d(d,GradxyzIn);
		array_free2_d(d,GradxyzOut);

		// Add FACET contributions to RHS and LHS

		// Interior FACET

		mm_d(CBCM,CBNT,CBNT,NvnSIn,1,NfnI,1.0,1.0,OPSIn->ChiS_fI[VfIn],nqNum_fI,RHSIn);
		mm_d(CBRM,CBT,CBNT,NvnSIn,NvnSIn,NfnI,1.0,1.0,OPSIn->ChiS_fI[VfIn],dnqNumduhatIn_fI,LHSInIn);

		// Exterior FACET
		if (!Boundary) {
			// RHS

			// Use -ve normal for opposite VOLUME
			for (n = 0; n < NfnI; n++)
				nqNum_fI[n] *= -1.0;

			// Rearrange nqNum to match node ordering from VOut
			array_rearrange_d(NfnI,1,nOrdInOut,'C',nqNum_fI);

			mm_d(CBCM,CBNT,CBNT,NvnSOut,1,NfnI,1.0,1.0,OPSOut->ChiS_fI[VfOut],nqNum_fI,RHSOut);

			// LHS

			// OutIn
			mm_d(CBRM,CBT,CBNT,NvnSIn,NvnSOut,NfnI,1.0,1.0,OPSIn->ChiS_fI[VfIn],dnqNumduhatOut_fI,LHSOutIn);

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
			mm_d(CBRM,CBT,CBNT,NvnSOut,NvnSIn,NfnI,1.0,1.0,OPSOut->ChiS_fI[VfOut],dnqNumduhatIn_fI,LHSInOut);

			// OutOut
			mm_d(CBRM,CBT,CBNT,NvnSOut,NvnSOut,NfnI,1.0,1.0,OPSOut->ChiS_fI[VfOut],dnqNumduhatOut_fI,LHSOutOut);
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

	implicit_info_Poisson();
	maxRHS = finalize_LHS(&A,&b,&x,0);

//	MatView(A,PETSC_VIEWER_STDOUT_SELF);
//	VecView(b,PETSC_VIEWER_STDOUT_SELF);
//	EXIT_MSG;

	// Solve linear system
	KSPCreate(MPI_COMM_WORLD,&ksp);
	setup_KSP(A,ksp);

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

//	if (!DB.Testing)
		printf("KSP iterations (cond, reason): %5d (% .3e, %d)\n",iteration_ksp,emax/emin,reason);

	// Note: maxRHS is meaningless as it is based on the initial (zero) solution.
	if (0)
		printf("% .3e\n",maxRHS);

	free(fNameOut);
	free(string);

	// Potentially adaptation option (ToBeDeleted)
}
