// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "solver_Poisson.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>

#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "S_FACE.h"
#include "Test.h"

#include "solver_functions.h"
#include "update_VOLUMEs.h"
#include "implicit_GradW.h"

#include "element_functions.h"
#include "matrix_functions.h"
#include "exact_solutions.h"
#include "setup_geom_factors.h"
#include "array_swap.h"
#include "array_free.h"
#include "finalize_LHS.h"
#include "solver_implicit.h"
#include "output_to_paraview.h"
#include "setup_Curved.h"

#include "array_print.h"
#include "array_norm.h" // ToBeDeleted

/*
 *	Purpose:
 *		Perform the implicit solve for the Poisson equation.
 *
 *	Comments:
 *		Many of the RHS terms computed are 0; they are included as they are used to check the linearization. Further,
 *		the computational cost is dominated by the global system solve making this additional cost negligible.
 *
 *	Notation:
 *
 *	References:
 */

struct S_OPERATORS {
	unsigned int NvnS, NfnI, NvnC,
	             *nOrdOutIn, *nOrdInOut;
	double       *w_fI, *ChiS_vI, **ChiS_fI, ***GradChiS_fI, **I_vC_fI, **D_Weak;
};

static void init_opsF(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const struct S_FACE *FACE,
                      const unsigned int IndFType)
{
	// Standard datatypes
	unsigned int PV, PF, Vtype, Vcurved, FtypeInt, IndOrdOutIn, IndOrdInOut;

	struct S_ELEMENT *ELEMENT, *ELEMENT_FACE;

	PV      = VOLUME->P;
	PF      = FACE->P;
	Vtype   = VOLUME->type;
	Vcurved = VOLUME->curved;

	FtypeInt = FACE->typeInt;
	IndOrdOutIn = FACE->IndOrdOutIn;
	IndOrdInOut = FACE->IndOrdInOut;

	ELEMENT       = get_ELEMENT_type(Vtype);
	ELEMENT_FACE = get_ELEMENT_FACE(Vtype,IndFType);

	OPS->NvnS = ELEMENT->NvnS[PV];
	if (FtypeInt == 's') {
		// Straight FACE Integration
		OPS->NfnI = ELEMENT->NfnIs[PF][IndFType];

		OPS->w_fI = ELEMENT->w_fIs[PF][IndFType];

		OPS->ChiS_fI     = ELEMENT->ChiS_fIs[PV][PF];
		OPS->GradChiS_fI = ELEMENT->GradChiS_fIs[PV][PF];

		OPS->nOrdInOut = ELEMENT_FACE->nOrd_fIs[PF][IndOrdInOut];
		OPS->nOrdOutIn = ELEMENT_FACE->nOrd_fIs[PF][IndOrdOutIn];

		if (!Vcurved) {
			OPS->NvnC = ELEMENT->NvnCs[PV];
			OPS->I_vC_fI = ELEMENT->I_vCs_fIs[PV][PF];
		} else {
			OPS->NvnC = ELEMENT->NvnCc[PV];
			OPS->I_vC_fI = ELEMENT->I_vCc_fIs[PV][PF];
		}
	} else {
		// Curved FACE Integration
		OPS->NfnI = ELEMENT->NfnIc[PF][IndFType];

		OPS->w_fI = ELEMENT->w_fIc[PF][IndFType];

		OPS->ChiS_fI     = ELEMENT->ChiS_fIc[PV][PF];
		OPS->GradChiS_fI = ELEMENT->GradChiS_fIc[PV][PF];

		OPS->nOrdInOut = ELEMENT_FACE->nOrd_fIc[PF][IndOrdInOut];
		OPS->nOrdOutIn = ELEMENT_FACE->nOrd_fIc[PF][IndOrdOutIn];

		if (!Vcurved) {
			OPS->NvnC = ELEMENT->NvnCs[PV];
			OPS->I_vC_fI = ELEMENT->I_vCs_fIc[PV][PF];
		} else {
			OPS->NvnC = ELEMENT->NvnCc[PV];
			OPS->I_vC_fI = ELEMENT->I_vCc_fIc[PV][PF];
		}
	}
}

static void compute_Qhat(void)
{
	/*
	 *	Purpose:
	 *		Compute the weak gradients and store LHSQ for use below.
	 */

	implicit_GradW_VOLUME();

	unsigned int d = DB.d;

	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		unsigned int const NvnS = VOLUME->NvnS;

		// Store contribution before multiplication by MInv
		for (size_t dim = 0; dim < d; dim++) {
			for (size_t i = 0; i < NvnS*NvnS; i++)
				VOLUME->LHSQ[dim][i] = VOLUME->QhatV_What[dim][i];
		}
	}

	implicit_GradW_FACE();
	implicit_GradW_finalize();
}

void project_to_sphere(const unsigned int Nn, double *XYZIn, double *XYZOut, const unsigned int BCcurved)
{
	printf("Error: Use compute_normal_displacement instead.\n"), EXIT_MSG;
	/*
	 *	Purpose:
	 *		Project coordinates to the sphere surface.
	 */

	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int n;
	double       r, r2, t, p, rIn, rOut, norm_rIn, norm_rOut, *XIn, *YIn, *ZIn, *XOut, *YOut, *ZOut;

	rIn  = DB.rIn;
	rOut = DB.rOut;

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

	if (BCcurved == 1) {
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
			array_print_d(Nn,d,XYZIn,'C');
			printf("Error: Found unsupported curved BC.\n"), EXIT_MSG;
		}

		for (n = 0; n < Nn; n++) {
			t = atan2(YIn[n],XIn[n]);
			XOut[n] = r*cos(t);
			YOut[n] = r*sin(t);

			if (d == 3) {
				// Note: Interpolation to FACE cubature nodes sometimes resulted in Z > r resulting in the computed p
				//       being imaginary.
				if (ZIn[n]-r > 0.0)
					p = acos(1.0);
				else
					p = acos(ZIn[n]/r);

				XOut[n] *= sin(p);
				YOut[n] *= sin(p);
				ZOut[n]  = r*cos(p);
			}
		}
	}
}

void boundary_Poisson(const unsigned int Nn, const unsigned int Nel, double *XYZ, double *normals, double *uL,
                      double *uR, double *graduL, double *graduR, const unsigned int BC, const unsigned int BCcurved)
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
	XYZproj = calloc(Nn*d , sizeof *XYZproj); // free
	if (d > 1) {
		if (d == 2) {
			if (BCcurved == 2)
				compute_normal_displacement(Nn,1,XYZ,normals,XYZproj,BC);
			for (n = 0; n < d*Nn; n++)
				XYZproj[n] += XYZ[n];
		} else if (d == 3) {
			// Using normal projection in 3D will result in gaps in the projected mesh. Problem?
			// Imagine projection of EDGE node of two adjacent TETs to the spherical boundary.
			printf("Try with normal projection.\n"), EXIT_MSG;
			project_to_sphere(Nn,XYZ,XYZproj,BCcurved);
		}
	} else {
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}

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

static void compute_What_VOLUME(void)
{
	unsigned int const d    = DB.d,
	                   Neq  = 1,
	                   Nvar = 1;

	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		unsigned int const NvnS = VOLUME->NvnS;

		// Compute RHS and LHS terms
		double **LHSQ = VOLUME->LHSQ;

		// RHS
		// Note: Using CBCM and CBNT for LHSQ results in LHSQ'
		memset(VOLUME->RHS,0.0,NvnS*Neq * sizeof *(VOLUME->RHS));
		for (size_t dim = 0; dim < d; dim++)
			mm_d(CBCM,CBNT,CBNT,NvnS,1,NvnS,-1.0,1.0,LHSQ[dim],VOLUME->QhatV[dim],VOLUME->RHS);
// Potentially modify this to be Qhat instead of QhatV. (ToBeDeleted)

		// LHS
		memset(VOLUME->LHS,0.0,NvnS*NvnS*Neq*Nvar * sizeof *(VOLUME->LHS));
		for (size_t dim = 0; dim < d; dim++)
			mm_d(CBRM,CBT,CBNT,NvnS,NvnS,NvnS,-1.0,1.0,LHSQ[dim],VOLUME->QhatV_What[dim],VOLUME->LHS);
	}
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
	 *			Hesthaven(2008)-Nodal_Discontinuous_Galerkin_Methods (section 7.2)
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
			tau = 1e1*(P+1)*(P+1)/h[n];

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

static void compute_What_FACE()
{
	// Initialize DB Parameters
	unsigned int d               = DB.d,
	             ViscousFluxType = DB.ViscousFluxType;

	// Standard datatypes
	unsigned int i, j, n, dim, dim1,
	             NvnSIn, NvnSOut,
	             BC, Boundary, VfIn, VfOut, fIn, EclassIn,
	             *nOrdOutIn, *nOrdInOut;
	double       **GradxyzIn, **GradxyzOut, *ChiS_fI, *ChiS_fI_std, *ChiSMInv_fIIn, *ChiSMInv_fIOut, *I_jump,
	             *LHSInIn, *LHSOutIn, *LHSInOut, *LHSOutOut,
	             *detJVIn_fI, *detJVOut_fI, *h, *n_fI, *detJF_fI, *w_fI, VolL, VolR,
	             *uIn_fI, *grad_uIn_fI, *uOut_fIIn, *grad_uOut_fIIn, *uOut_fI,
	             *nqNum_fI, *dnqNumduhatIn_fI, *dnqNumduhatOut_fI, *duOutduIn,
	             *gradu_avg, *u_jump, *u_jump_part,
	             *RHSIn, *RHSOut;

	unsigned int const Nvar = DB.Nvar,
	                   Neq  = DB.Neq;

	struct S_OPERATORS *OPSIn, *OPSOut;
	struct S_VOLUME    *VIn, *VOut;

	OPSIn  = malloc(sizeof *OPSIn);  // free
	OPSOut = malloc(sizeof *OPSOut); // free

	struct S_OPERATORS_F *OPSL[2], *OPSR[2];

	struct S_FDATA *FDATAL = malloc(sizeof *FDATAL), // free
	               *FDATAR = malloc(sizeof *FDATAR); // free
	FDATAL->OPS = (struct S_OPERATORS_F const *const *) OPSL;
	FDATAR->OPS = (struct S_OPERATORS_F const *const *) OPSR;

	struct S_NumericalFlux *NFluxData = malloc(sizeof *NFluxData); // free
	FDATAL->NFluxData = NFluxData;
	FDATAR->NFluxData = NFluxData;

	for (size_t i = 0; i < 2; i++) {
		OPSL[i] = malloc(sizeof *OPSL[i]); // free
		OPSR[i] = malloc(sizeof *OPSR[i]); // free
	}

	for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
		init_FDATA(FDATAL,FACE,'L');
		init_FDATA(FDATAR,FACE,'R');

		// Compute WL_fIL, WR_fIL, QpL_fIL, and QpR_fIL (i.e. as seen from the (L)eft VOLUME)
//		unsigned int const IndFType = FDATAL->IndFType,
		unsigned int       IndFType = FDATAL->IndFType,
		                   NfnI     = OPSL[IndFType]->NfnI;

		FDATAL->W_fIL = malloc(NfnI*Nvar * sizeof *(FDATAL->W_fIL)), // free
		FDATAR->W_fIL = malloc(NfnI*Nvar * sizeof *(FDATAR->W_fIL)); // free

		double **const QpL_fIL = malloc(d * sizeof *QpL_fIL), // free
		       **const QpR_fIL = malloc(d * sizeof *QpR_fIL); // free

		for (size_t dim = 0; dim < d; dim++) {
			QpL_fIL[dim] = malloc(NfnI*Nvar * sizeof *QpL_fIL[dim]); // free
			QpR_fIL[dim] = malloc(NfnI*Nvar * sizeof *QpR_fIL[dim]); // free
		}

		FDATAL->Qp_fIL = QpL_fIL;
		FDATAR->Qp_fIL = QpR_fIL;

		coef_to_values_fI(FDATAL,'W','I');
		coef_to_values_fI(FDATAL,'Q','I');
		compute_WR_QpR_fIL(FDATAR,FDATAL->W_fIL,FDATAR->W_fIL,(double const *const *const) QpL_fIL,QpR_fIL,'I');

		// Compute numerical flux and its Jacobian as seen from the left VOLUME
		NFluxData->WL_fIL              = FDATAL->W_fIL;
		NFluxData->WR_fIL              = FDATAR->W_fIL;
		NFluxData->nFluxViscNum_fI     = malloc(NfnI*Neq      * sizeof *(NFluxData->nFluxViscNum_fI));     // free
		NFluxData->dnFluxViscNumdWL_fI = malloc(NfnI*Neq*Nvar * sizeof *(NFluxData->dnFluxViscNumdWL_fI)); // free
		NFluxData->dnFluxViscNumdWR_fI = malloc(NfnI*Neq*Nvar * sizeof *(NFluxData->dnFluxViscNumdWR_fI)); // free
		NFluxData->dnFluxViscNumdQL_fI = malloc(d             * sizeof *(NFluxData->dnFluxViscNumdQL_fI)); // free
		NFluxData->dnFluxViscNumdQR_fI = malloc(d             * sizeof *(NFluxData->dnFluxViscNumdQR_fI)); // free

		for (size_t dim = 0; dim < d; dim++) {
			NFluxData->dnFluxViscNumdQL_fI[dim] = malloc(NfnI*Neq*Nvar * sizeof *(NFluxData->dnFluxViscNumdQL_fI[dim])); // free
			NFluxData->dnFluxViscNumdQR_fI[dim] = malloc(NfnI*Neq*Nvar * sizeof *(NFluxData->dnFluxViscNumdQR_fI[dim])); // free
		}

		compute_numerical_flux_viscous(FDATAL,FDATAR,'I');




		setup_geom_factors_highorder(FACE);

		VIn  = FACE->VIn;
		VfIn = FACE->VfIn;
		fIn  = VfIn/NFREFMAX;
if (0) {
		EclassIn = VIn->Eclass;
		IndFType = get_IndFType(EclassIn,fIn);
}
		init_opsF(OPSIn,VIn,FACE,IndFType);

		VOut  = FACE->VOut;
		VfOut = FACE->VfOut;

		init_opsF(OPSOut,VOut,FACE,IndFType);

		BC       = FACE->BC;
		Boundary = FACE->Boundary;

//		NfnI    = OPSIn->NfnI;
		NvnSIn  = OPSIn->NvnS;
		NvnSOut = OPSOut->NvnS;

		w_fI = OPSIn->w_fI;

		nOrdOutIn = OPSIn->nOrdOutIn;
		nOrdInOut = OPSIn->nOrdInOut;

		n_fI = FACE->n_fI;

		detJVIn_fI = FACE->detJVIn_fI;
		if (!Boundary) {
			detJVOut_fI = malloc(NfnI * sizeof *detJVOut_fI); // free
			for (n = 0; n < NfnI; n++)
				detJVOut_fI[n] = FACE->detJVOut_fI[nOrdOutIn[n]];
		} else {
			detJVOut_fI = detJVIn_fI;
		}

		detJF_fI = FACE->detJF_fI;
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
		// Alternative computation of GradxyzIn/Out commented below (No difference for straight VOLUMEs)
		mm_d(CBRM,CBNT,CBNT,NfnI,NvnSIn, NvnSIn, 1.0,0.0,OPSIn->ChiS_fI[VfIn],  VIn->MInv, ChiSMInv_fIIn);
		mm_d(CBRM,CBNT,CBNT,NfnI,NvnSOut,NvnSOut,1.0,0.0,OPSOut->ChiS_fI[VfOut],VOut->MInv,ChiSMInv_fIOut);
		array_rearrange_d(NfnI,NvnSOut,nOrdOutIn,'R',ChiSMInv_fIOut);

		// Note: GradxyzOut is arranged according to the FACE node ordering of VIn.
		for (dim1 = 0; dim1 < d; dim1++) {
			GradxyzIn[dim1]  = malloc(NfnI*NvnSIn  * sizeof **GradxyzIn);  // free
			GradxyzOut[dim1] = malloc(NfnI*NvnSOut * sizeof **GradxyzOut); // free

			mm_d(CBRM,CBNT,CBNT,NfnI,NvnSIn, NvnSIn, 1.0,0.0,ChiSMInv_fIIn, VIn->LHSQ[dim1], GradxyzIn[dim1]);
			mm_d(CBRM,CBNT,CBNT,NfnI,NvnSOut,NvnSOut,1.0,0.0,ChiSMInv_fIOut,VOut->LHSQ[dim1],GradxyzOut[dim1]);
		}

		mm_CTN_d(NfnI,1,NvnSIn,OPSIn->ChiS_fI[VfIn],VIn->What,uIn_fI);

		for (dim = 0; dim < d; dim++)
			mm_CTN_d(NfnI,1,NvnSIn,GradxyzIn[dim],VIn->What,&grad_uIn_fI[NfnI*dim]);
// matches ChiS_fI*QhatV (OK)

		// Compute_uOut_fI (Taking BCs into account if applicable)
		uOut_fIIn      = malloc(NfnI   * sizeof *uOut_fIIn);      // free
		grad_uOut_fIIn = calloc(NfnI*d , sizeof *grad_uOut_fIIn); // free

		if (!Boundary) {
			uOut_fI = malloc(NfnI * sizeof *uOut_fI); // free

			mm_CTN_d(NfnI,1,NvnSOut,OPSOut->ChiS_fI[VfOut],VOut->What,uOut_fI);

			// Reorder uOut_fI to correspond to inner VOLUME ordering
			for (n = 0; n < NfnI; n++)
				uOut_fIIn[n] = uOut_fI[nOrdOutIn[n]];
			free(uOut_fI);

			for (dim = 0; dim < d; dim++)
				mm_CTN_d(NfnI,1,NvnSOut,GradxyzOut[dim],VOut->What,&grad_uOut_fIIn[NfnI*dim]);
		} else {
// In solver_functions, grad_uOut_fIIn is computed based on the corrected grad_uIn_fI
			boundary_Poisson(NfnI,1,FACE->XYZ_fI,n_fI,uIn_fI,uOut_fIIn,grad_uIn_fI,grad_uOut_fIIn,
			                 BC % BC_STEP_SC,BC / BC_STEP_SC);
		}

		// Compute numerical flux and its Jacobians
		nqNum_fI = calloc(NfnI   , sizeof *nqNum_fI); // free

		gradu_avg = malloc(NfnI*d * sizeof *gradu_avg); // free
		u_jump    = malloc(NfnI*d * sizeof *u_jump);    // free

		jacobian_flux_coef(NfnI,1,n_fI,h,FACE->P,gradu_avg,u_jump,d,ViscousFluxType);

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
/*
double *Wdiff = malloc(NfnI * sizeof *Wdiff); // ToBeDeleted
for (size_t n = 0; n < NfnI; n++)
	Wdiff[n] = uIn_fI[n]-uOut_fIIn[n];

double **qF = malloc(d * sizeof *qF); // Delete
double *qpF = malloc(NfnI * sizeof *qpF); // Delete
double *QpF = malloc(NfnI*d * sizeof *QpF);
for (size_t dim = 0; dim < d; dim++) {
	for (size_t n = 0; n < NfnI; n++)
		qpF[n] = w_fI[n]*detJF_fI[n]*n_fI[n*d+dim]*-0.5*Wdiff[n];

	qF[dim] = malloc(NvnSIn * sizeof *qF[dim]); // Delete
	mm_d(CBCM,CBNT,CBNT,NvnSIn,1,NfnI,1.0,0.0,ChiS_fI,qpF,qF[dim]);

	double *tmp_d = mm_Alloc_d(CBCM,CBT,CBNT,NvnSIn,1,NvnSIn,1.0,VIn->MInv,qF[dim]);

	printf("619: %zu % .3e \n",dim,array_norm_diff_d(NvnSIn,tmp_d,FACE->QhatL[dim],"Inf"));
//	array_print_d(NvnSIn,1,tmp_d,'C');
//	array_print_d(NvnSIn,1,FACE->QhatL[dim],'C');

	mm_d(CBCM,CBT,CBNT,NfnI,1,NfnI,1.0*DB.NfMax,0.0,I_jump,qpF,&QpF[dim*NfnI]);
}

array_print_d(NfnI,d,QpF,'C');
for (size_t dim = 0; dim < d; dim++) {
	mm_d(CBCM,CBT,CBNT,NfnI,1,NvnSIn,1.0*DB.NfMax,0.0,ChiS_fI,FACE->QhatL[dim],&QpF[dim*NfnI]);
}
array_print_d(NfnI,d,QpF,'C');
for (size_t dim = 0; dim < d; dim++)
	array_print_d(NfnI,1,FDATAL->Qp_fIL[dim],'C');


//EXIT_BASIC;
*/
				ChiS_fI_std = OPSOut->ChiS_fI[VfOut];
				ChiS_fI = malloc(NfnI*NvnSOut * sizeof *ChiS_fI); // free

				// Reordering
				for (i = 0; i < NfnI; i++) {
				for (j = 0; j < NvnSOut; j++) {
					ChiS_fI[i*NvnSOut+j] = ChiS_fI_std[nOrdInOut[i]*NvnSOut+j];
				}}
				mm_d(CBRM,CBNT,CBT,NfnI,NfnI,NvnSOut,1.0,1.0,ChiSMInv_fIOut,ChiS_fI,I_jump);
/*
for (size_t dim = 0; dim < d; dim++) {
	for (size_t n = 0; n < NfnI; n++)
//		qpF[n] = w_fI[n]*detJF_fI[n]*n_fI[n*d+dim]*-0.5*Wdiff[n];
		qpF[n] = w_fI[n]*detJF_fI[n]*n_fI[n*d+dim]*-0.5;

	qF[dim] = malloc(NvnSOut * sizeof *qF[dim]); // Delete
	mm_d(CBCM,CBNT,CBNT,NvnSOut,1,NfnI,1.0,0.0,ChiS_fI,qpF,qF[dim]);

	double *tmp_d = mm_Alloc_d(CBCM,CBT,CBNT,NvnSIn,1,NvnSIn,1.0,VOut->MInv,qF[dim]);

	printf("640: %zu % .3e\n",dim,array_norm_diff_d(NvnSOut,tmp_d,FACE->QhatR[dim],"Inf"));
//	array_print_d(NvnSOut,1,tmp_d,'C');
//	array_print_d(NvnSOut,1,FACE->QhatR[dim],'C');

	mm_d(CBCM,CBT,CBNT,NfnI,1,NfnI,1.0*DB.NfMax,0.0,I_jump,qpF,&QpF[dim*NfnI]);
}

//EXIT_BASIC;

printf("Used\n");
array_print_d(NfnI,d,QpF,'C');
for (size_t dim = 0; dim < d; dim++) {
	mm_d(CBCM,CBT,CBNT,NfnI,1,NvnSOut,1.0*DB.NfMax,1.0,ChiS_fI,FACE->QhatR[dim],&QpF[dim*NfnI]);
//	mm_d(CBCM,CBT,CBNT,NfnI,1,NvnSOut,1.0*DB.NfMax,0.0,ChiS_fI,FACE->QhatR[dim],&QpF[dim*NfnI]);
}
array_print_d(NfnI,d,QpF,'C');
for (size_t dim = 0; dim < d; dim++)
	array_print_d(NfnI,1,FDATAR->Qp_fIL[dim],'C');

for (size_t dim = 0; dim < d; dim++) {
for (size_t n = 0; n < NfnI; n++)
//	QpF[dim*NfnI+n] = 0.5*(FDATAL->Qp_fIL[dim][n]+FDATAR->Qp_fIL[dim][n]);
	QpF[dim*NfnI+n] = 1.0*(FDATAL->Qp_fIL[dim][n]+FDATAR->Qp_fIL[dim][n]);
}
array_print_d(NfnI,d,QpF,'C');
*/

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
//printf("u_jump\n");
//array_print_d(NfnI,d,u_jump,'C');

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
/*
double *QpF2 = malloc(NfnI*d * sizeof *QpF2);
for (size_t dim = 0; dim < d; dim++) {
	for (size_t n = 0; n < NfnI; n++)
//		QpF2[dim*NfnI+n] = u_jump[NfnI*dim+n]*(uIn_fI[n]-uOut_fIIn[n]);
		QpF2[dim*NfnI+n] = u_jump[NfnI*dim+n];
}
array_print_d(NfnI,d,QpF2,'C');
*/
		for (n = 0; n < NfnI; n++) {
		for (dim = 0; dim < d; dim++) {
// this should be what is returned by compute_numerical_flux_viscous
			nqNum_fI[n] += n_fI[n*d+dim] *
				(  gradu_avg[NfnI*dim+n] * (grad_uIn_fI[NfnI*dim+n] + grad_uOut_fIIn[NfnI*dim+n])
				 + u_jump[NfnI*dim+n]    * (uIn_fI[n] - uOut_fIIn[n]));
		}}
/*
printf("725: %d %d %d\n",FACE->indexg,Boundary,BC);
array_print_d(NfnI,1,nqNum_fI,'C');
array_print_d(NfnI,1,NFluxData->nFluxViscNum_fI,'C');

if (!Boundary)
EXIT_BASIC;
*/

bool useAlternate = 0;

if (useAlternate) {
for (size_t n = 0; n < NfnI; n++)
	nqNum_fI[n] = NFluxData->nFluxViscNum_fI[n];
}

		free(grad_uIn_fI);
		free(grad_uOut_fIIn);
		free(uIn_fI);
		free(uOut_fIIn);

		dnqNumduhatIn_fI  = calloc(NfnI*NvnSIn  , sizeof *dnqNumduhatIn_fI);  // free
		dnqNumduhatOut_fI = calloc(NfnI*NvnSOut , sizeof *dnqNumduhatOut_fI); // free

		// InIn contributions
		ChiS_fI = OPSIn->ChiS_fI[VfIn];

		for (dim = 0; dim < d; dim++) {
if (!useAlternate) {
			for (n = 0; n < NfnI; n++) {
			for (j = 0; j < NvnSIn; j++) {
				// row-major
				dnqNumduhatIn_fI[n*NvnSIn+j] += n_fI[n*d+dim]*
					(gradu_avg[NfnI*dim+n]*GradxyzIn[dim][n*NvnSIn+j] + u_jump[NfnI*dim+n]*ChiS_fI[n*NvnSIn+j]);
			}}
} else {
// Can likely use compute_LHS_FACE_Q_Weak instead (ToModified)
mm_diag_d(NfnI,NvnSIn,NFluxData->dnFluxViscNumdQL_fI[dim],FDATAL->Qp_WhatL[dim],dnqNumduhatIn_fI,1.0,1.0,'L','R');
}
		}
		array_free2_d(d,GradxyzIn);

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
if (!useAlternate) {
				for (n = 0; n < NfnI; n++) {
				for (j = 0; j < NvnSOut; j++) {
					dnqNumduhatOut_fI[n*NvnSOut+j] += n_fI[n*d+dim] *
						(gradu_avg[NfnI*dim+n]*GradxyzOut[dim][n*NvnSOut+j] + u_jump[NfnI*dim+n]*ChiS_fI[n*NvnSOut+j]);
				}}
} else {
// Can likely use compute_LHS_FACE_Q_Weak instead (ToModified)
//mm_diag_d(NfnI,NvnSOut,NFluxData->dnFluxViscNumdQR_fI[dim],FDATAR->Qp_WhatL[dim],dnqNumduhatOut_fI,1.0,1.0,'L','R');
mm_diag_d(NfnI,NvnSOut,NFluxData->dnFluxViscNumdQR_fI[dim],FDATAL->Qp_WhatR[dim],dnqNumduhatOut_fI,1.0,1.0,'L','R');
}
			}
			free(ChiS_fI);
} else if (!useAlternate) {
//		} else {
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
		array_free2_d(d,GradxyzOut);

		// Multiply by area element and cubature weights
		for (n = 0; n < NfnI; n++) {
			nqNum_fI[n] *= detJF_fI[n]*w_fI[n];
			for (j = 0; j < NvnSIn; j++)
				dnqNumduhatIn_fI[n*NvnSIn+j] *= detJF_fI[n]*w_fI[n];
			for (j = 0; j < NvnSOut; j++)
				dnqNumduhatOut_fI[n*NvnSOut+j] *= detJF_fI[n]*w_fI[n];
		}

		// Finalize FACE RHS and LHS terms

		// Add VOLUME contributions to RHS and LHS
		RHSIn   = FACE->RHSIn;
		LHSInIn = FACE->LHSInIn;

		memset(FACE->RHSIn,0.0,NvnSIn * sizeof *FACE->RHSIn);
		memset(FACE->LHSInIn,0.0,NvnSIn*NvnSIn * sizeof *FACE->LHSInIn);

		RHSOut = LHSOutIn = LHSInOut = LHSOutOut = NULL;
		if (!Boundary) {
			RHSOut    = FACE->RHSOut;
			LHSOutIn  = FACE->LHSOutIn;
			LHSInOut  = FACE->LHSInOut;
			LHSOutOut = FACE->LHSOutOut;

			memset(FACE->RHSOut,0.0,NvnSOut * sizeof *FACE->RHSOut);
			memset(FACE->LHSOutIn,0.0,NvnSIn*NvnSOut * sizeof *FACE->LHSOutIn);
			memset(FACE->LHSInOut,0.0,NvnSOut*NvnSIn * sizeof *FACE->LHSInOut);
			memset(FACE->LHSOutOut,0.0,NvnSOut*NvnSOut * sizeof *FACE->LHSOutOut);
		}

		for (dim = 0; dim < d; dim++) {
			mm_d(CBCM,CBNT,CBNT,NvnSIn,1,NvnSIn,-1.0,1.0,VIn->LHSQ[dim], FACE->QhatL[dim], RHSIn);
			mm_d(CBRM,CBT,CBNT,NvnSIn,NvnSIn,NvnSIn,-1.0,1.0,VIn->LHSQ[dim],FACE->Qhat_WhatLL[dim],LHSInIn);
			if (!Boundary) {
				mm_d(CBCM,CBNT,CBNT,NvnSOut,1,NvnSOut,-1.0,1.0,VOut->LHSQ[dim],FACE->QhatR[dim],RHSOut);
				mm_d(CBRM,CBT,CBNT,NvnSIn, NvnSOut,NvnSIn,-1.0,1.0,VIn->LHSQ[dim],FACE->Qhat_WhatRL[dim],LHSOutIn);
				mm_d(CBRM,CBT,CBNT,NvnSOut,NvnSIn,NvnSOut,-1.0,1.0,VOut->LHSQ[dim],FACE->Qhat_WhatLR[dim],LHSInOut);
				mm_d(CBRM,CBT,CBNT,NvnSOut,NvnSOut,NvnSOut,-1.0,1.0,VOut->LHSQ[dim],FACE->Qhat_WhatRR[dim],LHSOutOut);
			}
		}

		// Add FACE contributions to RHS and LHS

		// Interior FACE

		mm_d(CBCM,CBNT,CBNT,NvnSIn,1,NfnI,1.0,1.0,OPSIn->ChiS_fI[VfIn],nqNum_fI,RHSIn);
		mm_d(CBRM,CBT,CBNT,NvnSIn,NvnSIn,NfnI,1.0,1.0,OPSIn->ChiS_fI[VfIn],dnqNumduhatIn_fI,LHSInIn);

		// Exterior FACE
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

	for (size_t i = 0; i < 2; i++) {
		free(OPSL[i]);
		free(OPSR[i]);
	}

	free(FDATAL);
	free(FDATAR);
	free(NFluxData);
}

void implicit_info_Poisson(void)
{
	update_VOLUME_Ops();
	update_memory_VOLUMEs();

	compute_Qhat();

	compute_What_VOLUME();
	compute_What_FACE();
}

void solver_Poisson(bool PrintEnabled)
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
		uhat = VOLUME->What;

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

	// Update Qhat based on computed solution
	implicit_GradW();

	// Output to paraview
	if (TestDB.ML <= 1 || (TestDB.PGlobal == 1) || (TestDB.PGlobal == 5 && TestDB.ML <= 4)) {
		strcpy(fNameOut,"SolFinal_");
		sprintf(string,"%dD_",DB.d);   strcat(fNameOut,string);
		                               strcat(fNameOut,DB.MeshType);
		if (DB.Adapt == ADAPT_0) {
			sprintf(string,"_ML%d",DB.ML); strcat(fNameOut,string);
			sprintf(string,"P%d_",DB.PGlobal); strcat(fNameOut,string);
		} else {
			sprintf(string,"_ML%d",TestDB.ML); strcat(fNameOut,string);
			sprintf(string,"P%d_",TestDB.PGlobal); strcat(fNameOut,string);
		}
		output_to_paraview(fNameOut);
	}

//	if (!DB.Testing)
	if (1||PrintEnabled)
		printf("KSP iterations (cond, reason): %5d (% .3e, %d)\n",iteration_ksp,emax/emin,reason);

	// Note: maxRHS is meaningless as it is based on the initial (zero) solution.
	if (0)
		printf("% .3e\n",maxRHS);

	free(fNameOut);
	free(string);

	// Potentially adaptation option (ToBeDeleted)
}
