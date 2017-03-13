// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "explicit_FACE_info.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "S_FACE.h"
#include "S_OpCSR.h"

#include "element_functions.h"
#include "sum_factorization.h"
#include "matrix_functions.h"
#include "boundary_conditions.h"
#include "fluxes_inviscid.h"
#include "array_swap.h"

/*
 *	Purpose:
 *		Evaluate the FACE contributions to the RHS term.
 *
 *	Comments:
 *		When adaptation is enabled, the method used here to lift the FACE information to the VOLUME remains quite
 *		similar to that of the conforming case. The normal numerical flux is computed on each of the FACEs of the mesh
 *		and is used directly to form RHSIn/Out terms to be added to the appropriate VOLUME. This is in contrast to the
 *		traditional mortar element method (based on my current understanding), which first uses an L2 projection of the
 *		normal numerical flux of all non-conforming FACEs to standard VOLUME FACEs and then computes RHSIn/Out exactly
 *		as if the mesh were conforming. As no L2 projection is required for the current approach, compute_FACE_RHS_EFE
 *		can in fact be called even for the case of non-conforming discretizations allowing for:
 *			1) Reduced cost because the interpolation from FACE solution to FACE cubature nodes is not required;
 *			2) Reduced aliasing because the normal numerical flux can be computed at the cubature nodes.
 *		Based on this discussion, it is unclear why this approach is not adopted instead of the traditional mortar
 *		method (Kopriva(1996)) as this alternative seems to satisfy both the conservation and outflow condition
 *		requirements which motivated the use of the mortar element method.
 *		If the adaptation method currently implemented is used, it GfS_fI (or L2fS_fI) should be changed to IfS_fI as
 *		there is no L2 projection currently being used. Note that IfS_fI == I_vS_vI of the FACE ELEMENT (ToBeDeleted).
 *		=> Perform convergence order verification for both methods and finalize conclusions (ToBeDeleted).
 *
 *		When adaptivity is used, exact flux evaluation on the FACEs cannot be used in 3D whenever TRIs are present
 *		because the cubature nodes do not form a basis for the mortar element. This implies that adaptivity necessarily
 *		introduces aliasing errors unless the WSH nodes are used for the FACE integration. Of course, aliasing will
 *		always be introduced for the lower order FACE because of the projection, but this seemed not to be significant
 *		based on previous results using PF = P+1 in the VOLUME and interpolating to integration nodes giving optimal
 *		convergence orders. (ToBeModified)
 *
 *		For WEDGE ELEMENTs, the nOrd arrays for QUAD FACEs are stored with the TRI OPs, while those for TRI FACEs are
 *		stored with the LINE OPs. While this is not logical, it precludes the need for an additional OP structure.
 *
 *		When adding in MPI functionality, will not have access to VIn/VOut => Store what is needed for each FACE in an
 *		MPI communication routine before running this function. (ToBeDeleted)
 *
 *		Vectorization is more involved for FACE terms as there are many more possible combinations than for VOLUMEs.
 *		Given that the VOLUME vectorization seems not to have a significant impact on performance, the FACE
 *		vectorization may not be pursued. (ToBeDeleted)
 *
 *	Notation:
 *
 *	References:
 *		Kopriva(1996)-A_Conservative_Staggered-Grid_Chebyshev_Multidomain_Method_for_Compressible_Flows_II._A_Semi-Structured_Method
 */

struct S_OPERATORS {
	unsigned int NvnS, NfnI, NfnS, NvnS_SF, NfnI_SF, NvnI_SF, *nOrdInOut, *nOrdOutIn;
	double       **ChiS_fI, **ChiS_vI, **I_Weak_FF, **I_Weak_VV, **ChiS_fS, **GfS_fI;

	struct S_OpCSR **ChiS_fI_sp, **I_Weak_FF_sp;
};

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const struct S_FACE *FACE,
                     const unsigned int IndClass);
static void compute_FACE_RHS_EFE    (void);
static void compute_FACE_RHS        (void);
//static void compute_FACEVec_RHS_EFE (void);
//static void compute_FACEVec_RHS     (void);

void explicit_FACE_info(void)
{
	// Initialize DB Parameters
	unsigned int Vectorized = DB.Vectorized,
	             Adapt      = DB.Adapt;

	switch (Adapt) {
case ADAPT_P: // ToBeModified (Also change setup_normals and output_to_paraview)
case ADAPT_H:
case ADAPT_HP:
	case ADAPT_0:
		switch (Vectorized) {
		case 0:
			compute_FACE_RHS_EFE();
			break;
		default:
			compute_FACE_RHS_EFE();
//			compute_FACEVec_RHS_EFE();
			break;
		}
		break;
	default: // ADAPT_P, ADAPT_H, ADAPT_HP
printf("Error: Should not be entering default in explicit_FACE_info.\n"), EXIT_MSG;
		switch (Vectorized) {
		case 0:
			compute_FACE_RHS();
			/* Difference: Interpolate to FACE basis, evaluate numerical flux, interpolate the cubature nodes on each
			 *             side. Requires n_Sf (Solution face) for curved elements. Make sure that L2 projections
			 *             are used.
			 * Note:       On TP FACEs, this routine should give identical results to EFE for conforming meshes.
			 */
			break;
		default:
			compute_FACE_RHS();
//			compute_FACEVec_RHS();
			break;
		}
		break;
	}
}

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const struct S_FACE *FACE,
                     const unsigned int IndClass)
{
	// Initialize DB Parameters
	unsigned int Adapt    = DB.Adapt,
	             ***SF_BE = DB.SF_BE;

	// Standard datatypes
	unsigned int PV, PF, Vtype, Eclass, FtypeInt, IndOrdInOut, IndOrdOutIn;

	struct S_ELEMENT *ELEMENT, *ELEMENT_OPS, *ELEMENT_FACE;

	// silence
	ELEMENT_OPS = NULL;

	PV       = VOLUME->P;
	PF       = FACE->P;
	Vtype    = VOLUME->type;
	Eclass   = VOLUME->Eclass;

	FtypeInt    = FACE->typeInt;
	IndOrdInOut = FACE->IndOrdInOut;
	IndOrdOutIn = FACE->IndOrdOutIn;

	ELEMENT       = get_ELEMENT_type(Vtype);
	ELEMENT_FACE = get_ELEMENT_FACE(Vtype,IndClass);
	if ((Eclass == C_TP && SF_BE[PF][0][1]) || (Eclass == C_WEDGE && SF_BE[PF][1][1]))
		ELEMENT_OPS = ELEMENT->ELEMENTclass[IndClass];
	else
		ELEMENT_OPS = ELEMENT;

	OPS->NvnS    = ELEMENT->NvnS[PV];
	OPS->NvnS_SF = ELEMENT_OPS->NvnS[PV];

	OPS->NfnS      = ELEMENT_FACE->NvnS[PF];
	OPS->ChiS_fS   = ELEMENT->ChiS_fS[PV][PF];
	if (FtypeInt == 's') {
		// Straight FACE Integration
		OPS->NfnI    = ELEMENT->NfnIs[PF][IndClass];
		OPS->NfnI_SF = ELEMENT_OPS->NfnIs[PF][0];
		OPS->NvnI_SF = ELEMENT_OPS->NvnIs[PF];

		OPS->ChiS_fI   = ELEMENT_OPS->ChiS_fIs[PV][PF];
		OPS->ChiS_vI   = ELEMENT_OPS->ChiS_vIs[PV][PF];
		OPS->I_Weak_FF = ELEMENT_OPS->Is_Weak_FF[PV][PF];
		OPS->I_Weak_VV = ELEMENT_OPS->Is_Weak_VV[PV][PF];

		OPS->ChiS_fI_sp   = ELEMENT->ChiS_fIs_sp[PV][PF];
		OPS->I_Weak_FF_sp = ELEMENT->Is_Weak_FF_sp[PV][PF];

		OPS->GfS_fI    = ELEMENT->GfS_fIs[PF][PF];

		OPS->nOrdInOut = ELEMENT_FACE->nOrd_fIs[PF][IndOrdInOut];
		switch (Adapt) {
		default: // ADAPT_P, ADAPT_H, ADAPT_HP
printf("Error: Should not be entering default in explicit_FACE_info.\n"), EXIT_MSG;
			OPS->nOrdOutIn = ELEMENT_FACE->nOrd_fS[PF][IndOrdOutIn];
			break;
case ADAPT_P: // ToBeModified (Also change setup_normals and output_to_paraview)
case ADAPT_H:
case ADAPT_HP:
		case ADAPT_0:
			OPS->nOrdOutIn = ELEMENT_FACE->nOrd_fIs[PF][IndOrdOutIn];
			break;
		}
	} else {
		// Curved FACE Integration
		OPS->NfnI    = ELEMENT->NfnIc[PF][IndClass];
		OPS->NfnI_SF = ELEMENT_OPS->NfnIc[PF][0];
		OPS->NvnI_SF = ELEMENT_OPS->NvnIc[PF];

		OPS->ChiS_fI   = ELEMENT_OPS->ChiS_fIc[PV][PF];
		OPS->ChiS_vI   = ELEMENT_OPS->ChiS_vIc[PV][PF];
		OPS->I_Weak_FF = ELEMENT_OPS->Ic_Weak_FF[PV][PF];
		OPS->I_Weak_VV = ELEMENT_OPS->Ic_Weak_VV[PV][PF];

		OPS->ChiS_fI_sp   = ELEMENT->ChiS_fIc_sp[PV][PF];
		OPS->I_Weak_FF_sp = ELEMENT->Ic_Weak_FF_sp[PV][PF];

		OPS->GfS_fI    = ELEMENT->GfS_fIc[PF][PF];

		OPS->nOrdInOut = ELEMENT_FACE->nOrd_fIc[PF][IndOrdInOut];
		switch (Adapt) {
		default: // ADAPT_P, ADAPT_H, ADAPT_HP
printf("Error: Should not be entering default in explicit_FACE_info.\n"), EXIT_MSG;
			OPS->nOrdOutIn = ELEMENT_FACE->nOrd_fS[PF][IndOrdOutIn];
			break;
case ADAPT_P: // ToBeModified (Also change setup_normals and output_to_paraview)
case ADAPT_H:
case ADAPT_HP:
		case ADAPT_0:
			OPS->nOrdOutIn = ELEMENT_FACE->nOrd_fIc[PF][IndOrdOutIn];
			break;
		}
	}
}

static void compute_FACE_RHS_EFE(void)
{
	// Initialize DB Parameters
	char         *Form            = DB.Form;
	unsigned int d                = DB.d,
	             NfrefMax         = DB.NfrefMax,
	             Nvar             = DB.Nvar,
	             Neq              = DB.Neq,
	             InviscidFluxType = DB.InviscidFluxType,
	             Collocated       = DB.Collocated,
	             *VFPartUnity     = DB.VFPartUnity,
	             ***SF_BE         = DB.SF_BE;

	// Standard datatypes
	unsigned int i, j, iInd, dim, P,
	             VfIn, VfOut, fIn, fOut, EclassIn, EclassOut, IndFType, Boundary, BC, BC_trail, SpOpIn, SpOpOut,
	             RowInd, RowSub, ReOrder, *RowTracker,
	             NfnI, NvnSIn, NvnSOut, *nOrdOutIn, *nOrdInOut,
	             NIn[3], NOut[3], Diag[3], NOut0, NOut1, NIn0, NIn1;
	double       *WIn_fI, *WOut_fI, *WOut_fIIn, *nFluxNum_fI, *RHSIn, *RHSOut, *n_fI, *detJF_fI, *OP[3], **OPF0, **OPF1;

	struct S_OPERATORS *OPSIn[2], *OPSOut[2];
	struct S_VOLUME    *VIn, *VOut;
	struct S_FACE     *FACE;

	// silence
	WOut_fI   = NULL;
	WOut_fIIn = NULL;
	NvnSOut   = 0;

	for (i = 0; i < 2; i++) {
		OPSIn[i]  = malloc(sizeof *OPSIn[i]);  // free
		OPSOut[i] = malloc(sizeof *OPSOut[i]); // free
	}

	for (FACE = DB.FACE; FACE; FACE = FACE->next) {
		P = FACE->P;

		// Obtain operators
		VIn    = FACE->VIn;
		VfIn   = FACE->VfIn;
		fIn    = VfIn/NfrefMax;
		SpOpIn = Collocated && (VfIn % NFREFMAX == 0 && VIn->P == P);

		EclassIn = VIn->Eclass;
		IndFType = get_IndFType(EclassIn,fIn);
		init_ops(OPSIn[0],VIn,FACE,0);
		if (VIn->type == WEDGE || VIn->type == PYR)
			init_ops(OPSIn[1],VIn,FACE,1);

		VOut    = FACE->VOut;
		VfOut   = FACE->VfOut;
		fOut    = VfOut/NfrefMax;
		SpOpOut = Collocated && (VfOut % NFREFMAX == 0 && VOut->P == P);

		EclassOut = VOut->Eclass;
		init_ops(OPSOut[0],VOut,FACE,0);
		if (VOut->type == WEDGE || VOut->type == PYR)
			init_ops(OPSOut[1],VOut,FACE,1);

		BC = FACE->BC;
		Boundary = !((VIn->indexg != VOut->indexg) || (VIn->indexg == VOut->indexg && fIn != fOut));
		// The second condition is for periodic elements which are connected to themselves
		if (Boundary != FACE->Boundary)
			printf("Error: Incorrect Boundary flag.\n"), EXIT_MSG;
		// ToBeDeleted: Replace with Boundary = FACE->Boundary

		// Compute WIn_fI
		NfnI   = OPSIn[IndFType]->NfnI;
		NvnSIn = OPSIn[IndFType]->NvnS;

		WIn_fI = malloc(NfnI*Nvar * sizeof *WIn_fI); // free
		if (EclassIn == C_TP && SF_BE[P][0][1]) {
			get_sf_parametersF(OPSIn[0]->NvnS_SF,OPSIn[0]->NvnI_SF,OPSIn[0]->ChiS_vI,
			                   OPSIn[0]->NvnS_SF,OPSIn[0]->NfnI_SF,OPSIn[0]->ChiS_fI,NIn,NOut,OP,d,VfIn,C_TP);

			if (SpOpIn) {
				for (dim = 0; dim < d; dim++)
					Diag[dim] = 2;
				Diag[fIn/2] = 0;
			} else {
				for (dim = 0; dim < d; dim++)
					Diag[dim] = 0;
			}

			sf_apply_d(VIn->What,WIn_fI,NIn,NOut,Nvar,OP,Diag,d);
		} else if (EclassIn == C_WEDGE && SF_BE[P][1][1]) {
			if (fIn < 3) { OPF0  = OPSIn[0]->ChiS_fI, OPF1  = OPSIn[1]->ChiS_vI;
			               NOut0 = OPSIn[0]->NfnI_SF, NOut1 = OPSIn[1]->NvnI_SF;
			} else {       OPF0  = OPSIn[0]->ChiS_vI, OPF1  = OPSIn[1]->ChiS_fI;
			               NOut0 = OPSIn[0]->NvnI_SF, NOut1 = OPSIn[1]->NfnI_SF; }
			get_sf_parametersF(OPSIn[0]->NvnS_SF,NOut0,OPF0,OPSIn[1]->NvnS_SF,NOut1,OPF1,NIn,NOut,OP,d,VfIn,C_WEDGE);

			if (SpOpIn) {
				for (dim = 0; dim < d; dim++)
					Diag[dim] = 2;
				if (fIn < 3)
					Diag[0] = 0;
				else
					Diag[2] = 0;
			} else {
				for (dim = 0; dim < d; dim++)
					Diag[dim] = 0;
				Diag[1] = 2;
			}

			sf_apply_d(VIn->What,WIn_fI,NIn,NOut,Nvar,OP,Diag,d);
		} else if ((SpOpIn && (EclassIn == C_TP || EclassIn == C_WEDGE)) || (VFPartUnity[EclassIn])) {
			mm_CTN_CSR_d(NfnI,Nvar,NvnSIn,OPSIn[0]->ChiS_fI_sp[VfIn],VIn->What,WIn_fI);
		} else {
			mm_CTN_d(NfnI,Nvar,NvnSIn,OPSIn[0]->ChiS_fI[VfIn],VIn->What,WIn_fI);
		}

		// Compute WOut_fI (Taking BCs into account if applicable)
		n_fI     = FACE->n_fI;

		nOrdInOut = OPSIn[IndFType]->nOrdInOut;
		nOrdOutIn = OPSIn[IndFType]->nOrdOutIn;

		NvnSOut = OPSOut[0]->NvnS;
		WOut_fIIn = malloc(NfnI*Nvar * sizeof *WOut_fIIn); // free
		if (BC == 0 || (BC % BC_STEP_SC > 50)) { // Internal/Periodic FACE
			WOut_fI = malloc(NfnI*Nvar * sizeof *WOut_fI); // free
			if (EclassOut == C_TP && SF_BE[P][0][1]) {
				get_sf_parametersF(OPSOut[0]->NvnS_SF,OPSOut[0]->NvnI_SF,OPSOut[0]->ChiS_vI,
				                   OPSOut[0]->NvnS_SF,OPSOut[0]->NfnI_SF,OPSOut[0]->ChiS_fI,NIn,NOut,OP,d,VfOut,C_TP);

// Note: Needs modification for h-adaptation (ToBeDeleted)
				if (SpOpOut) {
					for (dim = 0; dim < d; dim++)
						Diag[dim] = 2;
					Diag[fOut/2] = 0;
				} else {
					for (dim = 0; dim < d; dim++)
						Diag[dim] = 0;
				}

				sf_apply_d(VOut->What,WOut_fI,NIn,NOut,Nvar,OP,Diag,d);
			} else if (EclassOut == C_WEDGE && SF_BE[P][1][1]) {
				if (fOut < 3) { OPF0  = OPSOut[0]->ChiS_fI, OPF1  = OPSOut[1]->ChiS_vI;
				                NOut0 = OPSOut[0]->NfnI_SF, NOut1 = OPSOut[1]->NvnI_SF;
				} else {        OPF0  = OPSOut[0]->ChiS_vI, OPF1  = OPSOut[1]->ChiS_fI;
				                NOut0 = OPSOut[0]->NvnI_SF, NOut1 = OPSOut[1]->NfnI_SF; }
				get_sf_parametersF(OPSOut[0]->NvnS_SF,NOut0,OPF0,OPSOut[1]->NvnS_SF,NOut1,OPF1,NIn,NOut,OP,d,VfOut,C_WEDGE);

				if (SpOpOut) {
					for (dim = 0; dim < d; dim++)
						Diag[dim] = 2;
					if (fOut < 3)
						Diag[0] = 0;
					else
						Diag[2] = 0;
				} else {
					for (dim = 0; dim < d; dim++)
						Diag[dim] = 0;
					Diag[1] = 2;
				}

				sf_apply_d(VOut->What,WOut_fI,NIn,NOut,Nvar,OP,Diag,d);
			} else if ((SpOpOut && (EclassOut == C_TP || EclassOut == C_WEDGE)) || (VFPartUnity[EclassOut])) {
				mm_CTN_CSR_d(NfnI,Nvar,NvnSOut,OPSOut[0]->ChiS_fI_sp[VfOut],VOut->What,WOut_fI);
			} else {
				mm_CTN_d(NfnI,Nvar,NvnSOut,OPSOut[0]->ChiS_fI[VfOut],VOut->What,WOut_fI);
			}

			// Reorder WOut_fI to correspond to WIn_fI
			for (i = 0; i < Nvar; i++) {
				iInd = i*NfnI;
				for (j = 0; j < NfnI; j++) {
					BC_trail = BC % BC_STEP_SC;
					if (BC_trail > 50 || BC_trail == 0)
						WOut_fIIn[iInd+j] = WOut_fI[iInd+nOrdOutIn[j]];
					else
						WOut_fIIn[iInd+j] = WOut_fI[iInd+j];
				}
			}
			free(WOut_fI);
		} else { // Boundary FACE
			if (BC % BC_STEP_SC == BC_RIEMANN) {
				boundary_Riemann(NfnI,1,FACE->XYZ_fI,WIn_fI,NULL,WOut_fIIn,n_fI,d);
			} else if (BC % BC_STEP_SC == BC_SLIPWALL) {
				boundary_SlipWall(NfnI,1,WIn_fI,WOut_fIIn,n_fI,d);
			} else if (BC % BC_STEP_SC == BC_BACKPRESSURE) {
				boundary_BackPressure(NfnI,1,WIn_fI,WOut_fIIn,n_fI,d,Nvar);
			} else if (BC % BC_STEP_SC == BC_TOTAL_TP) {
				boundary_Total_TP(NfnI,1,FACE->XYZ_fI,WIn_fI,WOut_fIIn,n_fI,d,Nvar);
			} else if (BC % BC_STEP_SC == BC_SUPERSONIC_IN) {
				boundary_SupersonicInflow(NfnI,1,FACE->XYZ_fI,WIn_fI,WOut_fIIn,n_fI,d,Nvar);
			} else if (BC % BC_STEP_SC == BC_SUPERSONIC_OUT) {
				boundary_SupersonicOutflow(NfnI,1,FACE->XYZ_fI,WIn_fI,WOut_fIIn,n_fI,d,Nvar);
			} else {
				printf("Error: Unsupported (%d).\n",BC), EXIT_MSG;
			}
		}

		// Compute numerical flux
		nFluxNum_fI = malloc(NfnI*Neq * sizeof *nFluxNum_fI); // free
		detJF_fI = FACE->detJF_fI;

		switch (InviscidFluxType) {
		case FLUX_LF:
			flux_LF(NfnI,1,WIn_fI,WOut_fIIn,nFluxNum_fI,n_fI,d,Neq);
			break;
		case FLUX_ROE:
			flux_Roe(NfnI,1,WIn_fI,WOut_fIIn,nFluxNum_fI,n_fI,d,Neq);
			break;
		default:
			printf("Error: Unsupported InviscidFluxType used in explicit_FACE_info.\n"), EXIT_MSG;
			break;
		}

		// Multiply n dot FNum by the area element
		for (i = 0; i < Neq; i++) {
			iInd = i*NfnI;
			for (j = 0; j < NfnI; j++)
				nFluxNum_fI[iInd+j] *= detJF_fI[j];
		}

		// Compute FACE RHS terms
		RHSIn  = calloc(NvnSIn*Neq  , sizeof *RHSIn);  // keep (requires external free)
		RHSOut = calloc(NvnSOut*Neq , sizeof *RHSOut); // keep (requires external free)
		FACE->RHSIn  = RHSIn;
		FACE->RHSOut = RHSOut;

		RowTracker = malloc(NfnI * sizeof *RowTracker); // free

		if (strstr(Form,"Weak")) {
			// Interior FACE
			if (EclassIn == C_TP && SF_BE[P][0][1]) {
				get_sf_parametersF(OPSIn[0]->NvnI_SF,OPSIn[0]->NvnS_SF,OPSIn[0]->I_Weak_VV,
				                   OPSIn[0]->NfnI_SF,OPSIn[0]->NvnS_SF,OPSIn[0]->I_Weak_FF,NIn,NOut,OP,d,VfIn,C_TP);

				if (SpOpIn) {
					for (dim = 0; dim < d; dim++)
						Diag[dim] = 2;
					Diag[fIn/2] = 0;
				} else {
					for (dim = 0; dim < d; dim++)
						Diag[dim] = 0;
				}

				sf_apply_d(nFluxNum_fI,RHSIn,NIn,NOut,Neq,OP,Diag,d);
			} else if (EclassIn == C_WEDGE && SF_BE[P][1][1]) {
				if (fIn < 3) { OPF0 = OPSIn[0]->I_Weak_FF, OPF1 = OPSIn[1]->I_Weak_VV;
				               NIn0 = OPSIn[0]->NfnI_SF,   NIn1 = OPSIn[1]->NvnI_SF;
				} else {       OPF0 = OPSIn[0]->I_Weak_VV, OPF1 = OPSIn[1]->I_Weak_FF;
				               NIn0 = OPSIn[0]->NvnI_SF,   NIn1 = OPSIn[1]->NfnI_SF; }
				get_sf_parametersF(NIn0,OPSIn[0]->NvnS_SF,OPF0,NIn1,OPSIn[1]->NvnS_SF,OPF1,NIn,NOut,OP,d,VfIn,C_WEDGE);

				if (SpOpIn) {
					for (dim = 0; dim < d; dim++)
						Diag[dim] = 2;
					if (fIn < 3)
						Diag[0] = 0;
					else
						Diag[2] = 0;
				} else {
					for (dim = 0; dim < d; dim++)
						Diag[dim] = 0;
					Diag[1] = 2;
				}

				sf_apply_d(nFluxNum_fI,RHSIn,NIn,NOut,Neq,OP,Diag,d);
			} else if ((SpOpIn && (EclassIn == C_TP || EclassIn == C_WEDGE)) || (VFPartUnity[EclassIn])) {
				mm_CTN_CSR_d(NvnSIn,Neq,NfnI,OPSIn[0]->I_Weak_FF_sp[VfIn],nFluxNum_fI,RHSIn);
			} else  {
				mm_CTN_d(NvnSIn,Neq,NfnI,OPSIn[0]->I_Weak_FF[VfIn],nFluxNum_fI,RHSIn);
			}

			// Exterior FACE
			if (!Boundary) {
				// Use -ve normal for opposite FACE
				for (i = 0; i < Neq; i++) {
					iInd = i*NfnI;
					for (j = 0; j < NfnI; j++)
						nFluxNum_fI[iInd+j] *= -1.0;
				}

				// Re-arrange nFluxNum to match node ordering from opposite VOLUME
				for (i = 0; i < NfnI; i++)
					RowTracker[i] = i;

				for (RowInd = 0; RowInd < NfnI; RowInd++) {
					ReOrder = nOrdInOut[RowInd];
					for (RowSub = ReOrder; RowTracker[RowSub] != ReOrder; RowSub = RowTracker[RowSub])
						;

					if (RowInd != RowSub) {
						array_swap_d(&nFluxNum_fI[RowInd],&nFluxNum_fI[RowSub],Neq,NfnI);
						array_swap_ui(&RowTracker[RowInd],&RowTracker[RowSub],1,1);
					}
				}

				if (EclassOut == C_TP && SF_BE[P][0][1]) {
					get_sf_parametersF(OPSOut[0]->NvnI_SF,OPSOut[0]->NvnS_SF,OPSOut[0]->I_Weak_VV,
					                   OPSOut[0]->NfnI_SF,OPSOut[0]->NvnS_SF,OPSOut[0]->I_Weak_FF,NIn,NOut,OP,d,VfOut,C_TP);

					if (SpOpOut) {
						for (dim = 0; dim < d; dim++)
							Diag[dim] = 2;
						Diag[fOut/2] = 0;
					} else {
						for (dim = 0; dim < d; dim++)
							Diag[dim] = 0;
					}

					sf_apply_d(nFluxNum_fI,RHSOut,NIn,NOut,Neq,OP,Diag,d);
				} else if (EclassOut == C_WEDGE && SF_BE[P][1][1]) {
					if (fOut < 3) { OPF0 = OPSOut[0]->I_Weak_FF, OPF1 = OPSOut[1]->I_Weak_VV;
					                NIn0 = OPSOut[0]->NfnI_SF,   NIn1 = OPSOut[1]->NvnI_SF;
					} else {        OPF0 = OPSOut[0]->I_Weak_VV, OPF1 = OPSOut[1]->I_Weak_FF;
					                NIn0 = OPSOut[0]->NvnI_SF,   NIn1 = OPSOut[1]->NfnI_SF; }
					get_sf_parametersF(NIn0,OPSOut[0]->NvnS_SF,OPF0,NIn1,OPSOut[1]->NvnS_SF,OPF1,NIn,NOut,OP,d,VfOut,C_WEDGE);

					if (SpOpOut) {
						for (dim = 0; dim < d; dim++)
							Diag[dim] = 2;
						if (fOut < 3)
							Diag[0] = 0;
						else
							Diag[2] = 0;
					} else {
						for (dim = 0; dim < d; dim++)
							Diag[dim] = 0;
						Diag[1] = 2;
					}

					sf_apply_d(nFluxNum_fI,RHSOut,NIn,NOut,Neq,OP,Diag,d);
				} else if ((SpOpOut && (EclassOut == C_TP || EclassOut == C_WEDGE)) || (VFPartUnity[EclassOut])) {
					mm_CTN_CSR_d(NvnSOut,Neq,NfnI,OPSOut[0]->I_Weak_FF_sp[VfOut],nFluxNum_fI,RHSOut);
				} else {
					mm_CTN_d(NvnSOut,Neq,NfnI,OPSOut[0]->I_Weak_FF[VfOut],nFluxNum_fI,RHSOut);
				}
			}
		} else if (strstr(Form,"Strong")) {
			printf("Exiting: Implement the strong form in compute_FACE_RHS_EFE.\n"), EXIT_MSG;
		}

		free(RowTracker);
		free(WIn_fI);
		free(WOut_fIIn);
		free(nFluxNum_fI);
	}

	for (i = 0; i < 2; i++) {
		free(OPSIn[i]);
		free(OPSOut[i]);
	}
}

//static void compute_FACEVec_RHS_EFE(void) { }

static void compute_FACE_RHS(void)
{
	/*
	 *	Comments:
	 *		The part of this function to compute RHS terms is very similar (if not exactly the same) as that for
	 *		compute_FACE_RHS_EFE. Potentially combine these two functions to reduce redundant code (ToBeModified).
	 *		Possibly include Sum Factorization here after this is verified.
	 *		Think about whether there is any operator sparsity in this case.
	 */

	// Initialize DB Parameters
	char         *Form            = DB.Form;
	unsigned int d                = DB.d,
	             NfrefMax         = DB.NfrefMax,
	             Nvar             = DB.Nvar,
	             Neq              = DB.Neq,
	             Adapt            = DB.Adapt,
	             InviscidFluxType = DB.InviscidFluxType,
	             Collocated       = DB.Collocated,
	             *VFPartUnity     = DB.VFPartUnity,
	             ***SF_BE         = DB.SF_BE;

	// Standard datatypes
	unsigned int i, j, iInd, dim, P,
	             VfIn, VfOut, fIn, fOut, EclassIn, EclassOut, IndFType, Boundary, BC, BC_trail,
	             RowInd, RowSub, ReOrder, *RowTracker,
	             NfnS, NfnI, NvnSIn, NvnSOut, *nOrdOutIn, *nOrdInOut,
	             NIn[3], NOut[3], Diag[3], NIn0, NIn1;
	double       *WIn_fS, *WOut_fS, *WOut_fSIn, *nFluxNum_fS, *nFluxNum_fI, *RHSIn, *RHSOut, *n_fS, *detJF_fS,
	             *OP[3], **OPF0, **OPF1;

	struct S_OPERATORS *OPSIn[2], *OPSOut[2];
	struct S_VOLUME    *VIn, *VOut;
	struct S_FACE     *FACE;

	// silence
	nFluxNum_fI = NULL;

	for (i = 0; i < 2; i++) {
		OPSIn[i]  = malloc(sizeof *OPSIn[i]);  // free
		OPSOut[i] = malloc(sizeof *OPSOut[i]); // free
	}

	for (FACE = DB.FACE; FACE; FACE = FACE->next) {
		P = FACE->P;

		// Obtain operators
		VIn  = FACE->VIn;
		VfIn = FACE->VfIn;
		fIn  = VfIn/NfrefMax;

		EclassIn = VIn->Eclass;
		IndFType = get_IndFType(EclassIn,fIn);
		init_ops(OPSIn[0],VIn,FACE,0);
		if (VIn->type == WEDGE || VIn->type == PYR)
			init_ops(OPSIn[1],VIn,FACE,1);

		VOut  = FACE->VOut;
		VfOut = FACE->VfOut;
		fOut  = VfOut/NfrefMax;

		EclassOut = VOut->Eclass;
		init_ops(OPSOut[0],VOut,FACE,0);
		if (VOut->type == WEDGE || VOut->type == PYR)
			init_ops(OPSOut[1],VOut,FACE,1);

		BC = FACE->BC;
		Boundary = !((VIn->indexg != VOut->indexg) || (VIn->indexg == VOut->indexg && fIn != fOut));
		// The second condition is for periodic elements which are connected to themselves

		// Compute WIn_fS
		NfnS = OPSIn[IndFType]->NfnS;
		NvnSIn = OPSIn[IndFType]->NvnS;

		WIn_fS = malloc(NfnS*Nvar * sizeof *WIn_fS); // free
		// If VFPartUnity == 1, this operator would be sparse.
		mm_CTN_d(NfnS,Nvar,NvnSIn,OPSIn[0]->ChiS_fS[VfIn],VIn->What,WIn_fS);

		// Compute WOut_fI (Taking BCs into account if applicable)
		n_fS     = FACE->n_fS;

		nOrdInOut = OPSIn[IndFType]->nOrdInOut;
		nOrdOutIn = OPSIn[IndFType]->nOrdOutIn;

		NvnSOut = OPSOut[0]->NvnS;
		WOut_fSIn = malloc(NfnS*Nvar * sizeof *WOut_fSIn); // free
		if (BC == 0 || (BC % BC_STEP_SC > 50)) { // Internal/Periodic FACE
			WOut_fS = malloc(NfnS*Nvar * sizeof *WOut_fS); // free
			mm_CTN_d(NfnS,Nvar,NvnSOut,OPSOut[0]->ChiS_fS[VfOut],VOut->What,WOut_fS);

			// Reorder WOut_fS to correspond to WIn_fS
			for (i = 0; i < Nvar; i++) {
				iInd = i*NfnS;
				for (j = 0; j < NfnS; j++) {
					BC_trail = BC % BC_STEP_SC;
					if (BC_trail > 50 || BC_trail == 0)
						WOut_fSIn[iInd+j] = WOut_fS[iInd+nOrdOutIn[j]];
					else
						WOut_fSIn[iInd+j] = WOut_fS[iInd+j];
				}
			}
			free(WOut_fS);
		} else { // Boundary FACE
			if (BC % BC_STEP_SC == BC_RIEMANN) {
				boundary_Riemann(NfnS,1,FACE->XYZ_fS,WIn_fS,NULL,WOut_fSIn,n_fS,d);
			} else if (BC % BC_STEP_SC == BC_SLIPWALL) {
				boundary_SlipWall(NfnS,1,WIn_fS,WOut_fSIn,n_fS,d);
			} else {
				printf("Error: Unsupported BC in explicit_FACE_info.\n"), EXIT_MSG;
			}
		}
/*
if (FACE->indexg == 240) {
printf("%d\n",FACE->indexg);
//array_print_d(NfnS,Nvar,WIn_fS,'C');
//array_print_d(NfnS,Nvar,WOut_fSIn,'C');
}
*/
		// Compute numerical flux
		nFluxNum_fS = malloc(NfnS*Neq * sizeof *nFluxNum_fS); // free
		detJF_fS = FACE->detJF_fS;

		switch (InviscidFluxType) {
		case FLUX_LF:
			flux_LF(NfnS,1,WIn_fS,WOut_fSIn,nFluxNum_fS,n_fS,d,Neq);
			break;
		case FLUX_ROE:
			flux_Roe(NfnS,1,WIn_fS,WOut_fSIn,nFluxNum_fS,n_fS,d,Neq);
			break;
		default:
			printf("Error: Unsupported InviscidFluxType used in explicit_FACE_info.\n"), EXIT_MSG;
			break;
		}
		free(WIn_fS);
		free(WOut_fSIn);
/*
array_print_d(NfnS,d,n_fS,'R');
array_print_d(NfnS,Neq,nFluxNum_fS,'C');
*/
		// Multiply n dot FNum by the area element
		// Potentially do this after projection? Test and see. ToBeDeleted
		for (i = 0; i < Neq; i++) {
			iInd = i*NfnS;
			for (j = 0; j < NfnS; j++)
				nFluxNum_fS[iInd+j] *= detJF_fS[j];
		}

		// L2 projection to FACE cubature nodes
		NfnI   = OPSIn[IndFType]->NfnI;

// Don't forget to reorder nF_I for Outer FACE

		switch (Adapt) {
		default: // ADAPT_HP
			break;
		case ADAPT_P:
		case ADAPT_H:
			nFluxNum_fI = malloc(NfnI*Neq * sizeof *nFluxNum_fI); // free
			mm_CTN_d(NfnI,Neq,NfnS,OPSIn[0]->GfS_fI[VfIn],nFluxNum_fS,nFluxNum_fI);

			break;
		}
		free(nFluxNum_fS);

		RowTracker = malloc(NfnI * sizeof *RowTracker); // free


/*
if (FACE->indexg == 240) {
array_print_d(NfnI,Neq,nFluxNum_fI,'C');
//exit(1);
}
*/
		// Compute FACE RHS terms
		RHSIn  = calloc(NvnSIn*Neq  , sizeof *RHSIn);  // keep (requires external free)
		RHSOut = calloc(NvnSOut*Neq , sizeof *RHSOut); // keep (requires external free)
		FACE->RHSIn  = RHSIn;
		FACE->RHSOut = RHSOut;

		if (strstr(Form,"Weak")) {
			// Interior FACE
			if (EclassIn == C_TP && SF_BE[P][0][1]) {
				get_sf_parametersF(OPSIn[0]->NvnI_SF,OPSIn[0]->NvnS_SF,OPSIn[0]->I_Weak_VV,
								   OPSIn[0]->NfnI_SF,OPSIn[0]->NvnS_SF,OPSIn[0]->I_Weak_FF,NIn,NOut,OP,d,VfIn,C_TP);

				for (dim = 0; dim < d; dim++)
					Diag[dim] = 0;

				sf_apply_d(nFluxNum_fI,RHSIn,NIn,NOut,Neq,OP,Diag,d);
			} else if (EclassIn == C_WEDGE && SF_BE[P][1][1]) {
				if (fIn < 3) { OPF0 = OPSIn[0]->I_Weak_FF, OPF1 = OPSIn[1]->I_Weak_VV;
							   NIn0 = OPSIn[0]->NfnI_SF,   NIn1 = OPSIn[1]->NvnI_SF;
				} else {       OPF0 = OPSIn[0]->I_Weak_VV, OPF1 = OPSIn[1]->I_Weak_FF;
							   NIn0 = OPSIn[0]->NvnI_SF,   NIn1 = OPSIn[1]->NfnI_SF; }
				get_sf_parametersF(NIn0,OPSIn[0]->NvnS_SF,OPF0,NIn1,OPSIn[1]->NvnS_SF,OPF1,NIn,NOut,OP,d,VfIn,C_WEDGE);

				for (dim = 0; dim < d; dim++)
					Diag[dim] = 0;
				Diag[1] = 2;

				sf_apply_d(nFluxNum_fI,RHSIn,NIn,NOut,Neq,OP,Diag,d);
			} else if ((Collocated && (EclassIn == C_TP || EclassIn == C_WEDGE)) || (VFPartUnity[EclassIn])) {
				mm_CTN_CSR_d(NvnSIn,Neq,NfnI,OPSIn[0]->I_Weak_FF_sp[VfIn],nFluxNum_fI,RHSIn);
			} else  {
				mm_CTN_d(NvnSIn,Neq,NfnI,OPSIn[0]->I_Weak_FF[VfIn],nFluxNum_fI,RHSIn);
			}

			// Exterior FACE
			if (!Boundary) {
				// Use -ve normal for opposite FACE
				for (i = 0; i < Neq; i++) {
					iInd = i*NfnI;
					for (j = 0; j < NfnI; j++)
						nFluxNum_fI[iInd+j] *= -1.0;
				}

				// Re-arrange nFluxNum to match node ordering from opposite VOLUME
				for (i = 0; i < NfnI; i++)
					RowTracker[i] = i;

				for (RowInd = 0; RowInd < NfnI; RowInd++) {
					ReOrder = nOrdInOut[RowInd];
					for (RowSub = ReOrder; RowTracker[RowSub] != ReOrder; RowSub = RowTracker[RowSub])
						;

					if (RowInd != RowSub) {
						array_swap_d(&nFluxNum_fI[RowInd],&nFluxNum_fI[RowSub],Neq,NfnI);
						array_swap_ui(&RowTracker[RowInd],&RowTracker[RowSub],1,1);
					}
				}

				if (EclassOut == C_TP && SF_BE[P][0][1]) {
					get_sf_parametersF(OPSOut[0]->NvnI_SF,OPSOut[0]->NvnS_SF,OPSOut[0]->I_Weak_VV,
									   OPSOut[0]->NfnI_SF,OPSOut[0]->NvnS_SF,OPSOut[0]->I_Weak_FF,NIn,NOut,OP,d,VfOut,C_TP);

					for (dim = 0; dim < d; dim++)
						Diag[dim] = 0;

					sf_apply_d(nFluxNum_fI,RHSOut,NIn,NOut,Neq,OP,Diag,d);
				} else if (EclassOut == C_WEDGE && SF_BE[P][1][1]) {
					if (fOut < 3) { OPF0 = OPSOut[0]->I_Weak_FF, OPF1 = OPSOut[1]->I_Weak_VV;
					                NIn0 = OPSOut[0]->NfnI_SF,   NIn1 = OPSOut[1]->NvnI_SF;
					} else {        OPF0 = OPSOut[0]->I_Weak_VV, OPF1 = OPSOut[1]->I_Weak_FF;
					                NIn0 = OPSOut[0]->NvnI_SF,   NIn1 = OPSOut[1]->NfnI_SF; }
					get_sf_parametersF(NIn0,OPSOut[0]->NvnS_SF,OPF0,NIn1,OPSOut[1]->NvnS_SF,OPF1,NIn,NOut,OP,d,VfOut,C_WEDGE);

					for (dim = 0; dim < d; dim++)
						Diag[dim] = 0;
					Diag[1] = 2;

					sf_apply_d(nFluxNum_fI,RHSOut,NIn,NOut,Neq,OP,Diag,d);
				} else if ((Collocated && (EclassOut == C_TP || EclassOut == C_WEDGE)) || (VFPartUnity[EclassOut])) {
					mm_CTN_CSR_d(NvnSOut,Neq,NfnI,OPSOut[0]->I_Weak_FF_sp[VfOut],nFluxNum_fI,RHSOut);
				} else {
					mm_CTN_d(NvnSOut,Neq,NfnI,OPSOut[0]->I_Weak_FF[VfOut],nFluxNum_fI,RHSOut);
				}
			}
		} else if (strstr(Form,"Strong")) {
			printf("Exiting: Implement the strong form in compute_FACE_RHS_EFE.\n"), EXIT_MSG;
		}
		free(RowTracker);
		free(nFluxNum_fI);
//exit(1);
	}

	for (i = 0; i < 2; i++) {
		free(OPSIn[i]);
		free(OPSOut[i]);
	}
//exit(1);
}

//static void compute_FACEVec_RHS(void) { }
