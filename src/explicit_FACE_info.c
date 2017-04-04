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

#include "solver_functions.h"
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
 *
*	*** ToBeModified ***
 *
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

static void compute_FACE_RHS_EFE    (void);

void explicit_FACE_info(void)
{
	compute_FACE_RHS_EFE();
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
	             VfIn, VfOut, fIn, fOut, EclassIn, EclassOut, IndFType, Boundary, BC, SpOpIn, SpOpOut,
	             RowInd, RowSub, ReOrder, *RowTracker,
	             NfnI, NvnSIn, NvnSOut, *nOrdOutIn, *nOrdInOut,
	             NIn[3], NOut[3], Diag[3], NIn0, NIn1;
	double       *WIn_fI, *WOut_fI, *WOut_fIIn, *nFluxNum_fI, *RHSIn, *RHSOut, *n_fI, *detJF_fI, *OP[3], **OPF0, **OPF1;

	struct S_OPERATORS_F *OPSIn[2], *OPSOut[2];
	struct S_VOLUME    *VIn, *VOut;
	struct S_FACE     *FACE;

	struct S_FDATA *FDATAL, *FDATAR;
	FDATAL = malloc(sizeof *FDATAL); // free
	FDATAR = malloc(sizeof *FDATAR); // free

	// silence
	WOut_fI   = NULL;
	WOut_fIIn = NULL;
	NvnSOut   = 0;

	for (i = 0; i < 2; i++) {
		OPSIn[i]  = malloc(sizeof *OPSIn[i]);  // free
		OPSOut[i] = malloc(sizeof *OPSOut[i]); // free
	}
	for (FACE = DB.FACE; FACE; FACE = FACE->next) {
		FDATAL->OPS = OPSIn;
		FDATAR->OPS = OPSOut;

		init_FDATA(FDATAL,FACE,'L');
		init_FDATA(FDATAR,FACE,'R');

		P = FACE->P;

		// Obtain operators
		VIn    = FACE->VIn;
		VfIn   = FACE->VfIn;
		fIn    = VfIn/NfrefMax;
		SpOpIn = Collocated && (VfIn % NFREFMAX == 0 && VIn->P == P);

		EclassIn = VIn->Eclass;
		IndFType = get_IndFType(EclassIn,fIn);
		init_ops_FACE(OPSIn[0],VIn,FACE,0);
		if (VIn->type == WEDGE || VIn->type == PYR)
			init_ops_FACE(OPSIn[1],VIn,FACE,1);

		VOut    = FACE->VOut;
		VfOut   = FACE->VfOut;
		fOut    = VfOut/NfrefMax;
		SpOpOut = Collocated && (VfOut % NFREFMAX == 0 && VOut->P == P);

		EclassOut = VOut->Eclass;
		init_ops_FACE(OPSOut[0],VOut,FACE,0);
		if (VOut->type == WEDGE || VOut->type == PYR)
			init_ops_FACE(OPSOut[1],VOut,FACE,1);

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
		compute_W_fI(FDATAL,WIn_fI);

		// Compute WOut_fI (Taking BCs into account if applicable)
		n_fI     = FACE->n_fI;

		nOrdInOut = OPSIn[IndFType]->nOrdInOut;
		nOrdOutIn = OPSIn[IndFType]->nOrdOutIn;

		NvnSOut = OPSOut[0]->NvnS;
		WOut_fIIn = malloc(NfnI*Nvar * sizeof *WOut_fIIn); // free
		compute_WR_fIL(FDATAR,WIn_fI,WOut_fIIn);

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
