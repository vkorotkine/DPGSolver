// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "implicit_FACE_info.h"

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

// Remove what is no longer needed (ToBeDeleted)
#include "solver_functions.h"

#include "element_functions.h"
#include "sum_factorization.h"
#include "matrix_functions.h"
#include "boundary_conditions.h"
#include "jacobian_boundary_conditions.h"
#include "fluxes_inviscid.h"
#include "jacobian_fluxes_inviscid.h"
#include "array_swap.h"


/*
 *	Purpose:
 *		Evaluate the FACE contributions to the LHS term.
 *
 *	Comments:
 *		Check for relevant comments in implicit_VOLUME_info.
 *		After profiling, decide whether it would be important to include the sparse operators here. (ToBeDeleted)
 *
 *	Notation:
 *
 *	References:
 */

static void compute_FACE_EFE (void);

void implicit_FACE_info(void)
{
	compute_FACE_EFE();
}

static void compute_FACE_EFE(void)
{
	// Initialize DB Parameters
	char         *Form            = DB.Form;
	unsigned int d                = DB.d,
	             NfrefMax         = DB.NfrefMax,
	             Nvar             = DB.Nvar,
	             Neq              = DB.Neq,
	             Collocated       = DB.Collocated,
	             *VFPartUnity     = DB.VFPartUnity,
	             ***SF_BE         = DB.SF_BE;

	// Standard datatypes
	unsigned int i, j, n, dim, P, eq, var, iMax,
	             InddnFdWIn, InddnFdWOut, IndnF, IndI, Indeqvar, IndLHS,
	             VfIn, VfOut, fIn, fOut, EclassIn, EclassOut, IndFType, Boundary, BC, SpOpIn, SpOpOut,
	             RowInd, RowSub, ReOrder, *RowTracker,
	             NfnI, NvnSIn, NvnSOut, *nOrdOutIn, *nOrdInOut,
	             NIn[3], NOut[3], Diag[3], NIn0, NIn1;
	double       *WIn_fI, *WOut_fIIn,
	             *nFluxNum_fI, *dnFluxNumdWIn_fI, *dnFluxNumdWOut_fI,
	             *RHSIn, *RHSOut, *LHSInIn, *LHSOutIn, *LHSInOut, *LHSOutOut, *n_fI, *detJF_fI, *I_FF, *IdnFdW,
	             *ChiS_fI, *ChiS_fIOutIn, *ChiS_fIInOut,
	             *OP[3], **OPF0, **OPF1;

	struct S_OPERATORS_F *OPSIn[2], *OPSOut[2];
	struct S_VOLUME    *VIn, *VOut;
	struct S_FACE     *FACE;

	struct S_FDATA *FDATAL, *FDATAR;
	FDATAL = malloc(sizeof *FDATAL); // free
	FDATAR = malloc(sizeof *FDATAR); // free

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
		n_fI = FACE->n_fI;

		nOrdInOut = OPSIn[IndFType]->nOrdInOut;
		nOrdOutIn = OPSIn[IndFType]->nOrdOutIn;

		NvnSOut = OPSOut[0]->NvnS;
		WOut_fIIn = malloc(NfnI*Nvar * sizeof *WOut_fIIn); // free
		compute_WR_fIL(FDATAR,WIn_fI,WOut_fIIn);

		// Compute numerical flux and its Jacobian
		nFluxNum_fI       = malloc(NfnI*Neq      * sizeof *nFluxNum_fI);       // free
		dnFluxNumdWIn_fI  = malloc(NfnI*Neq*Nvar * sizeof *dnFluxNumdWIn_fI);  // free
		dnFluxNumdWOut_fI = malloc(NfnI*Neq*Nvar * sizeof *dnFluxNumdWOut_fI); // free

		detJF_fI = FACE->detJF_fI;

		struct S_NumericalFlux *NFluxData = malloc(sizeof *NFluxData); // tbd
		NFluxData->WL_fIL = WIn_fI;
		NFluxData->WR_fIL = WOut_fIIn;
		NFluxData->nFluxNum_fIL = nFluxNum_fI;
		NFluxData->dnFluxNumdWL_fIL = dnFluxNumdWIn_fI;
		NFluxData->dnFluxNumdWR_fIL = dnFluxNumdWOut_fI;
		compute_numerical_flux(FDATAL,NFluxData,'I');
/*
		switch (InviscidFluxType) {
		case FLUX_LF:
			flux_LF(NfnI,1,WIn_fI,WOut_fIIn,nFluxNum_fI,n_fI,d,Neq);
			jacobian_flux_LF(NfnI,1,WIn_fI,WOut_fIIn,dnFluxNumdWIn_fI,n_fI,d,Neq,'L');
			jacobian_flux_LF(NfnI,1,WIn_fI,WOut_fIIn,dnFluxNumdWOut_fI,n_fI,d,Neq,'R');
			break;
		case FLUX_ROE:
			flux_Roe(NfnI,1,WIn_fI,WOut_fIIn,nFluxNum_fI,n_fI,d,Neq);
			jacobian_flux_Roe(NfnI,1,WIn_fI,WOut_fIIn,dnFluxNumdWIn_fI,n_fI,d,Neq,'L');
			jacobian_flux_Roe(NfnI,1,WIn_fI,WOut_fIIn,dnFluxNumdWOut_fI,n_fI,d,Neq,'R');
			break;
		default:
			EXIT_UNSUPPORTED;
			break;
		}

		// Include the BC information in dnFluxNumWIn_fI if on a boundary
		if (Boundary) {
			dWOutdWIn = malloc(NfnI*Nvar*Nvar * sizeof *dWOutdWIn); // free
			if (BC % BC_STEP_SC == BC_RIEMANN)
				jacobian_boundary_Riemann(NfnI,1,FACE->XYZ_fI,WIn_fI,NULL,dWOutdWIn,n_fI,d,Neq);
			else if (BC % BC_STEP_SC == BC_SLIPWALL)
				jacobian_boundary_SlipWall(NfnI,1,WIn_fI,dWOutdWIn,n_fI,d,Neq);
			else if (BC % BC_STEP_SC == BC_BACKPRESSURE)
				jacobian_boundary_BackPressure(NfnI,1,WIn_fI,dWOutdWIn,n_fI,d,Neq);
			else if (BC % BC_STEP_SC == BC_TOTAL_TP)
				jacobian_boundary_Total_TP(NfnI,1,FACE->XYZ_fI,WIn_fI,dWOutdWIn,n_fI,d,Neq);
			else if (BC % BC_STEP_SC == BC_SUPERSONIC_IN)
				jacobian_boundary_SupersonicInflow(NfnI,1,FACE->XYZ_fI,WIn_fI,dWOutdWIn,n_fI,d,Neq);
			else if (BC % BC_STEP_SC == BC_SUPERSONIC_OUT)
				jacobian_boundary_SupersonicOutflow(NfnI,1,FACE->XYZ_fI,WIn_fI,dWOutdWIn,n_fI,d,Neq);
			else
				printf("Error: Unsupported BC.\n"), EXIT_MSG;

			for (eq = 0; eq < Neq; eq++) {
			for (var = 0; var < Nvar; var++) {
				InddnFdWIn = (eq*Nvar+var)*NfnI;

				for (i = 0; i < Nvar; i++) {
					InddnFdWOut  = (eq*Neq+i)*NfnI;
					InddWOutdWIn = (var*Nvar+i)*NfnI;
					for (n = 0; n < NfnI; n++)
						dnFluxNumdWIn_fI[InddnFdWIn+n] += dnFluxNumdWOut_fI[InddnFdWOut+n]*dWOutdWIn[InddWOutdWIn+n];
				}
			}}

			free(dWOutdWIn);
		}
*/
		// Multiply nFNum and its Jacobian by the area element
		for (eq = 0; eq < Neq; eq++) {
			IndnF = eq*NfnI;
			for (n = 0; n < NfnI; n++)
				nFluxNum_fI[IndnF+n] *= detJF_fI[n];

			for (var = 0; var < Nvar; var++) {
				InddnFdWIn = (eq*Nvar+var)*NfnI;
				for (n = 0; n < NfnI; n++) {
					dnFluxNumdWIn_fI[InddnFdWIn+n]  *= detJF_fI[n];
					dnFluxNumdWOut_fI[InddnFdWIn+n] *= detJF_fI[n];
				}
			}
		}

		// Compute FACE RHS and LHS terms
		RHSIn     = calloc(NvnSIn*Neq               , sizeof *RHSIn);     // keep (requires external free)
		RHSOut    = calloc(NvnSOut*Neq              , sizeof *RHSOut);    // keep (requires external free)
		LHSInIn   = calloc(NvnSIn*NvnSIn*Neq*Nvar   , sizeof *LHSInIn);   // keep (requires external free)
		LHSOutIn  = calloc(NvnSIn*NvnSOut*Neq*Nvar  , sizeof *LHSOutIn);  // keep (requires external free)
		LHSInOut  = calloc(NvnSOut*NvnSIn*Neq*Nvar  , sizeof *LHSInOut);  // keep (requires external free)
		LHSOutOut = calloc(NvnSOut*NvnSOut*Neq*Nvar , sizeof *LHSOutOut); // keep (requires external free)

		FACE->RHSIn     = RHSIn;
		FACE->RHSOut    = RHSOut;
		FACE->LHSInIn   = LHSInIn;
		FACE->LHSOutIn  = LHSOutIn;
		FACE->LHSInOut  = LHSInOut;
		FACE->LHSOutOut = LHSOutOut;

		RowTracker = malloc(NfnI * sizeof *RowTracker); // free

		if (strstr(Form,"Weak")) {
			// Interior FACE

			// RHS
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

			// LHS
//			if ((SpOpIn && (EclassIn == C_TP || EclassIn == C_WEDGE)) || (VFPartUnity[EclassIn])) {
//			} else  {
				I_FF = OPSIn[0]->I_Weak_FF[VfIn];
				IdnFdW = malloc(NvnSIn*NfnI * sizeof *IdnFdW); // free
				for (eq = 0; eq < Neq; eq++) {
				for (var = 0; var < Nvar; var++) {
					Indeqvar = eq*Nvar+var;

					InddnFdWIn = Indeqvar*NfnI;
					for (i = 0; i < NvnSIn; i++) {
						IndI = i*NfnI;
						for (j = 0; j < NfnI; j++)
							IdnFdW[IndI+j] = I_FF[IndI+j]*dnFluxNumdWIn_fI[InddnFdWIn+j];
					}

					IndLHS = Indeqvar*NvnSIn*NvnSIn;
					mm_d(CBRM,CBNT,CBNT,NvnSIn,NvnSIn,NfnI,1.0,0.0,IdnFdW,OPSIn[0]->ChiS_fI[VfIn],&LHSInIn[IndLHS]);
				}}

				free(IdnFdW);
//			}

			// Exterior FACE
			if (!Boundary) {
				// RHS

				// Use -ve normal for opposite FACE
				for (i = 0, iMax = Neq*NfnI; i < iMax; i++)
					nFluxNum_fI[i] *= -1.0;

				// Rearrange nFluxNum to match node ordering from opposite VOLUME
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

				// LHS (Currently not using sum factorization)
				if (SF_BE[P][0][1] || SF_BE[P][1][1])
					// Add support for sum factorization below if desired or ensure that correct operators are being
					// used as set in init_ops.
					EXIT_UNSUPPORTED;

				// OutIn (Effect of Out on In)
				IdnFdW       = malloc(NvnSIn*NfnI  * sizeof *IdnFdW);       // free
				ChiS_fIOutIn = malloc(NvnSOut*NfnI * sizeof *ChiS_fIOutIn); // free

				I_FF = OPSIn[0]->I_Weak_FF[VfIn];

				ChiS_fI = OPSOut[0]->ChiS_fI[VfOut];
				for (i = 0; i < NfnI; i++) {
				for (j = 0; j < NvnSOut; j++) {
					ChiS_fIOutIn[i*NvnSOut+j] = ChiS_fI[nOrdOutIn[i]*NvnSOut+j];
				}}

				for (eq = 0; eq < Neq; eq++) {
				for (var = 0; var < Nvar; var++) {
					Indeqvar = eq*Nvar+var;

					InddnFdWOut = Indeqvar*NfnI;
					for (i = 0; i < NvnSIn; i++) {
						IndI = i*NfnI;
						for (j = 0; j < NfnI; j++)
							IdnFdW[IndI+j] = I_FF[IndI+j]*dnFluxNumdWOut_fI[InddnFdWOut+j];
					}

					IndLHS = Indeqvar*NvnSIn*NvnSOut;
					mm_d(CBRM,CBNT,CBNT,NvnSIn,NvnSOut,NfnI,1.0,0.0,IdnFdW,ChiS_fIOutIn,&LHSOutIn[IndLHS]);
				}}

				free(ChiS_fIOutIn);
				free(IdnFdW);


				// Use -ve normal for opposite FACE
				for (i = 0, iMax = Neq*Nvar*NfnI; i < iMax; i++) {
					dnFluxNumdWIn_fI[i] *= -1.0;
					dnFluxNumdWOut_fI[i] *= -1.0;
				}

				// Rearrange dnFluxNumdWIn/Out to match node ordering from opposite VOLUME
				for (i = 0; i < NfnI; i++)
					RowTracker[i] = i;

				for (RowInd = 0; RowInd < NfnI; RowInd++) {
					ReOrder = nOrdInOut[RowInd];
					for (RowSub = ReOrder; RowTracker[RowSub] != ReOrder; RowSub = RowTracker[RowSub])
						;

					if (RowInd != RowSub) {
						array_swap_d(&dnFluxNumdWIn_fI[RowInd], &dnFluxNumdWIn_fI[RowSub], Neq*Nvar,NfnI);
						array_swap_d(&dnFluxNumdWOut_fI[RowInd],&dnFluxNumdWOut_fI[RowSub],Neq*Nvar,NfnI);
						array_swap_ui(&RowTracker[RowInd],&RowTracker[RowSub],1,1);
					}
				}

				I_FF = OPSOut[0]->I_Weak_FF[VfOut];

				// InOut (Effect of In on Out)
				IdnFdW       = malloc(NvnSOut*NfnI * sizeof *IdnFdW);       // free
				ChiS_fIInOut = malloc(NvnSIn*NfnI  * sizeof *ChiS_fIInOut); // free

				ChiS_fI = OPSIn[0]->ChiS_fI[VfIn];
				for (i = 0; i < NfnI; i++) {
				for (j = 0; j < NvnSIn; j++) {
					ChiS_fIInOut[i*NvnSIn+j] = ChiS_fI[nOrdInOut[i]*NvnSIn+j];
				}}

				for (eq = 0; eq < Neq; eq++) {
				for (var = 0; var < Nvar; var++) {
					Indeqvar = eq*Nvar+var;

					InddnFdWIn  = Indeqvar*NfnI;
					for (i = 0; i < NvnSOut; i++) {
						IndI = i*NfnI;
						for (j = 0; j < NfnI; j++)
							IdnFdW[IndI+j] = I_FF[IndI+j]*dnFluxNumdWIn_fI[InddnFdWIn+j];
					}

					IndLHS = Indeqvar*NvnSOut*NvnSIn;
					mm_d(CBRM,CBNT,CBNT,NvnSOut,NvnSIn,NfnI,1.0,0.0,IdnFdW,ChiS_fIInOut,&LHSInOut[IndLHS]);
				}}

				free(ChiS_fIInOut);
				free(IdnFdW);

				// OutOut
				IdnFdW       = malloc(NvnSOut*NfnI * sizeof *IdnFdW);       // free

				for (eq = 0; eq < Neq; eq++) {
				for (var = 0; var < Nvar; var++) {
					Indeqvar = eq*Nvar+var;

					InddnFdWOut = Indeqvar*NfnI;
					for (i = 0; i < NvnSOut; i++) {
						IndI = i*NfnI;
						for (j = 0; j < NfnI; j++)
							IdnFdW[IndI+j] = I_FF[IndI+j]*dnFluxNumdWOut_fI[InddnFdWOut+j];
					}

					IndLHS = (eq*Nvar+var)*NvnSOut*NvnSOut;
					mm_d(CBRM,CBNT,CBNT,NvnSOut,NvnSOut,NfnI,1.0,0.0,IdnFdW,OPSOut[0]->ChiS_fI[VfOut],&LHSOutOut[IndLHS]);
				}}

				free(IdnFdW);
			}
		} else if (strstr(Form,"Strong")) {
			printf("Exiting: Implement the strong form in compute_FACE_EFE.\n"), exit(1);
		}

		free(RowTracker);
		free(WIn_fI);
		free(WOut_fIIn);
		free(nFluxNum_fI);
		free(dnFluxNumdWIn_fI);
		free(dnFluxNumdWOut_fI);
	}

	for (i = 0; i < 2; i++) {
		free(OPSIn[i]);
		free(OPSOut[i]);
	}
}
