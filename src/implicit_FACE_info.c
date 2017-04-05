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
#include "array_print.h"


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
	unsigned int NfrefMax         = DB.NfrefMax,
	             Nvar             = DB.Nvar,
	             Neq              = DB.Neq,
	             Collocated       = DB.Collocated;

	// Standard datatypes
	unsigned int i, j, P, eq, var,
	             InddnFdWIn, InddnFdWOut, IndI, Indeqvar, IndLHS,
	             VfIn, VfOut, fIn, fOut, EclassIn, EclassOut, IndFType, Boundary, BC, SpOpIn, SpOpOut,
	             *RowTracker,
	             NfnI, NvnSIn, NvnSOut, *nOrdOutIn, *nOrdInOut;
	double       *WIn_fI, *WOut_fIIn,
	             *nFluxNum_fI, *dnFluxNumdWIn_fI, *dnFluxNumdWOut_fI,
	             *RHSIn, *RHSOut, *LHSInIn, *LHSOutIn, *LHSInOut, *LHSOutOut, *n_fI, *detJF_fI, *I_FF, *IdnFdW,
	             *ChiS_fI, *ChiS_fIOutIn, *ChiS_fIInOut;

	struct S_OPERATORS_F *OPSIn[2], *OPSOut[2];
	struct S_VOLUME    *VIn, *VOut;
	struct S_FACE     *FACE;

	struct S_FDATA *FDATAL, *FDATAR;
	FDATAL = malloc(sizeof *FDATAL); // free
	FDATAR = malloc(sizeof *FDATAR); // free

	struct S_NumericalFlux *NFluxData = malloc(sizeof *NFluxData); // free
	FDATAL->NFluxData = NFluxData;
	FDATAR->NFluxData = NFluxData;

	for (i = 0; i < 2; i++) {
		OPSIn[i]  = malloc(sizeof *OPSIn[i]);  // free
		OPSOut[i] = malloc(sizeof *OPSOut[i]); // free
	}

	for (FACE = DB.FACE; FACE; FACE = FACE->next) {
		FDATAL->OPS = OPSIn;
		FDATAR->OPS = OPSOut;

		init_FDATA(FDATAL,FACE,'L');
		init_FDATA(FDATAR,FACE,'R');


		// Compute WIn_fI
		NfnI   = OPSIn[IndFType]->NfnI;

		WIn_fI = malloc(NfnI*Nvar * sizeof *WIn_fI); // free
		compute_W_fI(FDATAL,WIn_fI);

		WOut_fIIn = malloc(NfnI*Nvar * sizeof *WOut_fIIn); // free
		compute_WR_fIL(FDATAR,WIn_fI,WOut_fIIn);

		// Compute numerical flux and its Jacobian
		nFluxNum_fI       = malloc(NfnI*Neq      * sizeof *nFluxNum_fI);       // free
		dnFluxNumdWIn_fI  = malloc(NfnI*Neq*Nvar * sizeof *dnFluxNumdWIn_fI);  // free
		dnFluxNumdWOut_fI = malloc(NfnI*Neq*Nvar * sizeof *dnFluxNumdWOut_fI); // free

		NFluxData->WL_fIL          = WIn_fI;
		NFluxData->WR_fIL          = WOut_fIIn;
		NFluxData->nFluxNum_fI     = nFluxNum_fI;
		NFluxData->dnFluxNumdWL_fI = dnFluxNumdWIn_fI;
		NFluxData->dnFluxNumdWR_fI = dnFluxNumdWOut_fI;

		compute_numerical_flux(FDATAL,'I');
		add_Jacobian_scaling_FACE(FDATAL,'I');

		// Compute FACE RHS and LHS terms

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


n_fI = FACE->n_fI;
nOrdInOut = OPSIn[IndFType]->nOrdInOut;
nOrdOutIn = OPSIn[IndFType]->nOrdOutIn;

detJF_fI = FACE->detJF_fI;

		NvnSIn = OPSIn[IndFType]->NvnS;
		NvnSOut = OPSOut[0]->NvnS;


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
			finalize_FACE_Inviscid_Weak(FDATAL,'L','E');
			finalize_FACE_Inviscid_Weak(FDATAL,'L','I');

			if (!FACE->Boundary) {
				// LHS (Currently not using sum factorization)

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

				// RHS
				swap_FACE_orientation(FDATAR,'I');
				finalize_FACE_Inviscid_Weak(FDATAR,'R','E');

				// LHS

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

	free(FDATAL);
	free(FDATAR);
	free(NFluxData);
}
