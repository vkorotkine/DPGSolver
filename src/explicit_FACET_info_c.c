// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "explicit_FACET_info_c.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "S_FACET.h"
#include "S_OpCSR.h"

#include "element_functions.h"
#include "sum_factorization.h"
#include "matrix_functions.h"
#include "boundary_conditions_c.h"
#include "fluxes_inviscid_c.h"
#include "array_swap.h"

#include "array_print.h"

/*
 *	Purpose:
 *		Identical to explicit_FACET_info using complex variables (for complex step verification).
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

struct S_OPERATORS {
	unsigned int NvnS, NfnI, NfnS, NvnS_SF, NfnI_SF, NvnI_SF, *nOrdInOut, *nOrdOutIn;
	double       **ChiS_fI, **ChiS_vI, **I_Weak_FF, **I_Weak_VV, **ChiS_fS, **GfS_fI;

	struct S_OpCSR **ChiS_fI_sp, **I_Weak_FF_sp;
};

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const struct S_FACET *FACET,
                     const unsigned int IndClass);
static void compute_FACET_RHS_EFE    (void);

void explicit_FACET_info_c(void)
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
			compute_FACET_RHS_EFE();
			break;
		default:
			compute_FACET_RHS_EFE();
			break;
		}
		break;
	default: // ADAPT_P, ADAPT_H, ADAPT_HP
printf("Error: Should not be entering default in explicit_FACET_info.\n"), exit(1);
		switch (Vectorized) {
		case 0:
//			compute_FACET_RHS();
			break;
		default:
//			compute_FACET_RHS();
			break;
		}
		break;
	}
}

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const struct S_FACET *FACET,
                     const unsigned int IndClass)
{
	// Initialize DB Parameters
	unsigned int Adapt    = DB.Adapt;

	// Standard datatypes
	unsigned int PV, PF, Vtype, FtypeInt, IndOrdInOut, IndOrdOutIn;

	struct S_ELEMENT *ELEMENT, *ELEMENT_OPS, *ELEMENT_FACET;

	// silence

	PV       = VOLUME->P;
	PF       = FACET->P;
	Vtype    = VOLUME->type;

	FtypeInt    = FACET->typeInt;
	IndOrdInOut = FACET->IndOrdInOut;
	IndOrdOutIn = FACET->IndOrdOutIn;

	ELEMENT       = get_ELEMENT_type(Vtype);
	ELEMENT_FACET = get_ELEMENT_FACET(Vtype,IndClass);
	ELEMENT_OPS   = ELEMENT;

	OPS->NvnS    = ELEMENT->NvnS[PV];
	OPS->NvnS_SF = ELEMENT_OPS->NvnS[PV];

	OPS->NfnS      = ELEMENT_FACET->NvnS[PF];
	OPS->ChiS_fS   = ELEMENT->ChiS_fS[PV][PF];
	if (FtypeInt == 's') {
		// Straight FACET Integration
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

		OPS->nOrdInOut = ELEMENT_FACET->nOrd_fIs[PF][IndOrdInOut];
		switch (Adapt) {
		default: // ADAPT_P, ADAPT_H, ADAPT_HP
printf("Error: Should not be entering default in explicit_FACET_info.\n"), exit(1);
			OPS->nOrdOutIn = ELEMENT_FACET->nOrd_fS[PF][IndOrdOutIn];
			break;
case ADAPT_P: // ToBeModified (Also change setup_normals and output_to_paraview)
case ADAPT_H:
case ADAPT_HP:
		case ADAPT_0:
			OPS->nOrdOutIn = ELEMENT_FACET->nOrd_fIs[PF][IndOrdOutIn];
			break;
		}
	} else {
		// Curved FACET Integration
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

		OPS->nOrdInOut = ELEMENT_FACET->nOrd_fIc[PF][IndOrdInOut];
		switch (Adapt) {
		default: // ADAPT_P, ADAPT_H, ADAPT_HP
printf("Error: Should not be entering default in explicit_FACET_info.\n"), exit(1);
			OPS->nOrdOutIn = ELEMENT_FACET->nOrd_fS[PF][IndOrdOutIn];
			break;
case ADAPT_P: // ToBeModified (Also change setup_normals and output_to_paraview)
case ADAPT_H:
case ADAPT_HP:
		case ADAPT_0:
			OPS->nOrdOutIn = ELEMENT_FACET->nOrd_fIc[PF][IndOrdOutIn];
			break;
		}
	}
}

static void compute_FACET_RHS_EFE(void)
{
	// Initialize DB Parameters
	char         *Form            = DB.Form;
	unsigned int d                = DB.d,
	             NfrefMax         = DB.NfrefMax,
	             Nvar             = DB.Nvar,
	             Neq              = DB.Neq,
	             InviscidFluxType = DB.InviscidFluxType;

	// Standard datatypes
	unsigned int   i, j, iInd,
	               VfIn, VfOut, fIn, fOut, EclassIn, IndFType, Boundary, BC, BC_trail,
	               RowInd, RowSub, ReOrder, *RowTracker,
	               NfnI, NvnSIn, NvnSOut, *nOrdOutIn, *nOrdInOut;
	double         *n_fI, *detJF_fI;
	double complex *WIn_fI, *WOut_fI, *WOut_fIIn, *nFluxNum_fI, *RHSIn, *RHSOut;

	struct S_OPERATORS *OPSIn[2], *OPSOut[2];
	struct S_VOLUME    *VIn, *VOut;
	struct S_FACET     *FACET;

	// silence
	WOut_fI   = NULL;
	WOut_fIIn = NULL;
	NvnSOut   = 0;

	for (i = 0; i < 2; i++) {
		OPSIn[i]  = malloc(sizeof *OPSIn[i]);  // free
		OPSOut[i] = malloc(sizeof *OPSOut[i]); // free
	}

	for (FACET = DB.FACET; FACET; FACET = FACET->next) {
		// Obtain operators
		VIn    = FACET->VIn;
		VfIn   = FACET->VfIn;
		fIn    = VfIn/NfrefMax;

		EclassIn = VIn->Eclass;
		IndFType = get_IndFType(EclassIn,fIn);
		init_ops(OPSIn[0],VIn,FACET,0);
		if (VIn->type == WEDGE || VIn->type == PYR)
			init_ops(OPSIn[1],VIn,FACET,1);

		VOut    = FACET->VOut;
		VfOut   = FACET->VfOut;
		fOut    = VfOut/NfrefMax;

		init_ops(OPSOut[0],VOut,FACET,0);
		if (VOut->type == WEDGE || VOut->type == PYR)
			init_ops(OPSOut[1],VOut,FACET,1);

		BC = FACET->BC;
		Boundary = !((VIn->indexg != VOut->indexg) || (VIn->indexg == VOut->indexg && fIn != fOut));
		// The second condition is for periodic elements which are connected to themselves

		// Compute WIn_fI
		NfnI   = OPSIn[IndFType]->NfnI;
		NvnSIn = OPSIn[IndFType]->NvnS;

		WIn_fI = malloc(NfnI*Nvar * sizeof *WIn_fI); // free
		mm_dcc(CBCM,CBT,CBNT,NfnI,Nvar,NvnSIn,1.0,OPSIn[0]->ChiS_fI[VfIn],VIn->What_c,WIn_fI);

		// Compute WOut_fI (Taking BCs into account if applicable)
		n_fI     = FACET->n_fI;

		nOrdInOut = OPSIn[IndFType]->nOrdInOut;
		nOrdOutIn = OPSIn[IndFType]->nOrdOutIn;

		NvnSOut = OPSOut[0]->NvnS;
		WOut_fI = malloc(NfnI*Nvar * sizeof *WOut_fI); // free
		WOut_fIIn = malloc(NfnI*Nvar * sizeof *WOut_fIIn); // free
		if (BC == 0 || (BC % BC_STEP_SC > 50)) { // Internal/Periodic FACET
			mm_dcc(CBCM,CBT,CBNT,NfnI,Nvar,NvnSOut,1.0,OPSOut[0]->ChiS_fI[VfOut],VOut->What_c,WOut_fI);

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
		} else { // Boundary FACET
			if (BC % BC_STEP_SC == BC_RIEMANN) {
				boundary_Riemann_c(NfnI,1,FACET->XYZ_fI,WIn_fI,NULL,WOut_fIIn,n_fI,d);
			} else if (BC % BC_STEP_SC == BC_SLIPWALL) {
				boundary_SlipWall_c(NfnI,1,WIn_fI,WOut_fIIn,n_fI,d);
			} else {
				printf("Error: Unsupported BC in explicit_FACET_info.\n"), exit(1);
			}
		}

		// Compute numerical flux
		nFluxNum_fI = malloc(NfnI*Neq * sizeof *nFluxNum_fI); // free
		detJF_fI = FACET->detJF_fI;

		switch (InviscidFluxType) {
		case FLUX_LF:
			flux_LF_c(NfnI,1,WIn_fI,WOut_fIIn,nFluxNum_fI,n_fI,d,Neq);
			break;
		case FLUX_ROE:
			flux_Roe_c(NfnI,1,WIn_fI,WOut_fIIn,nFluxNum_fI,n_fI,d,Neq);
			break;
		default:
			printf("Error: Unsupported InviscidFluxType used in explicit_FACET_info.\n"), exit(1);
			break;
		}
printf("expF: %d %d\n",FACET->indexg,Boundary);
array_print_cmplx(NfnI,Neq,nFluxNum_fI,'C');
EXIT_MSG;

		// Multiply n dot FNum by the area element
		for (i = 0; i < Neq; i++) {
			iInd = i*NfnI;
			for (j = 0; j < NfnI; j++)
				nFluxNum_fI[iInd+j] *= detJF_fI[j];
		}

		// Compute FACET RHS terms
		RHSIn  = calloc(NvnSIn*Neq  , sizeof *RHSIn);  // keep (requires external free)
		RHSOut = calloc(NvnSOut*Neq , sizeof *RHSOut); // keep (requires external free)
		if (FACET->RHSIn_c)
			free(FACET->RHSIn_c);
		FACET->RHSIn_c  = RHSIn;

		if (FACET->RHSOut_c)
			free(FACET->RHSOut_c);
		FACET->RHSOut_c = RHSOut;

		RowTracker = malloc(NfnI * sizeof *RowTracker); // free

		if (strstr(Form,"Weak")) {
			// Interior FACET
			mm_dcc(CBCM,CBT,CBNT,NvnSIn,Neq,NfnI,1.0,OPSIn[0]->I_Weak_FF[VfIn],nFluxNum_fI,RHSIn);

			// Exterior FACET
			if (!Boundary) {
				// Use -ve normal for opposite FACET
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
						array_swap_cmplx(&nFluxNum_fI[RowInd],&nFluxNum_fI[RowSub],Neq,NfnI);
						array_swap_ui(&RowTracker[RowInd],&RowTracker[RowSub],1,1);
					}
				}

				mm_dcc(CBCM,CBT,CBNT,NvnSOut,Neq,NfnI,1.0,OPSOut[0]->I_Weak_FF[VfOut],nFluxNum_fI,RHSOut);
			}
		} else if (strstr(Form,"Strong")) {
			printf("Exiting: Implement the strong form in compute_FACET_RHS_EFE.\n"), exit(1);
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
