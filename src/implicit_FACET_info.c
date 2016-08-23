// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "implicit_FACET_info.h"

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
#include "boundary_conditions.h"
#include "jacobian_boundary_conditions.h"
#include "fluxes_inviscid.h"
#include "jacobian_fluxes_inviscid.h"
#include "array_swap.h"

#include "array_print.h"

/*
 *	Purpose:
 *		Evaluate the FACET contributions to the LHS term.
 *
 *	Comments:
 *		Check for relevant comments in implicit_VOLUME_info.
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
static void compute_FACET_EFE    (void);

void implicit_FACET_info(void)
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
			compute_FACET_EFE();
			break;
		default:
			compute_FACET_EFE();
//			compute_FACETVec_EFE();
			break;
		}
		break;
	default: // ADAPT_P, ADAPT_H, ADAPT_HP
printf("Error: Should not be entering default in implicit_FACET_info.\n"), exit(1);
		switch (Vectorized) {
		case 0:
//			compute_FACET();
			break;
		default:
//			compute_FACET();
//			compute_FACETVec();
			break;
		}
		break;
	}
}

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const struct S_FACET *FACET,
                     const unsigned int IndClass)
{
	// Initialize DB Parameters
	unsigned int Adapt    = DB.Adapt,
	             ***SF_BE = DB.SF_BE;

	// Standard datatypes
	unsigned int PV, PF, Vtype, Eclass, FtypeInt, IndOrdInOut, IndOrdOutIn;

	struct S_ELEMENT *ELEMENT, *ELEMENT_OPS, *ELEMENT_FACET;

	// silence
	ELEMENT_OPS = NULL;

	PV       = VOLUME->P;
	PF       = FACET->P;
	Vtype    = VOLUME->type;
	Eclass   = VOLUME->Eclass;

	FtypeInt    = FACET->typeInt;
	IndOrdInOut = FACET->IndOrdInOut;
	IndOrdOutIn = FACET->IndOrdOutIn;

	ELEMENT       = get_ELEMENT_type(Vtype);
	ELEMENT_FACET = get_ELEMENT_FACET(Vtype,IndClass);
	if ((Eclass == C_TP && SF_BE[PF][0][1]) || (Eclass == C_WEDGE && SF_BE[PF][1][1]))
		ELEMENT_OPS = ELEMENT->ELEMENTclass[IndClass];
	else
		ELEMENT_OPS = ELEMENT;

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
printf("Error: Should not be entering default in implicit_FACET_info.\n"), exit(1);
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
printf("Error: Should not be entering default in implicit_FACET_info.\n"), exit(1);
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

static void compute_FACET_EFE(void)
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
	unsigned int i, j, n, iInd, dim, P, eq, var,
	             InddnFdWIn, InddnFdWOut, InddWOutdWIn,
	             VfIn, VfOut, fIn, fOut, EclassIn, EclassOut, IndFType, Boundary, BC, BC_trail, SpOpIn, SpOpOut,
	             RowInd, RowSub, ReOrder, *RowTracker,
	             NfnI, NvnSIn, NvnSOut, *nOrdOutIn, *nOrdInOut,
	             NIn[3], NOut[3], Diag[3], NOut0, NOut1, NIn0, NIn1;
	double       *WIn_fI, *WOut_fI, *WOut_fIIn,
	             *nFluxNum_fI, *dnFluxNumdWIn_fI, *dnFluxNumdWOut_fI, *dWOutdWIn,
	             *RHSIn, *RHSOut, *n_fI, *detJF_fI,
	             *OP[3], **OPF0, **OPF1;

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
		P = FACET->P;

		// Obtain operators
		VIn    = FACET->VIn;
		VfIn   = FACET->VfIn;
		fIn    = VfIn/NfrefMax;
		SpOpIn = Collocated && (VfIn % NFREFMAX == 0 && VIn->P == P);

		EclassIn = VIn->Eclass;
		IndFType = get_IndFType(EclassIn,fIn);
		init_ops(OPSIn[0],VIn,FACET,0);
		if (VIn->type == WEDGE || VIn->type == PYR)
			init_ops(OPSIn[1],VIn,FACET,1);

		VOut    = FACET->VOut;
		VfOut   = FACET->VfOut;
		fOut    = VfOut/NfrefMax;
		SpOpOut = Collocated && (VfOut % NFREFMAX == 0 && VOut->P == P);

		EclassOut = VOut->Eclass;
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
		n_fI     = FACET->n_fI;

		nOrdInOut = OPSIn[IndFType]->nOrdInOut;
		nOrdOutIn = OPSIn[IndFType]->nOrdOutIn;

		NvnSOut = OPSOut[0]->NvnS;
		WOut_fI = malloc(NfnI*Nvar * sizeof *WOut_fI); // free
		WOut_fIIn = malloc(NfnI*Nvar * sizeof *WOut_fIIn); // free
		if (BC == 0 || (BC % BC_STEP_SC > 50)) { // Internal/Periodic FACET
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
		} else { // Boundary FACET
			if (BC % BC_STEP_SC == BC_RIEMANN) {
				boundary_Riemann(NfnI,1,FACET->XYZ_fI,WIn_fI,NULL,WOut_fIIn,n_fI,d);
			} else if (BC % BC_STEP_SC == BC_SLIPWALL) {
				boundary_SlipWall(NfnI,1,WIn_fI,WOut_fIIn,n_fI,d);
			} else {
				printf("Error: Unsupported BC in implicit_FACET_info.\n"), exit(1);
			}
		}

		// Compute numerical flux and its jacobian
		nFluxNum_fI       = malloc(NfnI*Neq      * sizeof *nFluxNum_fI);       // free
		dnFluxNumdWIn_fI  = malloc(NfnI*Neq*Nvar * sizeof *dnFluxNumdWIn_fI);  // free
		dnFluxNumdWOut_fI = malloc(NfnI*Neq*Nvar * sizeof *dnFluxNumdWOut_fI); // free

		detJF_fI = FACET->detJF_fI;

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
			printf("Error: Unsupported InviscidFluxType.\n"), EXIT_MSG;
			break;
		}

		// Include the BC information in dnFluxNumWIn_fI if on a boundary
		if (Boundary) {
			dWOutdWIn = malloc(NfnI*Nvar*Nvar * sizeof *dWOutdWIn); // free
			if (BC % BC_STEP_SC == BC_RIEMANN)
				jacobian_boundary_Riemann(NfnI,1,FACET->XYZ_fI,WIn_fI,NULL,dWOutdWIn,n_fI,d,Neq);
			else if (BC % BC_STEP_SC == BC_SLIPWALL)
				jacobian_boundary_SlipWall(NfnI,1,WIn_fI,dWOutdWIn,n_fI,d,Neq);
			else
				printf("Error: Unsupported BC.\n"), EXIT_MSG;

			for (eq = 0; eq < Neq; eq++) {
			for (var = 0; var < Nvar; var++) {
				InddnFdWIn = (eq*Neq+var)*NfnI;

				for (i = 0; i < Nvar; i++) {
					InddnFdWOut  = (eq*Neq+i)*NfnI;
					InddWOutdWIn = (var*Nvar+i)*NfnI;
					for (n = 0; n < NfnI; n++)
						dnFluxNumdWIn_fI[InddnFdWIn+n] += dnFluxNumdWOut_fI[InddnFdWOut+n]*dWOutdWIn[InddWOutdWIn+n];
				}
			}}

			free(dWOutdWIn);
		}
if (FACET->indexg == 0) {
printf("impF: %d %d\n",FACET->indexg,Boundary);
array_print_d(NfnI,Neq,nFluxNum_fI,'C');
array_print_d(NfnI*Neq,Neq,dnFluxNumdWIn_fI,'C');
}
continue here.

		// Multiply n dot FNum by the area element
		for (i = 0; i < Neq; i++) {
			iInd = i*NfnI;
			for (j = 0; j < NfnI; j++)
				nFluxNum_fI[iInd+j] *= detJF_fI[j];
		}

		// Compute FACET RHS terms
		RHSIn  = calloc(NvnSIn*Neq  , sizeof *RHSIn);  // keep (requires external free)
		RHSOut = calloc(NvnSOut*Neq , sizeof *RHSOut); // keep (requires external free)
		FACET->RHSIn  = RHSIn;
		FACET->RHSOut = RHSOut;

		RowTracker = malloc(NfnI * sizeof *RowTracker); // free

		if (strstr(Form,"Weak")) {
			// Interior FACET
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
			printf("Exiting: Implement the strong form in compute_FACET_EFE.\n"), exit(1);
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
