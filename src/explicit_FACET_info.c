// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Evaluate the FACET contributions to the RHS term.
 *
 *	Comments:
 *		When adaptivity is used, exact flux evaluation on the FACETs cannot be used in 3D whenever TRIs are present
 *		because the cubature nodes do not form a basis for the mortar element. This implies that adaptivity necessarily
 *		introduces aliasing errors unless the WSH nodes are used for the FACET integration. Of course, aliasing will
 *		always be introduced for the lower order FACET because of the projection, but this seemed not to be significant
 *		based on previous results using PF = P+1 in the VOLUME and interpolating to integration nodes giving optimal
 *		convergence orders.
 *
 *		For WEDGE ELEMENTs, the nOrd arrays for QUAD FACETs are stored with the TRI OPs, while those for TRI FACETs are
 *		stored with the LINE OPs. While this is not logical, it precludes the need for an additional OP structure.
 *
 *		When adding in MPI functionality, will not have access to VIn/VOut => Store what is needed for each FACET in an
 *		MPI communication routine before running this function. (ToBeDeleted)
 *
 *		Vectorization is more involved for FACET terms as there are many more possible combinations than for VOLUME
 *		terms. Given that the VOLUME vectorization seems not to have a significant impact on performance, the FACET
 *		vectorization may not be pursued. (ToBeDeleted)
 *
 *	Notation:
 *
 *	References:
 */

struct S_OPERATORS {
	unsigned int NvnS, NfnI, NfnS, NvnS_SF, NfnI_SF, NvnI_SF, *nOrdInOut, *nOrdOutIn;
	double       **ChiS_fI, **ChiS_vI, **I_Weak_FF, **I_Weak_VV, **GvShat_fS, **GfS_fI;

	struct S_OpCSR **ChiS_fI_sp, **I_Weak_FF_sp;
};

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const struct S_FACET *FACET,
                     const unsigned int IndClass);
static void compute_FACET_RHS_EFE    (void);
static void compute_FACET_RHS        (void);
//static void compute_FACETVec_RHS_EFE (void);
//static void compute_FACETVec_RHS     (void);

void explicit_FACET_info(void)
{
	// Initialize DB Parameters
	unsigned int Vectorized = DB.Vectorized,
	             Adapt      = DB.Adapt;

	switch (Adapt) {
	case ADAPT_0:
		switch (Vectorized) {
		case 0:
			compute_FACET_RHS_EFE();
			break;
		default:
			compute_FACET_RHS_EFE();
//			compute_FACETVec_RHS_EFE();
			break;
		}
		break;
	case ADAPT_P:
	case ADAPT_H:
	case ADAPT_HP:
	default:
		switch (Vectorized) {
		case 0:
			compute_FACET_RHS();
			/* Difference: Interpolate to FACET basis, evaluate numerical flux, interpolate the cubature nodes on each
			 *             side. Requires n_Sf (Solution facet) for curved elements. Make sure that Galerkin projections
			 *             are used.
			 * Note:       On TP FACETs, this routine should give identical results to EFE for conforming meshes.
			 */
			break;
		default:
			compute_FACET_RHS();
//			compute_FACETVec_RHS();
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
	OPS->GvShat_fS = ELEMENT->GvShat_fS[PV][PF];
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

		switch (Adapt) {
		default: // ADAPT_P, ADAPT_H, ADAPT_HP
			OPS->nOrdInOut = ELEMENT_FACET->nOrd_fS[PF][IndOrdInOut];
			OPS->nOrdOutIn = ELEMENT_FACET->nOrd_fS[PF][IndOrdOutIn];
			break;
		case ADAPT_0:
			OPS->nOrdInOut = ELEMENT_FACET->nOrd_fIs[PF][IndOrdInOut];
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

		switch (Adapt) {
		default: // ADAPT_P, ADAPT_H, ADAPT_HP
			OPS->nOrdInOut = ELEMENT_FACET->nOrd_fS[PF][IndOrdInOut];
			OPS->nOrdOutIn = ELEMENT_FACET->nOrd_fS[PF][IndOrdOutIn];
			break;
		case ADAPT_0:
			OPS->nOrdInOut = ELEMENT_FACET->nOrd_fIc[PF][IndOrdInOut];
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
	             InviscidFluxType = DB.InviscidFluxType,
	             Collocated       = DB.Collocated,
	             *VFPartUnity     = DB.VFPartUnity,
	             ***SF_BE         = DB.SF_BE;

	// Standard datatypes
	unsigned int i, j, iInd, dim, P,
	             VfIn, VfOut, fIn, fOut, EclassIn, EclassOut, IndFType, Boundary, BC, BC_trail,
	             RowInd, RowSub, ReOrder, *RowTracker,
	             NfnI, NvnSIn, NvnSOut, *nOrdOutIn, *nOrdInOut,
	             NIn[3], NOut[3], Diag[3], NOut0, NOut1, NIn0, NIn1;
	double       *WIn_fI, *WOut_fI, *WOut_fIIn, *nFluxNum_fI, *RHSIn, *RHSOut, *n_fI, *detJF_fI, *OP[3], **OPF0, **OPF1;

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

	for (FACET = DB.FACET; FACET != NULL; FACET = FACET->next) {
		P = FACET->P;

		// Obtain operators
		VIn  = FACET->VIn;
		VfIn = FACET->VfIn;
		fIn  = VfIn/NfrefMax;

		EclassIn = VIn->Eclass;
		IndFType = get_IndFType(EclassIn,fIn);
		init_ops(OPSIn[0],VIn,FACET,0);
		if (VIn->type == WEDGE || VIn->type == PYR)
			init_ops(OPSIn[1],VIn,FACET,1);

		VOut  = FACET->VOut;
		VfOut = FACET->VfOut;
		fOut  = VfOut/NfrefMax;

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

// Note: Needs modification for h-adaptation (ToBeDeleted)
			if (Collocated) {
				for (dim = 0; dim < d; dim++)
					Diag[dim] = 2;
				Diag[fIn/2] = 0;
			} else {
				for (dim = 0; dim < d; dim++)
					Diag[dim] = 0;
				Diag[1] = 2;
			}

			sf_apply_d(VIn->What,WIn_fI,NIn,NOut,Nvar,OP,Diag,d);
		} else if (EclassIn == C_WEDGE && SF_BE[P][1][1]) {
			if (fIn < 3) { OPF0  = OPSIn[0]->ChiS_fI, OPF1  = OPSIn[1]->ChiS_vI;
			               NOut0 = OPSIn[0]->NfnI_SF, NOut1 = OPSIn[1]->NvnI_SF;
			} else {       OPF0  = OPSIn[0]->ChiS_vI, OPF1  = OPSIn[1]->ChiS_fI;
			               NOut0 = OPSIn[0]->NvnI_SF, NOut1 = OPSIn[1]->NfnI_SF; }
			get_sf_parametersF(OPSIn[0]->NvnS_SF,NOut0,OPF0,OPSIn[1]->NvnS_SF,NOut1,OPF1,NIn,NOut,OP,d,VfIn,C_WEDGE);

// Note: Needs modification for h-adaptation (ToBeDeleted)
			if (Collocated) {
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
		} else if ((Collocated && (EclassIn == C_TP || EclassIn == C_WEDGE)) || (VFPartUnity[EclassIn])) {
// Note: May need modification for h-adaptation (ToBeDeleted)
//       The operator is not necessarily sparse in this case.
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
				if (Collocated) {
					for (dim = 0; dim < d; dim++)
						Diag[dim] = 2;
					Diag[fOut/2] = 0;
				} else {
					for (dim = 0; dim < d; dim++)
						Diag[dim] = 0;
					Diag[1] = 2;
				}

				sf_apply_d(VOut->What,WOut_fI,NIn,NOut,Nvar,OP,Diag,d);
			} else if (EclassOut == C_WEDGE && SF_BE[P][1][1]) {
				if (fOut < 3) { OPF0  = OPSOut[0]->ChiS_fI, OPF1  = OPSOut[1]->ChiS_vI;
				                NOut0 = OPSOut[0]->NfnI_SF, NOut1 = OPSOut[1]->NvnI_SF;
				} else {        OPF0  = OPSOut[0]->ChiS_vI, OPF1  = OPSOut[1]->ChiS_fI;
				                NOut0 = OPSOut[0]->NvnI_SF, NOut1 = OPSOut[1]->NfnI_SF; }
				get_sf_parametersF(OPSOut[0]->NvnS_SF,NOut0,OPF0,OPSOut[1]->NvnS_SF,NOut1,OPF1,NIn,NOut,OP,d,VfOut,C_WEDGE);

// Note: Needs modification for h-adaptation (ToBeDeleted)
				if (Collocated) {
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
			} else if ((Collocated && (EclassOut == C_TP || EclassOut == C_WEDGE)) || (VFPartUnity[EclassOut])) {
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
				boundary_Riemann(NfnI,1,FACET->XYZ_fI,WIn_fI,NULL,WOut_fIIn,n_fI);
			} else if (BC % BC_STEP_SC == BC_SLIPWALL) {
				boundary_SlipWall(NfnI,1,WIn_fI,WOut_fIIn,n_fI);
			} else {
				printf("Error: Unsupported BC in explicit_FACET_info.\n"), exit(1);
			}
		}

/*
printf("%d\n",FACET->indexg);
array_print_d(NfnI,Nvar,WIn_fI,'C');
array_print_d(NfnI,Nvar,WOut_fIIn,'C');
//if (FACET->indexg == 2)
//	exit(1);
//exit(1);
*/

		// Compute numerical flux
		nFluxNum_fI = malloc(NfnI*Neq * sizeof *nFluxNum_fI); // free
		detJF_fI = FACET->detJF_fI;

		switch (InviscidFluxType) {
		case FLUX_LF:
			flux_LF(NfnI,1,WIn_fI,WOut_fIIn,nFluxNum_fI,n_fI,d,Neq);
			break;
		case FLUX_ROE:
			flux_ROE(NfnI,1,WIn_fI,WOut_fIIn,nFluxNum_fI,n_fI,d,Neq);
			break;
		default:
			printf("Error: Unsupported InviscidFluxType used in explicit_FACET_info.\n"), exit(1);
			break;
		}

/*
//if (FACET->indexg == 4) {
printf("%d %d %d %d %d\n",FACET->indexg,IndFType,VIn->indexg,VOut->indexg,VfIn);
array_print_d(NfnI,Neq,WIn_fI,'C');
array_print_d(NfnI,Neq,WOut_fIIn,'C');
array_print_d(NfnI,d,n_fI,'R');
array_print_d(NfnI,1,detJF_fI,'C');
array_print_d(NfnI,Neq,nFluxNum_fI,'C');
//}
*/

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

		if (strstr(Form,"Weak") != NULL) {
			// Interior FACET
			if (EclassIn == C_TP && SF_BE[P][0][1]) {
				get_sf_parametersF(OPSIn[0]->NvnI_SF,OPSIn[0]->NvnS_SF,OPSIn[0]->I_Weak_VV,
								   OPSIn[0]->NfnI_SF,OPSIn[0]->NvnS_SF,OPSIn[0]->I_Weak_FF,NIn,NOut,OP,d,VfIn,C_TP);

// Note: Needs modification for h-adaptation (ToBeDeleted)
				if (Collocated) {
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

// Note: Needs modification for h-adaptation (ToBeDeleted)
				if (Collocated) {
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
			} else if ((Collocated && (EclassIn == C_TP || EclassIn == C_WEDGE)) || (VFPartUnity[EclassIn])) {
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

// Note: Needs modification for h-adaptation (ToBeDeleted)
					if (Collocated) {
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

// Note: Needs modification for h-adaptation (ToBeDeleted)
					if (Collocated) {
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
				} else if ((Collocated && (EclassOut == C_TP || EclassOut == C_WEDGE)) || (VFPartUnity[EclassOut])) {
					mm_CTN_CSR_d(NvnSOut,Neq,NfnI,OPSOut[0]->I_Weak_FF_sp[VfOut],nFluxNum_fI,RHSOut);
				} else {
					mm_CTN_d(NvnSOut,Neq,NfnI,OPSOut[0]->I_Weak_FF[VfOut],nFluxNum_fI,RHSOut);
				}
			}
		} else if (strstr(Form,"Strong") != NULL) {
			printf("Exiting: Implement the strong form in compute_FACET_RHS_EFE.\n"), exit(1);
		}

/*
//if (FACET->indexg == 2) {
printf("%d %d %d %d %d %d %d\n",FACET->indexg,IndFType,VIn->indexg,VOut->indexg,VfIn,BC,Boundary);
array_print_d(NvnSIn,Neq,RHSIn,'C');
array_print_d(NvnSOut,Neq,RHSOut,'C');
//array_print_d(NfnI,d,n_fI,'R');
//array_print_d(NfnI,d,FACET->XYZ_fI,'C');
//array_print_d(NfnI,Neq,WIn_fI,'C');
//array_print_d(NfnI,Neq,WOut_fIIn,'C');
//array_print_d(NfnI,Neq,nFluxNum_fI,'R');
//exit(1);
//}
*/

		free(RowTracker);
		free(WIn_fI);
		free(WOut_fIIn);
		free(nFluxNum_fI);
	}
//exit(1);
//printf("\n\n\n\n\n");

	for (i = 0; i < 2; i++) {
		free(OPSIn[i]);
		free(OPSOut[i]);
	}
}

//static void compute_FACETVec_RHS_EFE(void) { }

static void compute_FACET_RHS(void)
{
	/*
	 *	Comments:
	 *		Possibly include Sum Factorization here after this is verified.
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
	             NIn[3], NOut[3], Diag[3], NOut0, NOut1, NIn0, NIn1;
	double       *WIn_fS, *WOut_fS, *WOut_fSIn, *nFluxNum_fS, *nFluxNum_fI, *RHSIn, *RHSOut, *n_fS, *detJF_fS,
	             *nFluxNum_fIIn, *nFluxNum_fIOut,
	             *OP[3], **OPF0, **OPF1;

	struct S_OPERATORS *OPSIn[2], *OPSOut[2];
	struct S_VOLUME    *VIn, *VOut;
	struct S_FACET     *FACET;

	// silence

	for (i = 0; i < 2; i++) {
		OPSIn[i]  = malloc(sizeof *OPSIn[i]);  // free
		OPSOut[i] = malloc(sizeof *OPSOut[i]); // free
	}

	for (FACET = DB.FACET; FACET != NULL; FACET = FACET->next) {
		P = FACET->P;

		// Obtain operators
		VIn  = FACET->VIn;
		VfIn = FACET->VfIn;
		fIn  = VfIn/NfrefMax;

		EclassIn = VIn->Eclass;
		IndFType = get_IndFType(EclassIn,fIn);
		init_ops(OPSIn[0],VIn,FACET,0);
		if (VIn->type == WEDGE || VIn->type == PYR)
			init_ops(OPSIn[1],VIn,FACET,1);

		VOut  = FACET->VOut;
		VfOut = FACET->VfOut;
		fOut  = VfOut/NfrefMax;

		EclassOut = VOut->Eclass;
		init_ops(OPSOut[0],VOut,FACET,0);
		if (VOut->type == WEDGE || VOut->type == PYR)
			init_ops(OPSOut[1],VOut,FACET,1);

		BC = FACET->BC;
		Boundary = !((VIn->indexg != VOut->indexg) || (VIn->indexg == VOut->indexg && fIn != fOut));
		// The second condition is for periodic elements which are connected to themselves

		// Compute WIn_fS
		NfnS = OPSIn[IndFType]->NfnS;
		NvnSIn = OPSIn[IndFType]->NvnS;

		WIn_fS = malloc(NfnS*Nvar * sizeof *WIn_fS); // free
		// If VFPartUnity == 1, this operator would be sparse.
		mm_CTN_d(NfnS,Nvar,NvnSIn,OPSIn[0]->GvShat_fS[VfIn],VIn->What,WIn_fS);

		// Compute WOut_fI (Taking BCs into account if applicable)
		n_fS     = FACET->n_fS;

		nOrdInOut = OPSIn[IndFType]->nOrdInOut;
		nOrdOutIn = OPSIn[IndFType]->nOrdOutIn;

		NvnSOut = OPSOut[0]->NvnS;
		WOut_fSIn = malloc(NfnS*Nvar * sizeof *WOut_fSIn); // free
		if (BC == 0 || (BC % BC_STEP_SC > 50)) { // Internal/Periodic FACET
			WOut_fS   = malloc(NfnS*Nvar * sizeof *WOut_fS); // free
			mm_CTN_d(NfnS,Nvar,NvnSOut,OPSOut[0]->GvShat_fS[VfOut],VOut->What,WOut_fS);

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
		} else { // Boundary FACET
			if (BC % BC_STEP_SC == BC_RIEMANN) {
				boundary_Riemann(NfnS,1,FACET->XYZ_fS,WIn_fS,NULL,WOut_fSIn,n_fS);
			} else if (BC % BC_STEP_SC == BC_SLIPWALL) {
				boundary_SlipWall(NfnS,1,WIn_fS,WOut_fSIn,n_fS);
			} else {
				printf("Error: Unsupported BC in explicit_FACET_info.\n"), exit(1);
			}
		}
/*
array_print_d(NfnS,Nvar,WIn_fS,'C');
array_print_d(NfnS,Nvar,WOut_fSIn,'C');
*/
		// Compute numerical flux
		nFluxNum_fS = malloc(NfnS*Neq * sizeof *nFluxNum_fS); // free
		detJF_fS = FACET->detJF_fS;

		switch (InviscidFluxType) {
		case FLUX_LF:
			flux_LF(NfnS,1,WIn_fS,WOut_fSIn,nFluxNum_fS,n_fS,d,Neq);
			break;
		case FLUX_ROE:
			flux_ROE(NfnS,1,WIn_fS,WOut_fSIn,nFluxNum_fS,n_fS,d,Neq);
			break;
		default:
			printf("Error: Unsupported InviscidFluxType used in explicit_FACET_info.\n"), exit(1);
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

		// Galerkin projection to FACET cubature nodes
		NfnI   = OPSIn[IndFType]->NfnI;

		RowTracker = malloc(NfnS * sizeof *RowTracker); // free
// Don't forget to reorder nF_I for Outer FACET

		switch (Adapt) {
		default: // ADAPT_HP
			break;
		case ADAPT_P:
		case ADAPT_H:
			nFluxNum_fIIn = malloc(NfnI*Neq * sizeof *nFluxNum_fIIn); // free
			mm_CTN_d(NfnI,Neq,NfnS,OPSIn[0]->GfS_fI[VfIn],nFluxNum_fS,nFluxNum_fIIn);

// Potentially perform the rearrangement after projection to fI nodes.
// In this case, make sure that nOrdInOut is for the cubature nodes.
			// Re-arrange nFluxNum to match node ordering from opposite VOLUME
			for (i = 0; i < NfnS; i++)
				RowTracker[i] = i;

			for (RowInd = 0; RowInd < NfnS; RowInd++) {
				ReOrder = nOrdInOut[RowInd];
				for (RowSub = ReOrder; RowTracker[RowSub] != ReOrder; RowSub = RowTracker[RowSub])
					;

				if (RowInd != RowSub) {
					array_swap_d(&nFluxNum_fS[RowInd],&nFluxNum_fS[RowSub],Neq,NfnS);
					array_swap_ui(&RowTracker[RowInd],&RowTracker[RowSub],1,1);
				}
			}
			nFluxNum_fIOut = malloc(NfnI*Neq * sizeof *nFluxNum_fIOut); // free
			mm_CTN_d(NfnI,Neq,NfnS,OPSOut[0]->GfS_fI[VfOut],nFluxNum_fS,nFluxNum_fIOut);
			// -ve sign (opposite normal vector) is added below

			break;
		// Note: Need not sum the contributions from all of the refined FACETs as they will be taken care of based on
		//       other FACETs. Need to ensure that RHSOut is only initialized once. ADAPT_P and ADAPT_H should then be
		//       extremely similar.
		}
		free(RowTracker);
		free(nFluxNum_fS);
/*
printf("%d\n",NfnI);
array_print_d(NfnI,Neq,nFluxNum_fIIn,'C');
array_print_d(NfnI,Neq,nFluxNum_fIOut,'C');
*/
		// Compute FACET RHS terms
		RHSIn  = calloc(NvnSIn*Neq  , sizeof *RHSIn);  // keep (requires external free)
		RHSOut = calloc(NvnSOut*Neq , sizeof *RHSOut); // keep (requires external free)
		FACET->RHSIn  = RHSIn;
		FACET->RHSOut = RHSOut;

		if (strstr(Form,"Weak") != NULL) {
			// Interior FACET
			if (EclassIn == C_TP && SF_BE[P][0][1]) {
				get_sf_parametersF(OPSIn[0]->NvnI_SF,OPSIn[0]->NvnS_SF,OPSIn[0]->I_Weak_VV,
								   OPSIn[0]->NfnI_SF,OPSIn[0]->NvnS_SF,OPSIn[0]->I_Weak_FF,NIn,NOut,OP,d,VfIn,C_TP);

				if (Collocated && Adapt == ADAPT_0) {
					for (dim = 0; dim < d; dim++)
						Diag[dim] = 2;
					Diag[fIn/2] = 0;
				} else {
					for (dim = 0; dim < d; dim++)
						Diag[dim] = 0;
				}

				sf_apply_d(nFluxNum_fIIn,RHSIn,NIn,NOut,Neq,OP,Diag,d);
			} else if (EclassIn == C_WEDGE && SF_BE[P][1][1]) {
				if (fIn < 3) { OPF0 = OPSIn[0]->I_Weak_FF, OPF1 = OPSIn[1]->I_Weak_VV;
							   NIn0 = OPSIn[0]->NfnI_SF,   NIn1 = OPSIn[1]->NvnI_SF;
				} else {       OPF0 = OPSIn[0]->I_Weak_VV, OPF1 = OPSIn[1]->I_Weak_FF;
							   NIn0 = OPSIn[0]->NvnI_SF,   NIn1 = OPSIn[1]->NfnI_SF; }
				get_sf_parametersF(NIn0,OPSIn[0]->NvnS_SF,OPF0,NIn1,OPSIn[1]->NvnS_SF,OPF1,NIn,NOut,OP,d,VfIn,C_WEDGE);

// Note: Needs modification for h-adaptation (ToBeDeleted)
				if (Collocated) {
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

				sf_apply_d(nFluxNum_fIIn,RHSIn,NIn,NOut,Neq,OP,Diag,d);
			} else if ((Collocated && (EclassIn == C_TP || EclassIn == C_WEDGE)) || (VFPartUnity[EclassIn])) {
				mm_CTN_CSR_d(NvnSIn,Neq,NfnI,OPSIn[0]->I_Weak_FF_sp[VfIn],nFluxNum_fIIn,RHSIn);
			} else  {
				mm_CTN_d(NvnSIn,Neq,NfnI,OPSIn[0]->I_Weak_FF[VfIn],nFluxNum_fIIn,RHSIn);
			}

			// Exterior FACET
			if (!Boundary) {
				// Use -ve normal for opposite FACET
				for (i = 0; i < Neq; i++) {
					iInd = i*NfnI;
					for (j = 0; j < NfnI; j++)
						nFluxNum_fIOut[iInd+j] *= -1.0;
				}

/*
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
*/

				if (EclassOut == C_TP && SF_BE[P][0][1]) {
					get_sf_parametersF(OPSOut[0]->NvnI_SF,OPSOut[0]->NvnS_SF,OPSOut[0]->I_Weak_VV,
									   OPSOut[0]->NfnI_SF,OPSOut[0]->NvnS_SF,OPSOut[0]->I_Weak_FF,NIn,NOut,OP,d,VfOut,C_TP);

// Note: Needs modification for h-adaptation (ToBeDeleted)
					if (Collocated) {
						for (dim = 0; dim < d; dim++)
							Diag[dim] = 2;
						Diag[fOut/2] = 0;
					} else {
						for (dim = 0; dim < d; dim++)
							Diag[dim] = 0;
					}

					sf_apply_d(nFluxNum_fIOut,RHSOut,NIn,NOut,Neq,OP,Diag,d);
				} else if (EclassOut == C_WEDGE && SF_BE[P][1][1]) {
					if (fOut < 3) { OPF0 = OPSOut[0]->I_Weak_FF, OPF1 = OPSOut[1]->I_Weak_VV;
					                NIn0 = OPSOut[0]->NfnI_SF,   NIn1 = OPSOut[1]->NvnI_SF;
					} else {        OPF0 = OPSOut[0]->I_Weak_VV, OPF1 = OPSOut[1]->I_Weak_FF;
					                NIn0 = OPSOut[0]->NvnI_SF,   NIn1 = OPSOut[1]->NfnI_SF; }
					get_sf_parametersF(NIn0,OPSOut[0]->NvnS_SF,OPF0,NIn1,OPSOut[1]->NvnS_SF,OPF1,NIn,NOut,OP,d,VfOut,C_WEDGE);

// Note: Needs modification for h-adaptation (ToBeDeleted)
					if (Collocated) {
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

					sf_apply_d(nFluxNum_fIOut,RHSOut,NIn,NOut,Neq,OP,Diag,d);
				} else if ((Collocated && (EclassOut == C_TP || EclassOut == C_WEDGE)) || (VFPartUnity[EclassOut])) {
					mm_CTN_CSR_d(NvnSOut,Neq,NfnI,OPSOut[0]->I_Weak_FF_sp[VfOut],nFluxNum_fIOut,RHSOut);
				} else {
					mm_CTN_d(NvnSOut,Neq,NfnI,OPSOut[0]->I_Weak_FF[VfOut],nFluxNum_fIOut,RHSOut);
				}
			}
		} else if (strstr(Form,"Strong") != NULL) {
			printf("Exiting: Implement the strong form in compute_FACET_RHS_EFE.\n"), exit(1);
		}
		free(nFluxNum_fIIn);
		free(nFluxNum_fIOut);
//exit(1);
	}

	for (i = 0; i < 2; i++) {
		free(OPSIn[i]);
		free(OPSOut[i]);
	}
//exit(1);
}

//static void compute_FACETVec_RHS(void) { }
