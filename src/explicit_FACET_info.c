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
 *		Need to add in support for sum factorization when working (ToBeDeleted).
 *		When adding in MPI functionality, will not have access to VIn/VOut => Store what is needed for each FACET in an
 *		MPI communication routine before running this function. (ToBeDeleted)
 *
 *	Notation:
 *
 *	References:
 */

struct S_OPERATORS {
	unsigned int NvnS, NfnI;
	double       **ChiS_fI;
};

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const struct S_FACET *FACET,
                     const unsigned int IndClass);
static void compute_FACET_RHS_EFE    (void);
static void compute_FACETVec_RHS_EFE (void);
static void compute_FACET_RHS        (void);
static void compute_FACETVec_RHS     (void);

void explicit_FACET_info(void)
{
	// Initialize DB Parameters
	unsigned int Vectorized = DB.Vectorized,
	             Adaptive   = DB.Adaptive;

	switch (Adaptive) {
	case 0:
		switch (Vectorized) {
		case 0:
			compute_FACET_RHS_EFE();
			break;
		default:
			compute_FACETVec_RHS_EFE();
			break;
		}
		break;
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
			compute_FACETVec_RHS();
			break;
		}
		break;
	}
}

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const struct S_FACET *FACET,
                     const unsigned int IndClass)
{
	/*
	 *	Comments:
	 *		VCurved may not be needed here in the if condition => check when done. (ToBeDeleted)
	 */

	// Standard datatypes
	unsigned int PV, PF, Vtype, Vcurved, FtypeInt;

	struct S_ELEMENT *ELEMENT, *ELEMENT_OPS;

	// silence
	ELEMENT_OPS = NULL;

	PV       = VOLUME->P;
	PF       = FACET->P;
	Vtype    = VOLUME->type;
	Vcurved  = VOLUME->curved;
	FtypeInt = FACET->typeInt;

	ELEMENT = get_ELEMENT_type(Vtype);
	if (Vtype == LINE || Vtype == QUAD || Vtype == HEX || Vtype == WEDGE)
		ELEMENT_OPS = ELEMENT->ELEMENTclass[IndClass];
	else if (Vtype == TRI || Vtype == TET || Vtype == PYR)
		ELEMENT_OPS = ELEMENT;

	OPS->NvnS = ELEMENT_OPS->NvnS[PV];
	if (!Vcurved) {
		// Straight VOLUME

		if (FtypeInt == 's') {
			// Straight FACET Integration
			OPS->NfnI = ELEMENT_OPS->NfnIs[PF][IndClass];

			OPS->ChiS_fI = ELEMENT_OPS->ChiS_fIs[PV][PF];
		} else {
			// Curved FACET Integration
			OPS->NfnI = ELEMENT_OPS->NfnIc[PF][IndClass];

			OPS->ChiS_fI = ELEMENT_OPS->ChiS_fIc[PV][PF];
		}
	} else {
		// Curved VOLUME

		if (FtypeInt == 's') {
			// Straight FACET Integration
			OPS->NfnI = ELEMENT_OPS->NfnIs[PF][IndClass];

			OPS->ChiS_fI = ELEMENT_OPS->ChiS_fIs[PV][PF];
		} else {
			// Curved FACET Integration
			OPS->NfnI = ELEMENT_OPS->NfnIc[PF][IndClass];

			OPS->ChiS_fI = ELEMENT_OPS->ChiS_fIc[PV][PF];
		}
	}
}

static void compute_FACET_RHS_EFE(void)
{
	// Initialize DB Parameters
	unsigned int NfrefMax = DB.NfrefMax,
	             Nvar     = DB.Nvar,
	             Neq      = DB.Neq;

	// Standard datatypes
	unsigned int i,
	             VfIn, VfOut, fIn, fOut, Eclass, IndFType, Boundary, BC,
	             NfnI, NvnS;
	double       *WIn_fI, *WOut_fI, *nFluxNum_fI, *RHSIn, *RHSOut;

	struct S_OPERATORS *OPSIn, *OPSOut;
	struct S_VOLUME    *VIn, *VOut;
	struct S_FACET     *FACET;

	OPSIn  = malloc(sizeof *OPSIn);  // free
	OPSOut = malloc(sizeof *OPSOut); // free

	for (FACET = DB.FACET; FACET != NULL; FACET = FACET->next) {
		// Obtain operators
		VIn  = FACET->VIn;
		VfIn = FACET->VfIn;
		fIn  = VfIn/NfrefMax;

		Eclass = get_Eclass(VIn->type);
		IndFType = get_IndFType(Eclass,fIn);
		init_ops(OPSIn,VIn,FACET,IndFType);

		VOut  = FACET->VOut;
		VfOut = FACET->VfOut;
		fOut  = VfOut/NfrefMax;

		Eclass = get_Eclass(VOut->type);
		IndFType = get_IndFType(Eclass,fOut);
		init_ops(OPSOut,VOut,FACET,IndFType);

		BC = FACET->BC;
		Boundary = !((VIn->indexg != VOut->indexg) || (VIn->indexg == VOut->indexg && fIn != fOut));
		// The second condition is for periodic elements which are connected to themselves

		// Compute WIn_fI
		if (0 && VIn->Eclass == C_TP) {
			; // update this with sum factorization
		} else if (1 || VIn->Eclass == C_SI || VIn->Eclass == C_PYR) {
			NfnI = OPSIn->NfnI;
			NvnS = OPSIn->NvnS;

			WIn_fI = malloc(NfnI*Nvar * sizeof *WIn_fI); // free
			mm_CTN_d(NfnI,Nvar,NvnS,OPSIn->ChiS_fI[VfIn],VIn->What,WIn_fI);
		} else if (VIn->Eclass == C_WEDGE) {
			; // update this with sum factorization
		}

		// Compute WOut_fI (Taking BCs into account if applicable)
		if (BC == 0 || (BC % BC_STEP_SC > 50)) { // Internal/Periodic FACET
			if (0 && VOut->Eclass == C_TP) {
				; // update this with sum factorization
			} else if (1 || VOut->Eclass == C_SI || VOut->Eclass == C_PYR) {
				WOut_fI = malloc(NfnI*Nvar * sizeof *WOut_fI); // free

printf("%d %d\n",VfIn,VfOut);
array_print_d(NfnI,NvnS,OPSOut->ChiS_fI[VfOut],'R');

				mm_CTN_d(NfnI,Nvar,NvnS,OPSOut->ChiS_fI[VfOut],VOut->What,WOut_fI);
			} else if (VOut->Eclass == C_WEDGE) {
				; // update this with sum factorization
			}

			// Reorder WOut_fI to correspond to WIn_fI
array_print_d(NfnI,Nvar,WIn_fI,'C');
array_print_d(NfnI,Nvar,WOut_fI,'C');

			printf("Working on the vOutInIO.\n");
			exit(1);
		} else { // Boundary FACET
			; // Add in support for boundary conditions
		}



		// Compute numerical flux
		nFluxNum_fI = malloc(NfnI*Neq * sizeof *nFluxNum_fI); // free












/*
		NvnSIn  = OPSIn->NvnS;
		NvnSOut = OPSOut->NvnS;

		RHSIn  = calloc(NvnSIn*Neq  , sizeof *RHSIn);  // keep (requires external free)
		RHSOut = calloc(NvnSOut*Neq , sizeof *RHSOut); // keep (requires external free)
		FACET->RHSIn  = RHSIn;
		FACET->RHSOut = RHSOut;
*/

		free(WIn_fI);
		free(WOut_fI);
		free(nFluxNum_fI);

	}

	free(OPSIn);
	free(OPSOut);
}

static void compute_FACETVec_RHS_EFE(void)
{

}

static void compute_FACET_RHS(void)
{

}

static void compute_FACETVec_RHS(void)
{

}
