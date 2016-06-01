// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

// Can likely be removed later (ToBeDeleted)
#include "math.h"

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
	unsigned int NvnS, NfnI, *nOrdInOut, *nOrdOutIn;
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
	unsigned int PV, PF, Vtype, Vcurved, FtypeInt, IndOrdInOut, IndOrdOutIn;

	struct S_ELEMENT *ELEMENT, *ELEMENT_OPS, *ELEMENT_FACET;

	// silence
	ELEMENT_OPS = NULL;

	PV       = VOLUME->P;
	PF       = FACET->P;
	Vtype    = VOLUME->type;
	Vcurved  = VOLUME->curved;

	FtypeInt    = FACET->typeInt;
	IndOrdInOut = FACET->IndOrdInOut;
	IndOrdOutIn = FACET->IndOrdOutIn;

	ELEMENT       = get_ELEMENT_type(Vtype);
	ELEMENT_FACET = get_ELEMENT_FACET(Vtype,IndClass);
	if (1 || Vtype == TRI || Vtype == TET || Vtype == PYR)
		ELEMENT_OPS = ELEMENT;
	else if (Vtype == LINE || Vtype == QUAD || Vtype == HEX || Vtype == WEDGE)
		ELEMENT_OPS = ELEMENT->ELEMENTclass[IndClass];

	OPS->NvnS = ELEMENT_OPS->NvnS[PV];
	if (!Vcurved) {
		// Straight VOLUME

		if (FtypeInt == 's') {
			// Straight FACET Integration
			OPS->NfnI = ELEMENT_OPS->NfnIs[PF][IndClass];

			OPS->ChiS_fI = ELEMENT_OPS->ChiS_fIs[PV][PF];

			OPS->nOrdInOut = ELEMENT_FACET->nOrd_fIs[PF][IndOrdInOut];
			OPS->nOrdOutIn = ELEMENT_FACET->nOrd_fIs[PF][IndOrdOutIn];
		} else {
			// Curved FACET Integration
			OPS->NfnI = ELEMENT_OPS->NfnIc[PF][IndClass];

			OPS->ChiS_fI = ELEMENT_OPS->ChiS_fIc[PV][PF];

			OPS->nOrdInOut = ELEMENT_FACET->nOrd_fIc[PF][IndOrdInOut];
			OPS->nOrdOutIn = ELEMENT_FACET->nOrd_fIc[PF][IndOrdOutIn];
		}
	} else {
		// Curved VOLUME

		if (FtypeInt == 's') {
			// Straight FACET Integration
			OPS->NfnI = ELEMENT_OPS->NfnIs[PF][IndClass];

			OPS->ChiS_fI = ELEMENT_OPS->ChiS_fIs[PV][PF];

			OPS->nOrdInOut = ELEMENT_FACET->nOrd_fIs[PF][IndOrdInOut];
			OPS->nOrdOutIn = ELEMENT_FACET->nOrd_fIs[PF][IndOrdOutIn];
		} else {
			// Curved FACET Integration
			OPS->NfnI = ELEMENT_OPS->NfnIc[PF][IndClass];

			OPS->ChiS_fI = ELEMENT_OPS->ChiS_fIc[PV][PF];

			OPS->nOrdInOut = ELEMENT_FACET->nOrd_fIc[PF][IndOrdInOut];
			OPS->nOrdOutIn = ELEMENT_FACET->nOrd_fIc[PF][IndOrdOutIn];
		}
	}
}

static void compute_FACET_RHS_EFE(void)
{
	// Initialize DB Parameters
	unsigned int d                = DB.d,
	             NfrefMax         = DB.NfrefMax,
	             Nvar             = DB.Nvar,
	             Neq              = DB.Neq,
	             InviscidFluxType = DB.InviscidFluxType;

	// Standard datatypes
	unsigned int i, j, iInd,
	             VfIn, VfOut, fIn, fOut, Eclass, IndFType, Boundary, BC, BC_trail,
	             NfnI, NvnS;
	double       *WIn_fI, *WOut_fI, *WOut_fIIn, *nFluxNum_fI, *RHSIn, *RHSOut, *n_fI;

	struct S_OPERATORS *OPSIn, *OPSOut;
	struct S_VOLUME    *VIn, *VOut;
	struct S_FACET     *FACET;
	struct S_ELEMENT   *ELEMENT;

	// silence
	WOut_fI   = NULL;
	WOut_fIIn = NULL;

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

				mm_CTN_d(NfnI,Nvar,NvnS,OPSOut->ChiS_fI[VfOut],VOut->What,WOut_fI);
			} else if (VOut->Eclass == C_WEDGE) {
				; // update this with sum factorization
			}

			// Reorder WOut_fI to correspond to WIn_fI
			ELEMENT = get_ELEMENT_type(VIn->type);

			WOut_fIIn = malloc(NfnI*Nvar * sizeof *WOut_fIIn); // free
			for (i = 0; i < Nvar; i++) {
				iInd = i*NfnI;
				for (j = 0; j < NfnI; j++) {
					BC_trail = BC % BC_STEP_SC;
					if (BC_trail > 50 || BC_trail == 0)
						WOut_fIIn[iInd+j] = WOut_fI[iInd+OPSIn->nOrdInOut[j]];
					else
						WOut_fIIn[iInd+j] = WOut_fI[iInd+j];
				}
			}

/*
printf("%d %d %d %d\n",FACET->indexg,BC,VfIn,VfOut);
printf("%d %d\n",VfIn,VfOut);
printf("%d %d\n",FACET->IndOrdInOut,FACET->IndOrdOutIn);
array_print_ui(1,NfnI,OPSIn->nOrdInOut,'R');

for (i = 0; i < NfnI; i++) {
	for (j = 0; j < Nvar; j++) {
		printf(" % .3e",WIn_fI[i+NfnI*j]-WOut_fIIn[i+NfnI*j]);
	}
	printf("\n");
}
printf("\n");

if (0 && FACET->indexg == 0) {
	array_print_d(NfnI,Nvar,WIn_fI,'C');
	array_print_d(NfnI,Nvar,WOut_fI,'C');
	array_print_d(NfnI,Nvar,WOut_fIIn,'C');
	exit(1);
}
*/
		} else { // Boundary FACET
			; // Add in support for boundary conditions
		}

		// Compute numerical flux
		nFluxNum_fI = malloc(NfnI*Neq * sizeof *nFluxNum_fI); // free
		n_fI = FACET->n;

for (i = 0; i < NfnI*Nvar; i++)
	WOut_fIIn[i] += 0.01*i;
for (i = 0; i < NfnI*d; i++)
	n_fI[i] += 0.02*i;

double tmp_d;

for (i = 0; i < NfnI; i++) {
	tmp_d = 0;
	for (j = 0; j < d; j++) {
		tmp_d += pow(n_fI[i*d+j],2.0);
	}
	tmp_d = sqrt(tmp_d);
	for (j = 0; j < d; j++) {
		n_fI[i*d+j] /= tmp_d;
	}
}

double *WIn2, *WOut2, *n2;

WIn2 = malloc(NfnI*Neq * sizeof *WIn2);
WOut2 = malloc(NfnI*Neq * sizeof *WOut2);
n2 = malloc(NfnI*d * sizeof *n2);

for (i = 0; i < NfnI; i++) {
	for (j = 0; j < Neq; j++) {
		WIn2[i*Neq+j] = (double) (((unsigned int) (WIn_fI[i*Neq+j]*1e3))/1e3);
		WOut2[i*Neq+j] = (double) (((unsigned int) (WOut_fIIn[i*Neq+j]*1e3))/1e3);
	}
	for (j = 0; j < d; j++) {
		n2[i*d+j] = (double) (((unsigned int) (n_fI[i*d+j]*1e4))/1e4);
	}
}

array_print_d(NfnI,Nvar,WIn2,'C');
array_print_d(NfnI,Nvar,WOut2,'C');
array_print_d(NfnI,d,n2,'R');

flux_LF(NfnI,1,WIn2,WOut2,nFluxNum_fI,n2,d,Neq);
array_print_d(NfnI,Neq,nFluxNum_fI,'C');

flux_ROE(NfnI,1,WIn2,WOut2,nFluxNum_fI,n2,d,Neq);
array_print_d(NfnI,Neq,nFluxNum_fI,'C');
exit(1);

		switch (InviscidFluxType) {
		case FLUX_LF:  flux_LF(NfnI,1,WIn_fI,WOut_fIIn,nFluxNum_fI,n_fI,d,Neq); break;
		case FLUX_ROE:
			; // write the Roe scheme
			break;
		}

array_print_d(NfnI,Neq,nFluxNum_fI,'C');

exit(1);











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
		free(WOut_fIIn);
		free(nFluxNum_fI);

	}
	exit(1);

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
