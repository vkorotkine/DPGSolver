// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

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

#include "element_functions.h"
#include "sum_factorization.h"
#include "matrix_functions.h"
#include "boundary_conditions.h"
#include "jacobian_boundary_conditions.h"
#include "fluxes_inviscid.h"
#include "jacobian_fluxes_inviscid.h"
#include "array_swap.h"

/*
#include "exact_solutions.h" // ToBeDeleted
#include "variable_functions.h" // ToBeDeleted
#include "array_print.h" // ToBeDeleted
#include <math.h> // ToBeDeleted
#include "array_norm.h" // ToBeDeleted
*/

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

struct S_OPERATORS {
	unsigned int NvnS, NfnI, NfnS, NvnS_SF, NfnI_SF, NvnI_SF, *nOrdInOut, *nOrdOutIn;
	double       **ChiS_fI, **ChiS_vI, **I_Weak_FF, **I_Weak_VV, **ChiS_fS, **GfS_fI;

	struct S_OpCSR **ChiS_fI_sp, **I_Weak_FF_sp;
};

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const struct S_FACE *FACE,
                     const unsigned int IndClass);
static void compute_FACE_EFE    (void);

void implicit_FACE_info(void)
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
			compute_FACE_EFE();
			break;
		default:
			compute_FACE_EFE();
//			compute_FACEVec_EFE();
			break;
		}
		break;
	default: // ADAPT_P, ADAPT_H, ADAPT_HP
printf("Error: Should not be entering default in implicit_FACE_info.\n"), exit(1);
		switch (Vectorized) {
		case 0:
//			compute_FACE();
			break;
		default:
//			compute_FACE();
//			compute_FACEVec();
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
printf("Error: Should not be entering default in implicit_FACE_info.\n"), exit(1);
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
printf("Error: Should not be entering default in implicit_FACE_info.\n"), exit(1);
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

static void compute_FACE_EFE(void)
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
	unsigned int i, j, n, iInd, dim, P, eq, var, iMax,
	             InddnFdWIn, InddnFdWOut, InddWOutdWIn, IndnF, IndI, Indeqvar, IndLHS,
	             VfIn, VfOut, fIn, fOut, EclassIn, EclassOut, IndFType, Boundary, BC, BC_trail, SpOpIn, SpOpOut,
	             RowInd, RowSub, ReOrder, *RowTracker,
	             NfnI, NvnSIn, NvnSOut, *nOrdOutIn, *nOrdInOut,
	             NIn[3], NOut[3], Diag[3], NOut0, NOut1, NIn0, NIn1;
	double       *WIn_fI, *WOut_fI, *WOut_fIIn,
	             *nFluxNum_fI, *dnFluxNumdWIn_fI, *dnFluxNumdWOut_fI, *dWOutdWIn,
	             *RHSIn, *RHSOut, *LHSInIn, *LHSOutIn, *LHSInOut, *LHSOutOut, *n_fI, *detJF_fI, *I_FF, *IdnFdW,
	             *ChiS_fI, *ChiS_fIOutIn, *ChiS_fIInOut,
	             *OP[3], **OPF0, **OPF1;

	struct S_OPERATORS *OPSIn[2], *OPSOut[2];
	struct S_VOLUME    *VIn, *VOut;
	struct S_FACE     *FACE;

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
/*
double *UEx;
UEx = malloc(NVAR3D*NfnI * sizeof *UEx); // free
compute_exact_solution(NfnI,FACE->XYZ_fI,UEx,0);
convert_variables(UEx,WOut_fIIn,3,d,NfnI,1,'p','c');
free(UEx);
*/
			} else {
				printf("Error: Unsupported BC in implicit_FACE_info.\n"), exit(1);
			}
		}

		// Compute numerical flux and its Jacobian
		nFluxNum_fI       = malloc(NfnI*Neq      * sizeof *nFluxNum_fI);       // free
		dnFluxNumdWIn_fI  = malloc(NfnI*Neq*Nvar * sizeof *dnFluxNumdWIn_fI);  // free
		dnFluxNumdWOut_fI = malloc(NfnI*Neq*Nvar * sizeof *dnFluxNumdWOut_fI); // free

		detJF_fI = FACE->detJF_fI;

		switch (InviscidFluxType) {
		case FLUX_LF:
			flux_LF(NfnI,1,WIn_fI,WOut_fIIn,nFluxNum_fI,n_fI,d,Neq);
			jacobian_flux_LF(NfnI,1,WIn_fI,WOut_fIIn,dnFluxNumdWIn_fI,n_fI,d,Neq,'L');
			jacobian_flux_LF(NfnI,1,WIn_fI,WOut_fIIn,dnFluxNumdWOut_fI,n_fI,d,Neq,'R');
			break;
		case FLUX_ROE:
			flux_Roe(NfnI,1,WIn_fI,WOut_fIIn,nFluxNum_fI,n_fI,d,Neq);
/*
if (BC % BC_STEP_SC == BC_SLIPWALL) {
// Modify computed nFluxNum_fI by using the solution to the Riemann problem for a reflection.
// See Toro(2009): eq. 6.23 to 6.25.

double *U_fI, *P_fI, *u_fI, *v_fI, *rho_fI;

U_fI = malloc(NfnI*Nvar * sizeof *U_fI); // free

convert_variables(WIn_fI,U_fI,d,d,NfnI,1,'c','p');
P_fI = &U_fI[NfnI*(Nvar-1)];
rho_fI = U_fI;
u_fI = &U_fI[NfnI*1];
v_fI = &U_fI[NfnI*2];

double n1, n2, uL, vL, pL, p, AL, BL, rhoL, aL, Vn;
for (n = 0; n < NfnI; n++) {
	n1 = n_fI[n*d+0];
	n2 = n_fI[n*d+1];

	rhoL = rho_fI[n];
	uL = u_fI[n];
	vL = v_fI[n];
	pL = P_fI[n];

	Vn = uL*n1+vL*n2;
	if (Vn < 0.0) { // Rarefaction
		aL = sqrt(GAMMA*pL/rhoL);
		p = pL*pow((1.0+0.5*GM1*Vn/aL),2.0*GAMMA/GM1);
	} else { // Shock
		AL = 2.0/((GAMMA+1)*rhoL);
		BL = GM1/(GAMMA+1)*pL;
		p = pL+0.5*Vn/AL*(Vn+sqrt(Vn*Vn+4*AL*(pL+BL)));
	}

	nFluxNum_fI[NfnI*0+n] = 0.0;
//	nFluxNum_fI[NfnI*1+n] = pL*n1;
//	nFluxNum_fI[NfnI*2+n] = pL*n2;
	nFluxNum_fI[NfnI*1+n] = p*n1;
	nFluxNum_fI[NfnI*2+n] = p*n2;
	nFluxNum_fI[NfnI*3+n] = 0.0;

	pL = p;
	p = pL;
}
free(U_fI);
//EXIT_MSG;
}
*/
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
				jacobian_boundary_Riemann(NfnI,1,FACE->XYZ_fI,WIn_fI,NULL,dWOutdWIn,n_fI,d,Neq);
			else if (BC % BC_STEP_SC == BC_SLIPWALL)
				jacobian_boundary_SlipWall(NfnI,1,WIn_fI,dWOutdWIn,n_fI,d,Neq);
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

				// LHS

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
