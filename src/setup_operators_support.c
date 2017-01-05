// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "setup_operators_support.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"

#include "select_functions.h"
#include "element_functions.h"
#include "matrix_functions.h"
#include "adaptation.h"
#include "plotting_element_info.h"
#include "bases.h"
#include "cubature.h"
#include "solver_functions.h"

#include "array_print.h"

/*
 *	Purpose:
 *		Provide support functions for setup_operators.
 *
 *	Comments:
 *
 *	Notation:
 *		See setup_operators for additional notation.
 *
 *	setup_ELEMENT_plotting:
 *		connectivity  : connectivity between nodes to form P1 sub-elements (see test_imp_plotting for visualization)
 *		connect_types : VTK element type numbering for each sub-element
 *		                Note: This is not trivial as TETs and PYRs are refined into a combination of TETs and PYRs.
 *		connect_NE    : (N)umber of sub (E)lements in connectivity array
 *
 *	setup_ELEMENT_normals:
 *		Theta_[] : Angles for conversion between reference and face coordinates (Zwanenburg(2016): Table 14)
 *		           Options: eta, zeta
 *		nr       : (n)ormal vector components for the (r)eference element
 *
 *	get_rst_vV: get (rst) (v)olume coordinates of reference ELEMENT (V)ertices.
 *
 *	get_BCoords_dEm(*): get (B)ary(C)entric (Coord)inates of (d)imension (E)lement (m)inus (1 or 2)
 *
 *		Note: See cubature_PYR for explanation of the method to obtain the barycentric coordinates.
 *
 *	setup_ELEMENT_VeV:
 *	setup_ELEMENT_VeF:
 *		VeV : (Ve)rtex to (V)OLUME operator used to project VOLUME vertex nodes to VOLUME vertex nodes for all supported
 *		      h-refinements.
 *		VeF : (Ve)rtex to (F)ACE operator used to project VOLUME vertex nodes to FACE vertex nodes for all supported
 *		      h-refinements.
 *
 *	get_L2_scaling: self-explanatory.
 *
 *	setup_ELEMENT_FACE_ordering:
 *		nOrd_f(1)(2) : (n)ode (Ord)ering of (f)ace of type (1) which is (2).
 *			(1): (S)olution, (I)ntegration
 *			(2): (s)traight, (c)urved
 *
 *	References:
 *		See setup_operators.
 */

void setup_ELEMENT_plotting(const unsigned int EType)
{
	// Standard datatypes
	unsigned int P, PSMin, PSMax, NvnP, NE, u1 = 1,
	             *connectivity, *types;
	double       *rst_vP;

	struct S_ELEMENT *ELEMENT;

	ELEMENT = get_ELEMENT_type(EType);

	get_PS_range(&PSMin,&PSMax);
	for (P = PSMin; P <= PSMax; P++) {
		plotting_element_info(&rst_vP,&connectivity,&types,&NvnP,&NE,max(P,u1),EType); // free
		ELEMENT->connectivity[P]  = connectivity;
		ELEMENT->connect_types[P] = types;
		ELEMENT->connect_NE[P]    = NE;
		ELEMENT->NvnP[P] = NvnP;

		free(rst_vP);
	}
}

void setup_ELEMENT_normals(const unsigned int EType)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int f, Nf;
	double Theta_eta[6], Theta_zeta[6], *nr;

	struct S_ELEMENT *ELEMENT;

	ELEMENT = get_ELEMENT_type(EType);

	Nf = ELEMENT->Nf;
	nr = ELEMENT->nr;

	if (EType == LINE) {
		Theta_eta[0] = 0.0; Theta_zeta[0] = PI;
		Theta_eta[1] = 0.0; Theta_zeta[1] = 0.0;
	} else if (EType == QUAD) {
		Theta_eta[0] = 0.0; Theta_zeta[0] = PI;
		Theta_eta[1] = 0.0; Theta_zeta[1] = 0.0;
		Theta_eta[2] = 0.0; Theta_zeta[2] = 1.5*PI;
		Theta_eta[3] = 0.0; Theta_zeta[3] = 0.5*PI;
	} else if (EType == HEX) {
		Theta_eta[0] = 0.0;    Theta_zeta[0] = PI;
		Theta_eta[1] = 0.0;    Theta_zeta[1] = 0.0;
		Theta_eta[2] = 0.0;    Theta_zeta[2] = 1.5*PI;
		Theta_eta[3] = 0.0;    Theta_zeta[3] = 0.5*PI;
		Theta_eta[4] = 0.5*PI; Theta_zeta[4] = 0.0;
		Theta_eta[5] = 1.5*PI; Theta_zeta[5] = PI;
	} else if (EType == TRI) {
		Theta_eta[0] = 0.0; Theta_zeta[0] = 1.0/6.0*PI;
		Theta_eta[1] = 0.0; Theta_zeta[1] = 5.0/6.0*PI;
		Theta_eta[2] = 0.0; Theta_zeta[2] = 9.0/6.0*PI;
	} else if (EType == TET) {
		double Theta_e = atan(2.0*sqrt(2.0));

		Theta_eta[0] = Theta_e - 0.5*PI; Theta_zeta[0] = 1.0/6.0*PI;
		Theta_eta[1] = Theta_e - 0.5*PI; Theta_zeta[1] = 5.0/6.0*PI;
		Theta_eta[2] = Theta_e - 0.5*PI; Theta_zeta[2] = 9.0/6.0*PI;
		Theta_eta[3] = 0.5*PI;           Theta_zeta[3] = 0.0;
	} else if (EType == WEDGE) {
		Theta_eta[0] = 0.0;    Theta_zeta[0] = 1.0/6.0*PI;
		Theta_eta[1] = 0.0;    Theta_zeta[1] = 5.0/6.0*PI;
		Theta_eta[2] = 0.0;    Theta_zeta[2] = 9.0/6.0*PI;
		Theta_eta[3] = 0.5*PI; Theta_zeta[3] = 0.0;
		Theta_eta[4] = 1.5*PI; Theta_zeta[4] = 0.0;
	} else if (EType == PYR) {
		double Theta_e = atan(sqrt(2.0));

		Theta_eta[0] = Theta_e - 0.5*PI; Theta_zeta[0] = PI;
		Theta_eta[1] = Theta_e - 0.5*PI; Theta_zeta[1] = 0.0;
		Theta_eta[2] = Theta_e - 0.5*PI; Theta_zeta[2] = 1.5*PI;
		Theta_eta[3] = Theta_e - 0.5*PI; Theta_zeta[3] = 0.5*PI;
		Theta_eta[4] = 0.5*PI;           Theta_zeta[4] = 0.0;
	} else {
		printf("Error: Unsupported EType.\n"), EXIT_MSG;
	}

	for (f = 0; f < Nf; f++) {
		           nr[f*d+0] =  cos(Theta_eta[f])*cos(Theta_zeta[f]);
		if (d > 1) nr[f*d+1] =  cos(Theta_eta[f])*sin(Theta_zeta[f]);
		if (d > 2) nr[f*d+2] = -sin(Theta_eta[f]);
	}
}

double *get_rst_vV(const struct S_ELEMENT *ELEMENT)
{
	unsigned int d, Nve, EType;
	double *rst_vV;

	d     = ELEMENT->d;
	Nve   = ELEMENT->Nve;
	EType = ELEMENT->type;

	rst_vV = malloc(Nve*d * sizeof *rst_vV); // keep

	switch (EType) {
	case LINE:
		rst_vV[0*Nve+0] = -1.0;
		rst_vV[0*Nve+1] =  1.0;
		break;
	case TRI:
		rst_vV[0*Nve+0] = -1.0; rst_vV[1*Nve+0] = -1.0/sqrt(3.0);
		rst_vV[0*Nve+1] =  1.0; rst_vV[1*Nve+1] = -1.0/sqrt(3.0);
		rst_vV[0*Nve+2] =  0.0; rst_vV[1*Nve+2] =  2.0/sqrt(3.0);
		break;
	case QUAD:
		rst_vV[0*Nve+0] = -1.0; rst_vV[1*Nve+0] = -1.0;
		rst_vV[0*Nve+1] =  1.0; rst_vV[1*Nve+1] = -1.0;
		rst_vV[0*Nve+2] = -1.0; rst_vV[1*Nve+2] =  1.0;
		rst_vV[0*Nve+3] =  1.0; rst_vV[1*Nve+3] =  1.0;
		break;
	case TET:
		rst_vV[0*Nve+0] = -1.0; rst_vV[1*Nve+0] = -1.0/sqrt(3.0); rst_vV[2*Nve+0] = -1.0/sqrt(6.0);
		rst_vV[0*Nve+1] =  1.0; rst_vV[1*Nve+1] = -1.0/sqrt(3.0); rst_vV[2*Nve+1] = -1.0/sqrt(6.0);
		rst_vV[0*Nve+2] =  0.0; rst_vV[1*Nve+2] =  2.0/sqrt(3.0); rst_vV[2*Nve+2] = -1.0/sqrt(6.0);
		rst_vV[0*Nve+3] =  0.0; rst_vV[1*Nve+3] =  0.0;           rst_vV[2*Nve+3] =  3.0/sqrt(6.0);
		break;
	case HEX:
		rst_vV[0*Nve+0] = -1.0; rst_vV[1*Nve+0] = -1.0; rst_vV[2*Nve+0] = -1.0;
		rst_vV[0*Nve+1] =  1.0; rst_vV[1*Nve+1] = -1.0; rst_vV[2*Nve+1] = -1.0;
		rst_vV[0*Nve+2] = -1.0; rst_vV[1*Nve+2] =  1.0; rst_vV[2*Nve+2] = -1.0;
		rst_vV[0*Nve+3] =  1.0; rst_vV[1*Nve+3] =  1.0; rst_vV[2*Nve+3] = -1.0;
		rst_vV[0*Nve+4] = -1.0; rst_vV[1*Nve+4] = -1.0; rst_vV[2*Nve+4] =  1.0;
		rst_vV[0*Nve+5] =  1.0; rst_vV[1*Nve+5] = -1.0; rst_vV[2*Nve+5] =  1.0;
		rst_vV[0*Nve+6] = -1.0; rst_vV[1*Nve+6] =  1.0; rst_vV[2*Nve+6] =  1.0;
		rst_vV[0*Nve+7] =  1.0; rst_vV[1*Nve+7] =  1.0; rst_vV[2*Nve+7] =  1.0;
		break;
	case WEDGE:
		rst_vV[0*Nve+0] = -1.0; rst_vV[1*Nve+0] = -1.0/sqrt(3.0); rst_vV[2*Nve+0] = -1.0;
		rst_vV[0*Nve+1] =  1.0; rst_vV[1*Nve+1] = -1.0/sqrt(3.0); rst_vV[2*Nve+1] = -1.0;
		rst_vV[0*Nve+2] =  0.0; rst_vV[1*Nve+2] =  2.0/sqrt(3.0); rst_vV[2*Nve+2] = -1.0;
		rst_vV[0*Nve+3] = -1.0; rst_vV[1*Nve+3] = -1.0/sqrt(3.0); rst_vV[2*Nve+3] =  1.0;
		rst_vV[0*Nve+4] =  1.0; rst_vV[1*Nve+4] = -1.0/sqrt(3.0); rst_vV[2*Nve+4] =  1.0;
		rst_vV[0*Nve+5] =  0.0; rst_vV[1*Nve+5] =  2.0/sqrt(3.0); rst_vV[2*Nve+5] =  1.0;
		break;
	case PYR:
		rst_vV[0*Nve+0] = -1.0; rst_vV[1*Nve+0] = -1.0; rst_vV[2*Nve+0] = -1.0/5.0*sqrt(2.0);
		rst_vV[0*Nve+1] =  1.0; rst_vV[1*Nve+1] = -1.0; rst_vV[2*Nve+1] = -1.0/5.0*sqrt(2.0);
		rst_vV[0*Nve+2] = -1.0; rst_vV[1*Nve+2] =  1.0; rst_vV[2*Nve+2] = -1.0/5.0*sqrt(2.0);
		rst_vV[0*Nve+3] =  1.0; rst_vV[1*Nve+3] =  1.0; rst_vV[2*Nve+3] = -1.0/5.0*sqrt(2.0);
		rst_vV[0*Nve+4] =  0.0; rst_vV[1*Nve+4] =  0.0; rst_vV[2*Nve+4] =  4.0/5.0*sqrt(2.0);
		break;
	default:
		printf("Error: Unsupported EType.\n"), EXIT_MSG;
		break;
	}
	return rst_vV;
}

struct S_BCOORDS *get_BCoords_dEm1(const struct S_ELEMENT *ELEMENT, const unsigned int IndFType)
{
	// Initialize DB Parameters
	unsigned int NP   = DB.NP,
	             PMax = DB.PMax;

	// Standard datatypes
	unsigned int P, EType, Nve,
	             *NfnG2, *NfnGc, *NfnS, *NfnIs, *NfnIc;
	double       **w_fIs, **w_fIc,
	             **BCoords_G2, **BCoords_Gc, **BCoords_S, **BCoords_Is, **BCoords_Ic, *one;

	struct S_BCOORDS *BCoords_dEm1;
	struct S_ELEMENT *ELEMENT_F;

	BCoords_dEm1 = malloc(sizeof *BCoords_dEm1);      // keep
	NfnG2        = calloc(NP , sizeof *(NfnG2));      // keep
	NfnGc        = malloc(NP * sizeof *(NfnGc));      // keep
	NfnS         = malloc(NP * sizeof *(NfnS));       // keep
	NfnIs        = malloc(NP * sizeof *(NfnIs));      // keep
	NfnIc        = malloc(NP * sizeof *(NfnIc));      // keep
	w_fIs        = malloc(NP * sizeof *(w_fIs));      // keep
	w_fIc        = malloc(NP * sizeof *(w_fIc));      // keep
	BCoords_G2   = malloc(NP * sizeof *(BCoords_G2)); // keep
	BCoords_Gc   = malloc(NP * sizeof *(BCoords_Gc)); // keep
	BCoords_S    = malloc(NP * sizeof *(BCoords_S));  // keep
	BCoords_Is   = malloc(NP * sizeof *(BCoords_Is)); // keep
	BCoords_Ic   = malloc(NP * sizeof *(BCoords_Ic)); // keep

	one = malloc(1 * sizeof *one); // free
	one[0] = 1.0;

	EType = ELEMENT->type;
	if (EType == LINE) {
// fix this to be consistent with treatment of other elements if possible (ToBeDeleted)
		ELEMENT_F = get_ELEMENT_type(POINT);

		Nve = 1;
		for (P = 0; P <= PMax; P++) {
			NfnG2[2] = 1;
			NfnGc[P] = 1;
			NfnS[P]  = 1;
			NfnIs[P] = 1;
			NfnIc[P] = 1;

			w_fIs[P] = mm_Alloc_d(CBCM,CBT,CBT,1,1,1,1.0,one,one); // keep
			w_fIc[P] = mm_Alloc_d(CBCM,CBT,CBT,1,1,1,1.0,one,one); // keep

			BCoords_G2[P] = malloc(1 * sizeof **(BCoords_G2));
			BCoords_Gc[P] = malloc(1 * sizeof **(BCoords_Gc));
			BCoords_S[P]  = malloc(1 * sizeof **(BCoords_S));
			BCoords_Is[P] = malloc(1 * sizeof **(BCoords_Is));
			BCoords_Ic[P] = malloc(1 * sizeof **(BCoords_Ic));

			BCoords_G2[P][0] = 1.0;
			BCoords_Gc[P][0] = 1.0;
			BCoords_S[P][0]  = 1.0;
			BCoords_Is[P][0] = 1.0;
			BCoords_Ic[P][0] = 1.0;
		}
	} else {
		// Use get_ELEMENT_F_type here (ToBeDeleted)
		if (EType == TRI || EType == QUAD) {
			ELEMENT_F = get_ELEMENT_type(LINE);
		} else if (EType == TET || (EType == WEDGE && IndFType == 1) || (EType == PYR && IndFType == 0)) {
			ELEMENT_F = get_ELEMENT_type(TRI);
		} else if (EType == HEX || (EType == WEDGE && IndFType == 0) || (EType == PYR && IndFType == 1)) {
			ELEMENT_F = get_ELEMENT_type(QUAD);
		} else {
			printf("Error: Unsupported EType/IndFType combination.\n"), EXIT_MSG;
		}

		// Initialize DB Parameters
		unsigned int *PGc           = DB.PGc,
		             **PIfs         = DB.PIfs,
		             **PIfc         = DB.PIfc;
		char         **NodeTypeG    = DB.NodeTypeG,
		             ***NodeTypeS   = DB.NodeTypeS,
		             ***NodeTypeIfs = DB.NodeTypeIfs,
		             ***NodeTypeIfc = DB.NodeTypeIfc;

		// Standard datatypes
		unsigned int dE, Nbf,
		             NfnGs, EType, EclassF,
		             dummy_ui, *dummyPtr_ui;
		double       *rst_vGs, *rst_fG2, *rst_fGc, *rst_fS, *rst_fIs, *rst_fIc,
		             *dummyPtr_d[2],
		             *IGs,
		             *ChiRefGs_vGs,
		             *ChiRefInvGs_vGs,
		             *ChiRefGs_fG2, *ChiRefGs_fGc, *ChiRefGs_fS, *ChiRefGs_fIs, *ChiRefGs_fIc;

		// Function pointers
		cubature_tdef   cubature;
		basis_tdef      basis;
		grad_basis_tdef grad_basis;

		EType = ELEMENT_F->type;
		select_functions(&basis,&grad_basis,&cubature,EType);

		// Only use EclassF here as setting up face nodes for TRIs in 3D results in invalid LINE NodeTypes.
		EclassF = get_Eclass(EType);


		dE = ELEMENT_F->d;

		Nve = ELEMENT_F->Nve;

		// It is important to use the nodes corresponding to the VeF ordering
		rst_vGs = get_rst_vV(ELEMENT_F); // free
		cubature(&dummyPtr_d[0],&dummyPtr_d[1],&dummyPtr_ui,&NfnGs,&dummy_ui,0,1,dE,NodeTypeG[EclassF]); // free
		free(dummyPtr_ui);
		free(dummyPtr_d[0]);

		IGs             = identity_d(NfnGs);                       // free
		ChiRefGs_vGs    = basis(1,rst_vGs,NfnGs,&Nbf,dE);          // free
		ChiRefInvGs_vGs = inverse_d(NfnGs,NfnGs,ChiRefGs_vGs,IGs); // free

		free(rst_vGs);

		free(IGs);
		free(ChiRefGs_vGs);

		for (P = 0; P <= PMax; P++) {
			cubature(&rst_fGc,&dummyPtr_d[0],&dummyPtr_ui,&NfnGc[P],&dummy_ui,0,PGc[P],          dE,NodeTypeG[EclassF]);      free(dummyPtr_ui); // free
			cubature(&rst_fS, &dummyPtr_d[0],&dummyPtr_ui,&NfnS[P], &dummy_ui,0,P,               dE,NodeTypeS[P][EclassF]);   free(dummyPtr_ui); // free
			cubature(&rst_fIs,&w_fIs[P],     &dummyPtr_ui,&NfnIs[P],&dummy_ui,1,PIfs[P][EclassF],dE,NodeTypeIfs[P][EclassF]); free(dummyPtr_ui); // free
			cubature(&rst_fIc,&w_fIc[P],     &dummyPtr_ui,&NfnIc[P],&dummy_ui,1,PIfc[P][EclassF],dE,NodeTypeIfc[P][EclassF]); free(dummyPtr_ui); // free

			ChiRefGs_fGc = basis(1,rst_fGc,NfnGc[P],&Nbf,dE); // free
			ChiRefGs_fS  = basis(1,rst_fS, NfnS[P], &Nbf,dE); // free
			ChiRefGs_fIs = basis(1,rst_fIs,NfnIs[P],&Nbf,dE); // free
			ChiRefGs_fIc = basis(1,rst_fIc,NfnIc[P],&Nbf,dE); // free

			if (P == 2) {
				cubature(&rst_fG2,&dummyPtr_d[0],&dummyPtr_ui,&NfnG2[P],&dummy_ui,0,P,dE,NodeTypeG[EclassF]); free(dummyPtr_ui); // free
				ChiRefGs_fG2 = basis(1,rst_fG2,NfnG2[P],&Nbf,dE); // free
				BCoords_G2[P] = mm_Alloc_d(CBCM,CBT,CBT,NfnG2[P],NfnGs,NfnGs,1.0,ChiRefGs_fG2,ChiRefInvGs_vGs); // keep

				free(rst_fG2);
				free(ChiRefGs_fG2);
			}

			BCoords_Gc[P] = mm_Alloc_d(CBCM,CBT,CBT,NfnGc[P],NfnGs,NfnGs,1.0,ChiRefGs_fGc,ChiRefInvGs_vGs); // keep
			BCoords_S[P]  = mm_Alloc_d(CBCM,CBT,CBT,NfnS[P], NfnGs,NfnGs,1.0,ChiRefGs_fS, ChiRefInvGs_vGs); // keep
			BCoords_Is[P] = mm_Alloc_d(CBCM,CBT,CBT,NfnIs[P],NfnGs,NfnGs,1.0,ChiRefGs_fIs,ChiRefInvGs_vGs); // keep
			BCoords_Ic[P] = mm_Alloc_d(CBCM,CBT,CBT,NfnIc[P],NfnGs,NfnGs,1.0,ChiRefGs_fIc,ChiRefInvGs_vGs); // keep

			free(rst_fGc);
			free(rst_fS);
			free(rst_fIs);
			free(rst_fIc);

			free(ChiRefGs_fGc);
			free(ChiRefGs_fS);
			free(ChiRefGs_fIs);
			free(ChiRefGs_fIc);
		}
		free(ChiRefInvGs_vGs);
	}

	free(one);

	BCoords_dEm1->Nve   = Nve;
	BCoords_dEm1->NfnG2 = NfnG2;
	BCoords_dEm1->NfnGc = NfnGc;
	BCoords_dEm1->NfnS  = NfnS;
	BCoords_dEm1->NfnIs = NfnIs;
	BCoords_dEm1->NfnIc = NfnIc;
	BCoords_dEm1->w_fIs = w_fIs;
	BCoords_dEm1->w_fIc = w_fIc;
	BCoords_dEm1->BCoords_G2 = BCoords_G2;
	BCoords_dEm1->BCoords_Gc = BCoords_Gc;
	BCoords_dEm1->BCoords_S  = BCoords_S;
	BCoords_dEm1->BCoords_Is = BCoords_Is;
	BCoords_dEm1->BCoords_Ic = BCoords_Ic;

	return BCoords_dEm1;
}

struct S_BCOORDS *get_BCoords_dEm2(const struct S_ELEMENT *ELEMENT)
{
	// Initialize DB Parameters
	char         **NodeTypeG = DB.NodeTypeG;
	unsigned int NP          = DB.NP,
	             PMax        = DB.PMax,
	             *PGc        = DB.PGc;

	// Standard datatypes
	unsigned int dE, Nbf, P, NenGs, Nve, EType, Eclass, *NenG2, *NenGc, dummy_ui, *dummyPtr_ui;
	double       **BCoords_G2, **BCoords_Gc, *rst_vGs, *rst_eG2, *rst_eGc,
	             *IGs, *ChiRefGs_vGs, *ChiRefInvGs_vGs, *ChiRefGs_eG2, *ChiRefGs_eGc, *dummyPtr_d[2];

	struct S_BCOORDS *BCoords_dEm2;
	struct S_ELEMENT *ELEMENT_E;

	// Function pointers
	cubature_tdef   cubature;
	basis_tdef      basis;

	// Note: This function is only needed for 3D ELEMENTs
	if (ELEMENT->d != DMAX)
		return NULL;

	BCoords_dEm2 = malloc(     sizeof *BCoords_dEm2); // keep
	NenG2        = calloc(NP , sizeof *NenG2);        // keep
	NenGc        = malloc(NP * sizeof *NenGc);        // keep
	BCoords_G2   = malloc(NP * sizeof *BCoords_G2);   // keep
	BCoords_Gc   = malloc(NP * sizeof *BCoords_Gc);   // keep

	ELEMENT_E = get_ELEMENT_type(LINE);

	EType = ELEMENT_E->type;

	select_functions_cubature(&cubature,EType);
	select_functions_basis(&basis,EType);
	Eclass = get_Eclass(EType);

	dE  = ELEMENT_E->d;
	Nve = ELEMENT_E->Nve;

	// It is important to use the nodes corresponding to the VeE ordering
	rst_vGs = get_rst_vV(ELEMENT_E); // free
	cubature(&dummyPtr_d[0],&dummyPtr_d[1],&dummyPtr_ui,&NenGs,&dummy_ui,0,1,dE,NodeTypeG[Eclass]); // free
	free(dummyPtr_ui);
	free(dummyPtr_d[0]);

	IGs             = identity_d(NenGs);                       // free
	ChiRefGs_vGs    = basis(1,rst_vGs,NenGs,&Nbf,dE);          // free
	ChiRefInvGs_vGs = inverse_d(NenGs,NenGs,ChiRefGs_vGs,IGs); // free

	free(rst_vGs);

	free(IGs);
	free(ChiRefGs_vGs);

	for (P = 0; P <= PMax; P++) {
		cubature(&rst_eGc,&dummyPtr_d[0],&dummyPtr_ui,&NenGc[P],&dummy_ui,0,PGc[P],   dE,NodeTypeG[Eclass]); free(dummyPtr_ui); // free

		ChiRefGs_eGc = basis(1,rst_eGc,NenGc[P],&Nbf,dE); // free

		if (P == 2) {
			cubature(&rst_eG2,&dummyPtr_d[0],&dummyPtr_ui,&NenG2[P],&dummy_ui,0,P,dE,NodeTypeG[Eclass]); free(dummyPtr_ui); // free
			ChiRefGs_eG2 = basis(1,rst_eG2,NenG2[P],&Nbf,dE); // free
			BCoords_G2[P] = mm_Alloc_d(CBCM,CBT,CBT,NenG2[P],NenGs,NenGs,1.0,ChiRefGs_eG2,ChiRefInvGs_vGs); // keep

			free(rst_eG2);
			free(ChiRefGs_eG2);
		}
		BCoords_Gc[P] = mm_Alloc_d(CBCM,CBT,CBT,NenGc[P],NenGs,NenGs,1.0,ChiRefGs_eGc,ChiRefInvGs_vGs); // keep

		free(rst_eGc);
		free(ChiRefGs_eGc);
	}
	free(ChiRefInvGs_vGs);

	BCoords_dEm2->Nve   = Nve;
	BCoords_dEm2->NenG2 = NenG2;
	BCoords_dEm2->NenGc = NenGc;
	BCoords_dEm2->BCoords_G2 = BCoords_G2;
	BCoords_dEm2->BCoords_Gc = BCoords_Gc;

	return BCoords_dEm2;
}

static void initialize_VeVref(const unsigned int EType, double *VeVref)
{
	// Initialize DB Parameters
	unsigned int TETrefineType = DB.TETrefineType;

	// Standard datatypes
	unsigned int i, j, k, Nve, Nvref, *Nvve;

	struct S_ELEMENT *ELEMENT;

	ELEMENT = get_ELEMENT_type(EType);

	Nve   = ELEMENT->Nve;
	Nvref = ELEMENT->Nvref;
	Nvve  = ELEMENT->Nvve;

	switch (EType) {
	case LINE: {
		double Barycentric_coords[3][2] = {{ 1.0 , 0.0 },
		                                   { 0.0 , 1.0 },
		                                   { 0.5 , 0.5 }};
		unsigned int Vh_coord_Ind[3][2] = {{ 0 , 1 },
		                                   { 0 , 2 },
		                                   { 2 , 1 }};

		for (k = 0; k < Nvref; k++) {
		for (i = 0; i < Nve; i++) {
		for (j = 0; j < Nve; j++) {
			*VeVref++ = Barycentric_coords[Vh_coord_Ind[k][i]][j];
		}}}
		break;
	} case TRI: {
		double Barycentric_coords[6][3] = {{ 1.0 , 0.0 , 0.0 },
		                                   { 0.0 , 1.0 , 0.0 },
		                                   { 0.0 , 0.0 , 1.0 },
		                                   { 0.5 , 0.5 , 0.0 },
		                                   { 0.5 , 0.0 , 0.5 },
		                                   { 0.0 , 0.5 , 0.5 }};
		unsigned int Vh_coord_Ind[5][3] = {{ 0 , 1 , 2 },
		                                   { 0 , 3 , 4 },
		                                   { 3 , 1 , 5 },
		                                   { 4 , 5 , 2 },
		                                   { 5 , 4 , 3 }};

		for (k = 0; k < Nvref; k++) {
		for (i = 0; i < Nve; i++) {
		for (j = 0; j < Nve; j++) {
			*VeVref++ = Barycentric_coords[Vh_coord_Ind[k][i]][j];
		}}}
		break;
	} case QUAD: {
		double Barycentric_coords[9][4] = {{ 1.0 , 0.0 , 0.0 , 0.0 },
		                                   { 0.0 , 1.0 , 0.0 , 0.0 },
		                                   { 0.0 , 0.0 , 1.0 , 0.0 },
		                                   { 0.0 , 0.0 , 0.0 , 1.0 },
		                                   { 0.5 , 0.5 , 0.0 , 0.0 },
		                                   { 0.5 , 0.0 , 0.5 , 0.0 },
		                                   { 0.0 , 0.5 , 0.0 , 0.5 },
		                                   { 0.0 , 0.0 , 0.5 , 0.5 },
		                                   { 0.25, 0.25, 0.25, 0.25}};
		unsigned int Vh_coord_Ind[9][4] = {{ 0 , 1 , 2 , 3 },
		                                   { 0 , 4 , 5 , 8 },
		                                   { 4 , 1 , 8 , 6 },
		                                   { 5 , 8 , 2 , 7 },
		                                   { 8 , 6 , 7 , 3 },
		                                   { 0 , 4 , 2 , 7 },
		                                   { 4 , 1 , 7 , 3 },
		                                   { 0 , 1 , 5 , 6 },
		                                   { 5 , 6 , 2 , 3 }};

		for (k = 0; k < Nvref; k++) {
		for (i = 0; i < Nve; i++) {
		for (j = 0; j < Nve; j++) {
			*VeVref++ = Barycentric_coords[Vh_coord_Ind[k][i]][j];
		}}}
		break;
	} case TET: {
		double Barycentric_coords[11][4] = {{ 1.0 , 0.0 , 0.0 , 0.0 },
		                                    { 0.0 , 1.0 , 0.0 , 0.0 },
		                                    { 0.0 , 0.0 , 1.0 , 0.0 },
		                                    { 0.0 , 0.0 , 0.0 , 1.0 },
		                                    { 0.5 , 0.5 , 0.0 , 0.0 },
		                                    { 0.5 , 0.0 , 0.5 , 0.0 },
		                                    { 0.0 , 0.5 , 0.5 , 0.0 },
		                                    { 0.5 , 0.0 , 0.0 , 0.5 },
		                                    { 0.0 , 0.5 , 0.0 , 0.5 },
		                                    { 0.0 , 0.0 , 0.5 , 0.5 },
		                                    { 0.25, 0.25, 0.25, 0.25}};
		if (TETrefineType == TET8) {
			unsigned int Vh_coord_Ind[9][4] = {{ 0 , 1 , 2 , 3 },
			                                   { 0 , 4 , 5 , 7 },
			                                   { 4 , 1 , 6 , 8 },
			                                   { 5 , 6 , 2 , 9 },
			                                   { 7 , 8 , 9 , 3 },
			                                   { 4 , 9 , 8 , 6 },
			                                   { 9 , 4 , 7 , 5 },
			                                   { 8 , 7 , 9 , 4 },
			                                   { 6 , 5 , 4 , 9 }};

			for (k = 0; k < Nvref; k++) {
			for (i = 0; i < Nvve[k]; i++) {
			for (j = 0; j < Nve; j++) {
				*VeVref++ = Barycentric_coords[Vh_coord_Ind[k][i]][j];
			}}}
		} else if (TETrefineType == TET12) {
			unsigned int Vh_coord_Ind[13][4] = {{ 0 , 1 , 2 , 3 },
			                                    { 0 , 4 , 5 , 7 },
			                                    { 4 , 1 , 6 , 8 },
			                                    { 5 , 6 , 2 , 9 },
			                                    { 7 , 8 , 9 , 3 },
			                                    { 10, 9 , 8 , 6 },
			                                    { 9 , 10, 7 , 5 },
			                                    { 8 , 7 , 10, 4 },
			                                    { 6 , 5 , 4 , 10},
			                                    { 10, 4 , 7 , 5 },
			                                    { 4 , 10, 8 , 6 },
			                                    { 6 , 5 , 10, 9 },
			                                    { 8 , 7 , 9 , 10}};

			for (k = 0; k < Nvref; k++) {
			for (i = 0; i < Nvve[k]; i++) {
			for (j = 0; j < Nve; j++) {
				*VeVref++ = Barycentric_coords[Vh_coord_Ind[k][i]][j];
			}}}
		} else if (TETrefineType == TET6) {
			unsigned int Vh_coord_Ind[7][5] = {{ 0 , 1 , 2 , 3 , -1},
			                                   { 0 , 4 , 5 , 7 , -1},
			                                   { 4 , 1 , 6 , 8 , -1},
			                                   { 5 , 6 , 2 , 9 , -1},
			                                   { 7 , 8 , 9 , 3 , -1},
			                                   { 8 , 7 , 6 , 5 , 4 },
			                                   { 7 , 8 , 5 , 6 , 9 }};

			for (k = 0; k < Nvref; k++) {
			for (i = 0; i < Nvve[k]; i++) {
			for (j = 0; j < Nve; j++) {
				*VeVref++ = Barycentric_coords[Vh_coord_Ind[k][i]][j];
			}}}
		} else {
			printf("Error: Unsupported.\n"), EXIT_MSG;
		}
		break;
	} case PYR: {
		double Barycentric_coords[14][5] = {{ 1.0 , 0.0 , 0.0 , 0.0 , 0.0 },
		                                    { 0.0 , 1.0 , 0.0 , 0.0 , 0.0 },
		                                    { 0.0 , 0.0 , 1.0 , 0.0 , 0.0 },
		                                    { 0.0 , 0.0 , 0.0 , 1.0 , 0.0 },
		                                    { 0.0 , 0.0 , 0.0 , 0.0 , 1.0 },
		                                    { 0.5 , 0.5 , 0.0 , 0.0 , 0.0 },
		                                    { 0.5 , 0.0 , 0.5 , 0.0 , 0.0 },
		                                    { 0.25, 0.25, 0.25, 0.25, 0.0 },
		                                    { 0.0 , 0.5 , 0.0 , 0.5 , 0.0 },
		                                    { 0.0 , 0.0 , 0.5 , 0.5 , 0.0 },
		                                    { 0.5 , 0.0 , 0.0 , 0.0 , 0.5 },
		                                    { 0.0 , 0.5 , 0.0 , 0.0 , 0.5 },
		                                    { 0.0 , 0.0 , 0.5 , 0.0 , 0.5 },
		                                    { 0.0 , 0.0 , 0.0 , 0.5 , 0.5 }};
		unsigned int Vh_coord_Ind[11][5] = {{ 0 , 1 , 2 , 3 , 4 },
		                                    { 0 , 5 , 6 , 7 , 10},
		                                    { 5 , 1 , 7 , 8 , 11},
		                                    { 6 , 7 , 2 , 9 , 12},
		                                    { 7 , 8 , 9 , 3 , 13},
		                                    { 6 , 7 , 12, 10, -1},
		                                    { 7 , 8 , 13, 11, -1},
		                                    { 11, 10, 7 , 5 , -1},
		                                    { 13, 12, 9 , 7 , -1},
		                                    { 11, 10, 13, 12, 7 },
		                                    { 10, 11, 12, 13, 4 }};

		for (k = 0; k < Nvref; k++) {
		for (i = 0; i < Nvve[k]; i++) {
		for (j = 0; j < Nve; j++) {
			*VeVref++ = Barycentric_coords[Vh_coord_Ind[k][i]][j];
		}}}
		break;
	} default:
		printf("Error: Unsupported EType.\n"), EXIT_MSG;
		break;
	}

}

void setup_ELEMENT_VeV(const unsigned int EType)
{
	/*
	 *	Comments:
	 *		For TET refinement, note that there are 3 options for the internal TET positioning. Here only one option is
	 *		supported but the others can be added if, for example, it is desired to refine TETs by splitting along the
	 *		longest physical diagonal. Similarly, if it is found that TETs are advantageous over PYRs, the top two PYRs
	 *		in the refined PYR can be split into four/eight TETs.
	 */

	// Initialize DB Parameters
	unsigned int TETrefineType = DB.TETrefineType;

	// Standard datatypes
	unsigned int i, j, jMax, VeVrefInd, Nve, *Nvve, Nvref, EcType, Nnodes;
	double       **VeV, *VeVref;

	struct S_ELEMENT *ELEMENT;

	ELEMENT = get_ELEMENT_type(EType);

	Nve   = ELEMENT->Nve;
	Nvve  = ELEMENT->Nvve;
	VeV   = ELEMENT->VeV;
	Nvref = ELEMENT->Nvref;

	if (EType == LINE || EType == TRI) {
		Nnodes = Nvref*Nve*Nve;
	} else if (EType == TET) {
		if      (TETrefineType == TET8)  Nnodes = 144;
		else if (TETrefineType == TET12) Nnodes = 208;
		else if (TETrefineType == TET6)  Nnodes = 120;
		else
			printf("Error: Unsupported.\n"), EXIT_MSG;
	} else if (EType == PYR) {
		Nnodes = 255;
	} else {
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}

	VeVref = malloc(Nnodes * sizeof *VeVref); // free

	if (EType == LINE || EType == TRI) {
		for (i = 0; i < Nvref; i++)
			Nvve[i] = Nve;
	} else {
		for (i = 0; i < Nvref; i++) {
			EcType = get_VOLUMEc_type(EType,i);
			if      (EcType == TET)  Nvve[i] = 4;
			else if (EcType == PYR)  Nvve[i] = 5;
			else
				printf("Error: Unsupported.\n"), EXIT_MSG;
		}
	}

	initialize_VeVref(EType,VeVref);

	VeVrefInd = 0;
	for (i = 0; i < Nvref; i++) {
		if (i)
			VeVrefInd += jMax;

		jMax = Nve*Nvve[i];
		VeV[i] = malloc(jMax * sizeof *VeV[i]); // keep
		for (j = 0; j < jMax; j++)
			VeV[i][j] = VeVref[VeVrefInd+j];
	}
	free(VeVref);
}

void setup_ELEMENT_VeF(const unsigned int EType)
{
	// Standard datatypes
	unsigned int i, j, k, l, iMax, jMax, kMax, lMax, iStep, Nve, *Nfve, *Nfref, Nf, *VeFcon, Nnodes;
	double       *VeF, *VeVref[3];

	struct S_ELEMENT *ELEMENT, *ELEMENT_F;

	// silence
	VeF = NULL;

	ELEMENT = get_ELEMENT_type(EType);
	Nve    = ELEMENT->Nve;
	Nfve   = ELEMENT->Nfve;
	Nf     = ELEMENT->Nf;
	VeFcon = ELEMENT->VeFcon;

	Nfref   = ELEMENT->Nfref;

	unsigned int size_VeF = 0;

	for (i = 0; i < 3; i++) {
		if      (i == 0) ELEMENT_F = get_ELEMENT_type(LINE);
		else if (i == 1) ELEMENT_F = get_ELEMENT_type(TRI);
		else if (i == 2) ELEMENT_F = get_ELEMENT_type(QUAD);

		Nnodes = (ELEMENT_F->Nvref)*pow(ELEMENT_F->Nve,2.0);
		VeVref[i] = malloc(Nnodes * sizeof *VeVref[i]); // free
		initialize_VeVref(ELEMENT_F->type,VeVref[i]);
	}

	switch(EType) {
	case LINE:
		VeF = calloc(Nve*Nfve[0]*Nfref[0]*Nf , sizeof *VeF); // free

		VeF[0] = 1.0; VeF[1] = 0.0;
		VeF[2] = 0.0; VeF[3] = 1.0;

		break;
	case TRI:
	case QUAD:
		VeF = calloc(Nve*Nfve[0]*Nfref[0]*Nf , sizeof *VeF); // free

		for (i = 0, iMax = Nf;               i < iMax; i++) {
		for (j = 0, jMax = Nfve[i]*Nfref[i]; j < jMax; j++) {
		for (k = 0, kMax = Nfve[i];          k < kMax; k++) {
			VeF[i*(Nfref[i]*Nfve[i]*Nve)+j*Nve+VeFcon[i*NFVEMAX+k]] = VeVref[0][j*kMax+k];
		}}}

		break;
	case TET:
		VeF = calloc(Nve*Nfve[0]*Nfref[0]*Nf , sizeof *VeF); // free

		for (i = 0, iMax = Nf;               i < iMax; i++) {
		for (j = 0, jMax = Nfve[i]*Nfref[i]; j < jMax; j++) {
		for (k = 0, kMax = Nfve[i];          k < kMax; k++) {
			VeF[i*(Nfref[i]*Nfve[i]*Nve)+j*Nve+VeFcon[i*NFVEMAX+k]] = VeVref[1][j*kMax+k];
		}}}

		break;
	case HEX:
		VeF = calloc(Nve*Nfve[0]*Nfref[0]*Nf , sizeof *VeF); // free

		for (i = 0, iMax = Nf;               i < iMax; i++) {
		for (j = 0, jMax = Nfve[i]*Nfref[i]; j < jMax; j++) {
		for (k = 0, kMax = Nfve[i];          k < kMax; k++) {
			VeF[i*(Nfref[i]*Nfve[i]*Nve)+j*Nve+VeFcon[i*NFVEMAX+k]] = VeVref[2][j*kMax+k];
		}}}

		break;
	case WEDGE:
	case PYR:
		for (i = 0; i < Nf; i++)
			size_VeF += Nfve[i]*Nfref[i];
		size_VeF *= Nve;

		VeF = calloc(size_VeF , sizeof *VeF); // free

		for (i = 0, iMax = Nf; i < iMax; i++) {
			for (j = iStep = 0; j < i; j++)
				iStep += Nfref[j]*Nfve[j];
			iStep *= Nve;

			for (j = 0, jMax = Nfve[i]*Nfref[i]; j < jMax; j++) {
			for (k = 0, kMax = Nfve[i];          k < kMax; k++) {
				if (Nfve[i] == 3)
					VeF[iStep+j*Nve+VeFcon[i*NFVEMAX+k]] = VeVref[1][j*kMax+k];
				else if (Nfve[i] == 4)
					VeF[iStep+j*Nve+VeFcon[i*NFVEMAX+k]] = VeVref[2][j*kMax+k];
			}
		}}

		break;
	default:
		// Already caught error above.
		break;
	}

	for (i = 0, iMax = Nf;       i < iMax; i++) {
		for (j = iStep = 0; j < i; j++)
			iStep += Nfref[j]*Nfve[j];
		iStep *= Nve;

		for (j = 0, jMax = Nfref[i]; j < jMax; j++) {
			ELEMENT->VeF[i*NFREFMAX+j] = malloc(Nfve[i]*Nve * sizeof **(ELEMENT->VeF)); // keep
			for (k = 0, kMax = Nfve[i];  k < kMax; k++) {
			for (l = 0, lMax = Nve;      l < lMax; l++) {
				ELEMENT->VeF[i*NFREFMAX+j][k*Nve+l] = VeF[iStep+j*(Nfve[i]*Nve)+k*(Nve)+l];
			}}
		}
	}
	free(VeF);

	for (i = 0; i < 3; i++)
		free(VeVref[i]);
}

void setup_ELEMENT_VeE(const unsigned int EType)
{
	// Standard datatypes
	unsigned int i, j, k, l, iMax, jMax, kMax, lMax, iStep, Nve, Neve, Ne, Neref, *VeEcon, Nnodes;
	double       *VeE, *VeVref;

	struct S_ELEMENT *ELEMENT, *ELEMENT_E;

	ELEMENT = get_ELEMENT_type(EType);

	if (ELEMENT->d != DMAX)
		printf("Error: Should not be entering.\n"), EXIT_MSG;

	Nve    = ELEMENT->Nve;
	Neve   = ELEMENT->Neve;
	Ne     = ELEMENT->Ne;
	Neref  = ELEMENT->Neref;
	VeEcon = ELEMENT->VeEcon;

	ELEMENT_E = get_ELEMENT_type(LINE);
	Nnodes = NREFMAXLINE*pow(ELEMENT_E->Nve,2.0);
	VeVref = malloc(Nnodes * sizeof *VeVref); // free
	initialize_VeVref(ELEMENT_E->type,VeVref);

	VeE = calloc(Nve*Neve*Neref*Ne , sizeof *VeE); // free

	for (i = 0, iMax = Ne;         i < iMax; i++) {
	for (j = 0, jMax = Neve*Neref; j < jMax; j++) {
	for (k = 0, kMax = Neve;       k < kMax; k++) {
		VeE[i*(Neref*Neve*Nve)+j*Nve+VeEcon[i*NEVEMAX+k]] = VeVref[j*kMax+k];
	}}}
	free(VeVref);

	for (i = 0, iMax = Ne;       i < iMax; i++) {
		for (j = iStep = 0; j < i; j++)
			iStep += Neref*Neve;
		iStep *= Nve;

		for (j = 0, jMax = Neref; j < jMax; j++) {
			ELEMENT->VeE[i*NEREFMAX+j] = malloc(Neve*Nve * sizeof **(ELEMENT->VeE)); // keep
			for (k = 0, kMax = Neve;  k < kMax; k++) {
			for (l = 0, lMax = Nve;   l < lMax; l++) {
				ELEMENT->VeE[i*NEREFMAX+j][k*Nve+l] = VeE[iStep+j*(Neve*Nve)+k*Nve+l];
			}}
		}
	}
	free(VeE);
}

double get_L2_scaling(const unsigned int EType, const unsigned int vref)
{
	/*
	 *	Purpose:
	 *		Returns the scaling for the L2 projection operator from fine to coarse.
	 */

	// Initialize DB Parameters
	unsigned int TETrefineType = DB.TETrefineType;

	if (vref == 0)
		return 1.0;

	switch (EType) {
	case LINE:
		return 0.5;
		break;
	case TRI:
		if (vref < 5)
			return 0.25;
		else
			printf("Error: Unsupported vref (%d) (anisotropic) for TRIs.\n",vref), EXIT_MSG;
		break;
	case TET:
		if (TETrefineType == TET8) {
			if (vref < 9)
				return 1.0/8.0;
			else
				printf("Error: Unsupported.\n"), EXIT_MSG;
		} else if (TETrefineType == TET12) {
			if (vref < 5)
				return 1.0/8.0;
			else if (vref < 13)
				return 1.0/16.0;
			else
				printf("Error: Unsupported.\n"), EXIT_MSG;
		} else if (TETrefineType == TET6) {
			if (vref < 5)
				return 1.0/8.0; // reference TET scaled by 0.5
			else if (vref < 7)
				return 1.0/8.0; // reference PYR scaled by 0.5
			else
				printf("Error: Unsupported.\n"), EXIT_MSG;
		}
		break;
	case PYR:
		if (vref < 5 || vref > 8)
			return 1.0/8.0; // reference PYR scaled by 0.5
		else
			return 1.0/8.0; // reference TET scaled by 0.5
		break;
	case QUAD:
	case HEX:
	case WEDGE:
	default:
		printf("Error: Unsupported EType (%d).\n",EType), EXIT_MSG;
		break;
	}

	// silence
	return -1e20;
}

void setup_ELEMENT_FACE_ordering(const unsigned int FType)
{
	// Returned operators
	unsigned int ***nOrd_fS, ***nOrd_fIs, ***nOrd_fIc;

	// Initialize DB Parameters
	const unsigned int d    = DB.d,
	                   PMax = DB.PMax;

	// Standard datatypes
	unsigned int P, NOrd, IndOrd;

	struct S_ELEMENT *ELEMENT;

	ELEMENT = get_ELEMENT_type(FType);

	if (FType != POINT) {
		// Initialize DB Parameters
		unsigned int **PIfs         = DB.PIfs,
		             **PIfc         = DB.PIfc;
		char         ***NodeTypeS   = DB.NodeTypeS,
		             ***NodeTypeIfs = DB.NodeTypeIfs,
		             ***NodeTypeIfc = DB.NodeTypeIfc;

		// Standard datatypes
		unsigned int dE, NfnS, NfnIs, NfnIc, NsS, NsIs, NsIc, Eclass,
		             *symms_fS, *symms_fIs, *symms_fIc;
		double       *rst_fS, *rst_fIs, *rst_fIc, *w;

		// Function pointers
		cubature_tdef   cubature;
		basis_tdef      basis;
		grad_basis_tdef grad_basis;

		select_functions(&basis,&grad_basis,&cubature,FType);

		// silence
		NfnS = NfnIs = NfnIc = NsS = NsIs = NsIc = 0;

		Eclass  = get_Eclass(FType);

		dE = d-1;

		if      (FType == LINE)  NOrd = 2;
		else if (FType == QUAD)  NOrd = 8;
		else if (FType == TRI)   NOrd = 6;
		else    printf("Error: Unsupported FType.\n"), EXIT_MSG;

		nOrd_fS  = ELEMENT->nOrd_fS;
		nOrd_fIs = ELEMENT->nOrd_fIs;
		nOrd_fIc = ELEMENT->nOrd_fIc;

		for (P = 0; P <= PMax; P++) {
			cubature(&rst_fS, &w,&symms_fS, &NfnS, &NsS, 0,P,              dE,NodeTypeS[P][Eclass]);   // free
			cubature(&rst_fIs,&w,&symms_fIs,&NfnIs,&NsIs,0,PIfs[P][Eclass],dE,NodeTypeIfs[P][Eclass]); // free
			cubature(&rst_fIc,&w,&symms_fIc,&NfnIc,&NsIc,0,PIfc[P][Eclass],dE,NodeTypeIfc[P][Eclass]); // free

			for (IndOrd = 0; IndOrd < NOrd; IndOrd++) {
				nOrd_fS[P][IndOrd]  = malloc(NfnS  * sizeof ***nOrd_fS);  //  keep
				nOrd_fIs[P][IndOrd] = malloc(NfnIs * sizeof ***nOrd_fIs); //  keep
				nOrd_fIc[P][IndOrd] = malloc(NfnIc * sizeof ***nOrd_fIc); //  keep

				get_face_ordering(d,IndOrd,FType,NfnS, NsS, symms_fS, rst_fS, nOrd_fS[P][IndOrd]);
				get_face_ordering(d,IndOrd,FType,NfnIs,NsIs,symms_fIs,rst_fIs,nOrd_fIs[P][IndOrd]);
				get_face_ordering(d,IndOrd,FType,NfnIc,NsIc,symms_fIc,rst_fIc,nOrd_fIc[P][IndOrd]);
			}
			free(rst_fS);
			free(rst_fIs);
			free(rst_fIc);
			free(symms_fS);
			free(symms_fIs);
			free(symms_fIc);
		}
	} else {
		NOrd = 1;

		nOrd_fS  = ELEMENT->nOrd_fS;
		nOrd_fIs = ELEMENT->nOrd_fIs;
		nOrd_fIc = ELEMENT->nOrd_fIc;

		for (P = 0; P <= PMax; P++) {
		for (IndOrd = 0; IndOrd < NOrd; IndOrd++) {
			nOrd_fS[P][IndOrd]  = malloc(1 * sizeof ***nOrd_fS);  //  keep
			nOrd_fIs[P][IndOrd] = malloc(1 * sizeof ***nOrd_fIs); //  keep
			nOrd_fIc[P][IndOrd] = malloc(1 * sizeof ***nOrd_fIc); //  keep

			get_face_ordering(d,0,0,0,0,NULL,NULL,nOrd_fS[P][IndOrd]);
			get_face_ordering(d,0,0,0,0,NULL,NULL,nOrd_fIs[P][IndOrd]);
			get_face_ordering(d,0,0,0,0,NULL,NULL,nOrd_fIc[P][IndOrd]);
		}}
	}
}
