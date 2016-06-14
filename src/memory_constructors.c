// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>

#include "database.h"
#include "functions.h"
#include "parameters.h"

/*
 *	Purpose:
 *		Allocate memory for and initialize new structures.
 *
 *	Comments:
 *		All vGs arrays only need a range of 1 for the first dimension, as in I_vGs_vP. Consider removing superfluous
 *		memory storage in future (ToBeDeleted)
 *
 *	Notation:
 *
 *	References:
 *
 */

struct S_ELEMENT *New_ELEMENT(void)
{
	// Initialize DB Parameters
	unsigned int d    = DB.d,
	             PP   = DB.PP,
	             PMax = DB.PMax,
	             NP   = DB.NP;

	// Standard datatypes
	unsigned int P, Pb, PbMin, PbMax;

	struct S_ELEMENT *ELEMENT;

	ELEMENT = malloc(sizeof *ELEMENT); // free

	// Mesh
	ELEMENT->present = 0;
	ELEMENT->type    = 0;
	ELEMENT->d       = 0;
	ELEMENT->Nve     = 0;
	ELEMENT->Nf      = 0;

	ELEMENT->Nfve    = calloc(NFMAX         , sizeof *(ELEMENT->Nfve));    // free
	ELEMENT->VeCGmsh = calloc(NVEMAX        , sizeof *(ELEMENT->VeCGmsh)); // free
	ELEMENT->VeFcon  = calloc(NFMAX*NFVEMAX , sizeof *(ELEMENT->VeFcon));  // free

	// Operators

	// h-refinement related
	ELEMENT->Nfref   = calloc(NFMAX          , sizeof *(ELEMENT->Nfref));   // free
	ELEMENT->NfMixed = calloc(NFMIXEDMAX     , sizeof *(ELEMENT->NfMixed)); // free
	ELEMENT->VeF     = calloc(NFREFMAX*NFMAX , sizeof *(ELEMENT->VeF));     // free
	ELEMENT->Nvve    = calloc(NVREFMAX       , sizeof *(ELEMENT->Nvve));    // free
	ELEMENT->VeV     = calloc(NVREFMAX       , sizeof *(ELEMENT->VeV));     // free

	// Normals
	ELEMENT->nr = calloc(NFMAX*DMAX , sizeof *(ELEMENT->nr)); // free

	// Plotting
	ELEMENT->connectivity  = NULL;
	ELEMENT->connect_types = NULL;
	ELEMENT->NvnP          = 0;
	ELEMENT->connect_NE    = 0;

	// Operators
	ELEMENT->NvnGs = calloc(1  , sizeof *(ELEMENT->NvnGs)); // free
	ELEMENT->NvnGc = calloc(NP , sizeof *(ELEMENT->NvnGc)); // free
	ELEMENT->NvnCs = calloc(NP , sizeof *(ELEMENT->NvnCs)); // free
	ELEMENT->NvnCc = calloc(NP , sizeof *(ELEMENT->NvnCc)); // free
	ELEMENT->NvnIs = calloc(NP , sizeof *(ELEMENT->NvnIs)); // free
	ELEMENT->NvnIc = calloc(NP , sizeof *(ELEMENT->NvnIc)); // free
	ELEMENT->NvnS  = calloc(NP , sizeof *(ELEMENT->NvnS));  // free
	ELEMENT->NfnIs = calloc(NP , sizeof *(ELEMENT->NfnIs)); // free
	ELEMENT->NfnIc = calloc(NP , sizeof *(ELEMENT->NfnIc)); // free

	ELEMENT->w_vIs = calloc(NP , sizeof *(ELEMENT->w_vIs)); // free
	ELEMENT->w_vIc = calloc(NP , sizeof *(ELEMENT->w_vIc)); // free

	ELEMENT->ChiS_vP    = calloc(NP , sizeof *(ELEMENT->ChiS_vP));    // free
	ELEMENT->ChiS_vIs   = calloc(NP , sizeof *(ELEMENT->ChiS_vIs));   // free
	ELEMENT->ChiS_vIc   = calloc(NP , sizeof *(ELEMENT->ChiS_vIc));   // free
	ELEMENT->ChiInvS_vS = calloc(NP , sizeof *(ELEMENT->ChiInvS_vS)); // free

	ELEMENT->ICs = calloc(NP , sizeof *(ELEMENT->ICs)); // free
	ELEMENT->ICc = calloc(NP , sizeof *(ELEMENT->ICc)); // free

	ELEMENT->I_vGs_vP  = calloc(1  , sizeof *(ELEMENT->I_vGs_vP));  // free
	ELEMENT->I_vGs_vGc = calloc(1  , sizeof *(ELEMENT->I_vGs_vGc)); // free
	ELEMENT->I_vGs_vCs = calloc(NP , sizeof *(ELEMENT->I_vGs_vCs)); // free
	ELEMENT->I_vGs_vIs = calloc(NP , sizeof *(ELEMENT->I_vGs_vIs)); // free
	ELEMENT->I_vGs_vIc = calloc(NP , sizeof *(ELEMENT->I_vGs_vIc)); // free
	ELEMENT->I_vGs_vS  = calloc(NP , sizeof *(ELEMENT->I_vGs_vS));  // free
	ELEMENT->I_vGc_vP  = calloc(NP , sizeof *(ELEMENT->I_vGc_vP));  // free
	ELEMENT->I_vGc_vCc = calloc(NP , sizeof *(ELEMENT->I_vGc_vCc)); // free
	ELEMENT->I_vGc_vIs = calloc(NP , sizeof *(ELEMENT->I_vGc_vIs)); // free
	ELEMENT->I_vGc_vIc = calloc(NP , sizeof *(ELEMENT->I_vGc_vIc)); // free
	ELEMENT->I_vGc_vS  = calloc(NP , sizeof *(ELEMENT->I_vGc_vS));  // free
	ELEMENT->I_vCs_vIs = calloc(NP , sizeof *(ELEMENT->I_vCs_vIs)); // free
	ELEMENT->I_vCs_vIc = calloc(NP , sizeof *(ELEMENT->I_vCs_vIc)); // free
	ELEMENT->I_vCc_vIc = calloc(NP , sizeof *(ELEMENT->I_vCc_vIc)); // free
	ELEMENT->I_vCc_vIs = calloc(NP , sizeof *(ELEMENT->I_vCc_vIs)); // free

	ELEMENT->D_vGs_vCs = calloc(NP , sizeof *(ELEMENT->D_vGs_vCs)); // free
	ELEMENT->D_vGs_vIs = calloc(NP , sizeof *(ELEMENT->D_vGs_vIs)); // free
	ELEMENT->D_vCs_vCs = calloc(NP , sizeof *(ELEMENT->D_vCs_vCs)); // free
	ELEMENT->D_vGc_vCc = calloc(NP , sizeof *(ELEMENT->D_vGc_vCc)); // free
	ELEMENT->D_vGc_vIc = calloc(NP , sizeof *(ELEMENT->D_vGc_vIc)); // free
	ELEMENT->D_vCc_vCc = calloc(NP , sizeof *(ELEMENT->D_vCc_vCc)); // free

	ELEMENT->ChiS_fIs = calloc(NP , sizeof *(ELEMENT->ChiS_fIs)); // free
	ELEMENT->ChiS_fIc = calloc(NP , sizeof *(ELEMENT->ChiS_fIc)); // free

	ELEMENT->I_vGs_fIs = calloc(NP , sizeof *(ELEMENT->I_vGs_fIs)); // free
	ELEMENT->I_vGs_fIc = calloc(NP , sizeof *(ELEMENT->I_vGs_fIc)); // free
	ELEMENT->I_vGc_fIs = calloc(NP , sizeof *(ELEMENT->I_vGc_fIs)); // free
	ELEMENT->I_vGc_fIc = calloc(NP , sizeof *(ELEMENT->I_vGc_fIc)); // free
	ELEMENT->I_vCs_fIs = calloc(NP , sizeof *(ELEMENT->I_vCs_fIs)); // free
	ELEMENT->I_vCs_fIc = calloc(NP , sizeof *(ELEMENT->I_vCs_fIc)); // free
	ELEMENT->I_vCc_fIs = calloc(NP , sizeof *(ELEMENT->I_vCc_fIs)); // free
	ELEMENT->I_vCc_fIc = calloc(NP , sizeof *(ELEMENT->I_vCc_fIc)); // free

	ELEMENT->Is_Weak_VV = calloc(NP , sizeof *(ELEMENT->Is_Weak_VV)); // free
	ELEMENT->Ic_Weak_VV = calloc(NP , sizeof *(ELEMENT->Ic_Weak_VV)); // free
	ELEMENT->Is_Weak_FF = calloc(NP , sizeof *(ELEMENT->Is_Weak_FF)); // free
	ELEMENT->Ic_Weak_FF = calloc(NP , sizeof *(ELEMENT->Ic_Weak_FF)); // free
	ELEMENT->Ds_Weak_VV = calloc(NP , sizeof *(ELEMENT->Ds_Weak_VV)); // free
	ELEMENT->Dc_Weak_VV = calloc(NP , sizeof *(ELEMENT->Dc_Weak_VV)); // free

	ELEMENT->nOrd_fIs  = calloc(NP , sizeof *(ELEMENT->nOrd_fIs)); // free
	ELEMENT->nOrd_fIc  = calloc(NP , sizeof *(ELEMENT->nOrd_fIc)); // free


	ELEMENT->I_vGs_vP[0]  = calloc(NP , sizeof **(ELEMENT->I_vGs_vP));
	ELEMENT->I_vGs_vGc[0] = calloc(NP , sizeof **(ELEMENT->I_vGs_vGc));
	for (P = 0; P < NP; P++) {
		ELEMENT->NfnIs[P] = calloc(NESUBCMAX , sizeof **(ELEMENT->NfnIs));
		ELEMENT->NfnIc[P] = calloc(NESUBCMAX , sizeof **(ELEMENT->NfnIc));

		ELEMENT->ChiS_vP[P]    = calloc(NP , sizeof **(ELEMENT->ChiS_vP));
		ELEMENT->ChiS_vIs[P]   = calloc(NP , sizeof **(ELEMENT->ChiS_vIs));
		ELEMENT->ChiS_vIc[P]   = calloc(NP , sizeof **(ELEMENT->ChiS_vIc));
		ELEMENT->ChiInvS_vS[P] = calloc(NP , sizeof **(ELEMENT->ChiInvS_vS));

		ELEMENT->ICs[P] = calloc(NP , sizeof **(ELEMENT->ICs));
		ELEMENT->ICc[P] = calloc(NP , sizeof **(ELEMENT->ICc));

		ELEMENT->I_vGs_vCs[P] = calloc(NP , sizeof **(ELEMENT->I_vGs_vCs));
		ELEMENT->I_vGs_vIs[P] = calloc(NP , sizeof **(ELEMENT->I_vGs_vIs));
		ELEMENT->I_vGs_vIc[P] = calloc(NP , sizeof **(ELEMENT->I_vGs_vIc));
		ELEMENT->I_vGs_vS[P]  = calloc(NP , sizeof **(ELEMENT->I_vGs_vS));
		ELEMENT->I_vGc_vP[P]  = calloc(NP , sizeof **(ELEMENT->I_vGc_vP));
		ELEMENT->I_vGc_vCc[P] = calloc(NP , sizeof **(ELEMENT->I_vGc_vCc));
		ELEMENT->I_vGc_vIs[P] = calloc(NP , sizeof **(ELEMENT->I_vGc_vIs));
		ELEMENT->I_vGc_vIc[P] = calloc(NP , sizeof **(ELEMENT->I_vGc_vIc));
		ELEMENT->I_vGc_vS[P]  = calloc(NP , sizeof **(ELEMENT->I_vGc_vS));
		ELEMENT->I_vCs_vIs[P] = calloc(NP , sizeof **(ELEMENT->I_vCs_vIs));
		ELEMENT->I_vCs_vIc[P] = calloc(NP , sizeof **(ELEMENT->I_vCs_vIc));
		ELEMENT->I_vCc_vIs[P] = calloc(NP , sizeof **(ELEMENT->I_vCc_vIs));
		ELEMENT->I_vCc_vIc[P] = calloc(NP , sizeof **(ELEMENT->I_vCc_vIc));

		ELEMENT->D_vGs_vCs[P] = calloc(NP , sizeof **(ELEMENT->D_vGs_vCs));
		ELEMENT->D_vGs_vIs[P] = calloc(NP , sizeof **(ELEMENT->D_vGs_vIs));
		ELEMENT->D_vGc_vCc[P] = calloc(NP , sizeof **(ELEMENT->D_vGc_vCc));
		ELEMENT->D_vGc_vIc[P] = calloc(NP , sizeof **(ELEMENT->D_vGc_vIc));
		ELEMENT->D_vCs_vCs[P] = calloc(NP , sizeof **(ELEMENT->D_vCs_vCs));
		ELEMENT->D_vCc_vCc[P] = calloc(NP , sizeof **(ELEMENT->D_vCc_vCc));

		ELEMENT->ChiS_fIs[P] = calloc(NP , sizeof **(ELEMENT->ChiS_fIs));
		ELEMENT->ChiS_fIc[P] = calloc(NP , sizeof **(ELEMENT->ChiS_fIc));

		ELEMENT->I_vGs_fIs[P] = calloc(NP , sizeof **(ELEMENT->I_vGs_fIs));
		ELEMENT->I_vGs_fIc[P] = calloc(NP , sizeof **(ELEMENT->I_vGs_fIc));
		ELEMENT->I_vGc_fIs[P] = calloc(NP , sizeof **(ELEMENT->I_vGc_fIs));
		ELEMENT->I_vGc_fIc[P] = calloc(NP , sizeof **(ELEMENT->I_vGc_fIc));
		ELEMENT->I_vCs_fIs[P] = calloc(NP , sizeof **(ELEMENT->I_vCs_fIs));
		ELEMENT->I_vCs_fIc[P] = calloc(NP , sizeof **(ELEMENT->I_vCs_fIc));
		ELEMENT->I_vCc_fIs[P] = calloc(NP , sizeof **(ELEMENT->I_vCc_fIs));
		ELEMENT->I_vCc_fIc[P] = calloc(NP , sizeof **(ELEMENT->I_vCc_fIc));

		ELEMENT->Is_Weak_VV[P] = calloc(NP , sizeof **(ELEMENT->Is_Weak_VV));
		ELEMENT->Ic_Weak_VV[P] = calloc(NP , sizeof **(ELEMENT->Ic_Weak_VV));
		ELEMENT->Is_Weak_FF[P] = calloc(NP , sizeof **(ELEMENT->Is_Weak_FF));
		ELEMENT->Ic_Weak_FF[P] = calloc(NP , sizeof **(ELEMENT->Ic_Weak_FF));
		ELEMENT->Ds_Weak_VV[P] = calloc(NP , sizeof **(ELEMENT->Ds_Weak_VV));
		ELEMENT->Dc_Weak_VV[P] = calloc(NP , sizeof **(ELEMENT->Dc_Weak_VV));

		if (P == PP) {
			ELEMENT->I_vGs_vP[0][PP] = calloc(1 , sizeof ***(ELEMENT->I_vGs_vP));
		}

		if      (P == 0)    PbMin = P,   PbMax = P+1;
		else if (P == PMax) PbMin = P-1, PbMax = PMax;
		else                PbMin = P-1, PbMax = P+1;
		for (Pb = PbMin; Pb <= PbMax; Pb++) {
			ELEMENT->ChiS_vIs[P][Pb]   = calloc(NVREFSFMAX , sizeof ***(ELEMENT->ChiS_vIs));
			ELEMENT->ChiS_vIc[P][Pb]   = calloc(NVREFSFMAX , sizeof ***(ELEMENT->ChiS_vIc));

			ELEMENT->I_vGs_vIs[P][Pb] = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vGs_vIs));
			ELEMENT->I_vGs_vIc[P][Pb] = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vGs_vIc));
			ELEMENT->I_vGc_vIs[P][Pb] = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vGc_vIs));
			ELEMENT->I_vGc_vIc[P][Pb] = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vGc_vIc));
			ELEMENT->I_vCs_vIs[P][Pb] = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vCs_vIs));
			ELEMENT->I_vCs_vIc[P][Pb] = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vCs_vIc));
			ELEMENT->I_vCc_vIs[P][Pb] = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vCc_vIs));
			ELEMENT->I_vCc_vIc[P][Pb] = calloc(NVREFSFMAX , sizeof ***(ELEMENT->I_vCc_vIc));

			if (P == Pb) {
				ELEMENT->ChiS_vP[P][PP]    = calloc(1          , sizeof ***(ELEMENT->ChiS_vP));
				ELEMENT->ChiInvS_vS[P][Pb] = calloc(1          , sizeof ***(ELEMENT->ChiInvS_vS));

				ELEMENT->ICs[P][Pb] = calloc(1 , sizeof ***(ELEMENT->ICs));
				ELEMENT->ICc[P][Pb] = calloc(1 , sizeof ***(ELEMENT->ICc));

				ELEMENT->I_vGs_vGc[0][Pb] = calloc(1          , sizeof ***(ELEMENT->I_vGs_vGc));
				ELEMENT->I_vGs_vCs[P][Pb] = calloc(1          , sizeof ***(ELEMENT->I_vGs_vCs));
				ELEMENT->I_vGs_vS[P][Pb]  = calloc(1          , sizeof ***(ELEMENT->I_vGs_vS));
				ELEMENT->I_vGc_vP[P][PP]  = calloc(1          , sizeof ***(ELEMENT->I_vGc_vP));
				ELEMENT->I_vGc_vCc[P][Pb] = calloc(1          , sizeof ***(ELEMENT->I_vGc_vCc));
				ELEMENT->I_vGc_vS[P][Pb]  = calloc(1          , sizeof ***(ELEMENT->I_vGc_vS));

				ELEMENT->D_vGs_vCs[P][Pb] = calloc(1 , sizeof ***(ELEMENT->D_vGs_vCs));
				ELEMENT->D_vGs_vIs[P][Pb] = calloc(1 , sizeof ***(ELEMENT->D_vGs_vIs));
				ELEMENT->D_vGc_vCc[P][Pb] = calloc(1 , sizeof ***(ELEMENT->D_vGc_vCc));
				ELEMENT->D_vGc_vIc[P][Pb] = calloc(1 , sizeof ***(ELEMENT->D_vGc_vIc));
				ELEMENT->D_vCs_vCs[P][Pb] = calloc(1 , sizeof ***(ELEMENT->D_vCs_vCs));
				ELEMENT->D_vCc_vCc[P][Pb] = calloc(1 , sizeof ***(ELEMENT->D_vCc_vCc));

				ELEMENT->D_vGs_vCs[P][Pb][0] = calloc(d , sizeof ****(ELEMENT->D_vGs_vCs));
				ELEMENT->D_vGs_vIs[P][Pb][0] = calloc(d , sizeof ****(ELEMENT->D_vGs_vIs));
				ELEMENT->D_vGc_vCc[P][Pb][0] = calloc(d , sizeof ****(ELEMENT->D_vGc_vCc));
				ELEMENT->D_vGc_vIc[P][Pb][0] = calloc(d , sizeof ****(ELEMENT->D_vGc_vIc));
				ELEMENT->D_vCs_vCs[P][Pb][0] = calloc(d , sizeof ****(ELEMENT->D_vCs_vCs));
				ELEMENT->D_vCc_vCc[P][Pb][0] = calloc(d , sizeof ****(ELEMENT->D_vCc_vCc));

				ELEMENT->Is_Weak_VV[P][Pb] = calloc(1 , sizeof ***(ELEMENT->Is_Weak_VV));
				ELEMENT->Ic_Weak_VV[P][Pb] = calloc(1 , sizeof ***(ELEMENT->Ic_Weak_VV));
				ELEMENT->Ds_Weak_VV[P][Pb] = calloc(1 , sizeof ***(ELEMENT->Ds_Weak_VV));
				ELEMENT->Dc_Weak_VV[P][Pb] = calloc(1 , sizeof ***(ELEMENT->Dc_Weak_VV));

				ELEMENT->Ds_Weak_VV[P][Pb][0] = calloc(d , sizeof ****(ELEMENT->Ds_Weak_VV));
				ELEMENT->Dc_Weak_VV[P][Pb][0] = calloc(d , sizeof ****(ELEMENT->Dc_Weak_VV));
			}

			ELEMENT->ChiS_fIs[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->ChiS_fIs));
			ELEMENT->ChiS_fIc[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->ChiS_fIc));

			ELEMENT->I_vGs_fIs[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_vGs_fIs));
			ELEMENT->I_vGs_fIc[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_vGs_fIc));
			ELEMENT->I_vGc_fIs[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_vGc_fIs));
			ELEMENT->I_vGc_fIc[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_vGc_fIc));
			ELEMENT->I_vCs_fIs[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_vCs_fIs));
			ELEMENT->I_vCs_fIc[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_vCs_fIc));
			ELEMENT->I_vCc_fIs[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_vCc_fIs));
			ELEMENT->I_vCc_fIc[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->I_vCc_fIc));

			ELEMENT->Is_Weak_FF[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->Is_Weak_FF));
			ELEMENT->Ic_Weak_FF[P][Pb] = calloc(NFREFMAX*NFMAX , sizeof ***(ELEMENT->Ic_Weak_FF));
		}

		ELEMENT->nOrd_fIs[P] = calloc(NFORDMAX, sizeof **(ELEMENT->nOrd_fIs));
		ELEMENT->nOrd_fIc[P] = calloc(NFORDMAX, sizeof **(ELEMENT->nOrd_fIc));
	}


	ELEMENT->next = NULL;
	ELEMENT->ELEMENTclass  = calloc(NESUBCMAX  , sizeof *(ELEMENT->ELEMENTclass)); // free
	ELEMENT->ELEMENT_FACET = calloc(NFMIXEDMAX , sizeof *(ELEMENT->ELEMENTclass)); // free

	return ELEMENT;
}


struct S_VOLUME *New_VOLUME(void)
{
	//unsigned int NP = DB.NP;

	struct S_VOLUME *VOLUME;
	VOLUME = malloc(sizeof *VOLUME); // free

	// Structures
	VOLUME->indexl = 0;
	VOLUME->indexg = 0;
	VOLUME->P      = 0;
	VOLUME->type   = 0;
	VOLUME->Eclass = 0;
	VOLUME->update = 0;
	VOLUME->curved = 0;

//	VOLUME->Vneigh = malloc(6*9 * sizeof *(VOLUME->Vneigh)); // tbd
//	VOLUME->Fneigh = malloc(6*9 * sizeof *(VOLUME->Fneigh)); // tbd

	VOLUME->XYZ_vC = NULL; // free

	// Geometry
	VOLUME->NvnG  = 0;
	VOLUME->XYZ_S = NULL; // free
	VOLUME->XYZ   = NULL; // free

	VOLUME->detJV_vI = NULL; // free
	VOLUME->C_vC     = NULL; // free (in setup_normals)
	VOLUME->C_vI     = NULL; // free
	VOLUME->C_vf     = calloc(NFMAX , sizeof *(VOLUME->C_vf)); // free

	// Initialization
	VOLUME->NvnS = 0;
	VOLUME->What = NULL; // free
	VOLUME->RES  = NULL; // free

	// Solving
	VOLUME->RHS       = NULL; // tbd
	VOLUME->wdetJV_vI = NULL; // free
	VOLUME->MInv      = NULL; // free

	// hp adaptivity
//	VOLUME->minRES = 0.0;
//	VOLUME->maxRES = 0.0;

	VOLUME->Vadapt     = 0;
	VOLUME->adapt_type = 0;
	VOLUME->PNew       = 0;

	// structs
	VOLUME->next    = NULL;
	VOLUME->grpnext = NULL;

	return VOLUME;
}

struct S_FACET *New_FACET(void)
{
	struct S_FACET *FACET;
	FACET = malloc(sizeof *FACET); // free

	// Structures
	FACET->indexg = 0;
	FACET->P      = 0;
	FACET->type   = 0;
	FACET->VfIn   = 0;
	FACET->VfOut  = 0;
	FACET->BC     = 0;

	FACET->VIn   = NULL; // free (in memory_destructor_V)
	FACET->VOut  = NULL; // free (in memory_destructor_V)
	FACET->VfIn  = 0;
	FACET->VfOut = 0;

	FACET->IndOrdInOut = 0;
	FACET->IndOrdOutIn = 0;

	// Geometry
	FACET->curved  = 0;
	FACET->typeInt = 0;

	FACET->n_fI     = NULL; // free
	FACET->detJF_fI = NULL; // free

	// Solving
	FACET->RHSIn  = NULL; // tbd
	FACET->RHSOut = NULL; // tbd

	FACET->next = NULL;

	return FACET;
}
