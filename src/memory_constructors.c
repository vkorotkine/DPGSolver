// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>

#include "database.h"
#include "functions.h"

/*
 *	Purpose:
 *		Allocate memory and initialize new structures.
 *
 *	Comments:
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

	ELEMENT->Nfve    = calloc(6       , sizeof *(ELEMENT->Nfve));    // free
	ELEMENT->VeCGmsh = calloc(8       , sizeof *(ELEMENT->VeCGmsh)); // free
	ELEMENT->VeFcon  = calloc(6*4     , sizeof *(ELEMENT->VeFcon));  // free

	// Operators

	// h-refinement related
	ELEMENT->Nfref   = calloc(8       , sizeof *(ELEMENT->Nfref));   // free
	ELEMENT->NfMixed = calloc(2       , sizeof *(ELEMENT->NfMixed)); // free
	ELEMENT->VeF     = calloc(8*4*9*6 , sizeof *(ELEMENT->VeF));     // free
	// 8 VOLUME vertices (HEX), 4 FACET vertices (QUAD),
	// 9 potential refinements (HEX), 6 FACETs (HEX)

	// Normals
	ELEMENT->nr = calloc(6*3 , sizeof *(ELEMENT->nr)); // free

	// Plotting
	ELEMENT->connectivity  = NULL;
	ELEMENT->connect_types = NULL;
	ELEMENT->NvnP          = 0;
	ELEMENT->connect_NE    = 0;

	// Operators
	ELEMENT->NvnGs  = calloc(1  , sizeof *(ELEMENT->NvnGs));  // free
	ELEMENT->NvnGc  = calloc(NP , sizeof *(ELEMENT->NvnGc));  // free
	ELEMENT->NvnCs  = calloc(NP , sizeof *(ELEMENT->NvnCs));  // free
	ELEMENT->NvnCc  = calloc(NP , sizeof *(ELEMENT->NvnCc));  // free
	ELEMENT->NvnIs  = calloc(NP , sizeof *(ELEMENT->NvnIs));  // free
	ELEMENT->NvnIc  = calloc(NP , sizeof *(ELEMENT->NvnIc));  // free
	ELEMENT->NvnS   = calloc(NP , sizeof *(ELEMENT->NvnS));   // free

	ELEMENT->NfnIs  = calloc(NP , sizeof *(ELEMENT->NfnIs));  // free
	ELEMENT->NfnIc  = calloc(NP , sizeof *(ELEMENT->NfnIc));  // free

	ELEMENT->ChiInvS_vS = calloc(NP , sizeof *(ELEMENT->ChiInvS_vS)); // free
	ELEMENT->ChiS_vP    = calloc(NP , sizeof *(ELEMENT->ChiS_vP));    // free
	ELEMENT->ChiS_vIs   = calloc(NP , sizeof *(ELEMENT->ChiS_vIs));   // free
	ELEMENT->ChiS_vIc   = calloc(NP , sizeof *(ELEMENT->ChiS_vIc));   // free

	ELEMENT->ChiS_fIs   = calloc(NP , sizeof *(ELEMENT->ChiS_fIs)); // free
	ELEMENT->ChiS_fIc   = calloc(NP , sizeof *(ELEMENT->ChiS_fIc)); // free

	ELEMENT->ICs       = calloc(NP , sizeof *(ELEMENT->ICs));       // free
	ELEMENT->ICc       = calloc(NP , sizeof *(ELEMENT->ICc));       // free
	ELEMENT->I_vGs_vP  = calloc(1  , sizeof *(ELEMENT->I_vGs_vP));  // free
	ELEMENT->I_vGs_vGc = calloc(NP , sizeof *(ELEMENT->I_vGs_vGc)); // free
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

	ELEMENT->I_vGs_fIs = calloc(NP , sizeof *(ELEMENT->I_vGs_fIs)); // free
	ELEMENT->I_vGs_fIc = calloc(NP , sizeof *(ELEMENT->I_vGs_fIc)); // free
	ELEMENT->I_vGc_fIs = calloc(NP , sizeof *(ELEMENT->I_vGc_fIs)); // free
	ELEMENT->I_vGc_fIc = calloc(NP , sizeof *(ELEMENT->I_vGc_fIc)); // free
	ELEMENT->I_vCs_fIs = calloc(NP , sizeof *(ELEMENT->I_vCs_fIs)); // free
	ELEMENT->I_vCs_fIc = calloc(NP , sizeof *(ELEMENT->I_vCs_fIc)); // free
	ELEMENT->I_vCc_fIs = calloc(NP , sizeof *(ELEMENT->I_vCc_fIs)); // free
	ELEMENT->I_vCc_fIc = calloc(NP , sizeof *(ELEMENT->I_vCc_fIc)); // free

	ELEMENT->D_vGs_vCs = calloc(NP , sizeof *(ELEMENT->D_vGs_vCs)); // free
	ELEMENT->D_vGs_vIs = calloc(NP , sizeof *(ELEMENT->D_vGs_vIs)); // free
	ELEMENT->D_vCs_vCs = calloc(NP , sizeof *(ELEMENT->D_vCs_vCs)); // free
	ELEMENT->D_vGc_vCc = calloc(NP , sizeof *(ELEMENT->D_vGc_vCc)); // free
	ELEMENT->D_vGc_vIc = calloc(NP , sizeof *(ELEMENT->D_vGc_vIc)); // free
	ELEMENT->D_vCc_vCc = calloc(NP , sizeof *(ELEMENT->D_vCc_vCc)); // free

	ELEMENT->Ds_Weak   = calloc(NP , sizeof *(ELEMENT->Ds_Weak));   // free
	ELEMENT->Dc_Weak   = calloc(NP , sizeof *(ELEMENT->Dc_Weak));   // free

	for (P = 0; P < NP; P++) {
		ELEMENT->NfnIs[P]     = calloc(2 , sizeof **(ELEMENT->NfnIs));
		ELEMENT->NfnIc[P]     = calloc(2 , sizeof **(ELEMENT->NfnIc));

		ELEMENT->D_vGs_vCs[P] = calloc(d , sizeof **(ELEMENT->D_vGs_vCs));
		ELEMENT->D_vGs_vIs[P] = calloc(d , sizeof **(ELEMENT->D_vGs_vIs));
		ELEMENT->D_vCs_vCs[P] = calloc(d , sizeof **(ELEMENT->D_vCs_vCs));
		ELEMENT->D_vGc_vCc[P] = calloc(d , sizeof **(ELEMENT->D_vGc_vCc));
		ELEMENT->D_vGc_vIc[P] = calloc(d , sizeof **(ELEMENT->D_vGc_vIc));
		ELEMENT->D_vCc_vCc[P] = calloc(d , sizeof **(ELEMENT->D_vCc_vCc));

		ELEMENT->Ds_Weak[P]   = calloc(d , sizeof **(ELEMENT->Ds_Weak));
		ELEMENT->Dc_Weak[P]   = calloc(d , sizeof **(ELEMENT->Dc_Weak));

		ELEMENT->ChiS_fIs[P]  = calloc(NP, sizeof **(ELEMENT->ChiS_fIs));
		ELEMENT->ChiS_fIc[P]  = calloc(NP, sizeof **(ELEMENT->ChiS_fIc));

		ELEMENT->I_vGs_fIs[P] = calloc(NP, sizeof **(ELEMENT->I_vGs_fIs));
		ELEMENT->I_vGs_fIc[P] = calloc(NP, sizeof **(ELEMENT->I_vGs_fIc));
		ELEMENT->I_vGc_fIs[P] = calloc(NP, sizeof **(ELEMENT->I_vGc_fIs));
		ELEMENT->I_vGc_fIc[P] = calloc(NP, sizeof **(ELEMENT->I_vGc_fIc));
		ELEMENT->I_vCs_fIs[P] = calloc(NP, sizeof **(ELEMENT->I_vCs_fIs));
		ELEMENT->I_vCs_fIc[P] = calloc(NP, sizeof **(ELEMENT->I_vCs_fIc));
		ELEMENT->I_vCc_fIs[P] = calloc(NP, sizeof **(ELEMENT->I_vCc_fIs));
		ELEMENT->I_vCc_fIc[P] = calloc(NP, sizeof **(ELEMENT->I_vCc_fIc));

		if      (P == 0)    PbMin = P,   PbMax = P+1;
		else if (P == PMax) PbMin = P-1, PbMax = PMax;
		else                PbMin = P-1, PbMax = P+1;

		for (Pb = PbMin; Pb <= PbMax; Pb++) {
			ELEMENT->ChiS_fIs[P][Pb]  = calloc(54, sizeof ***(ELEMENT->ChiS_fIs));
			ELEMENT->ChiS_fIc[P][Pb]  = calloc(54, sizeof ***(ELEMENT->ChiS_fIc));

			ELEMENT->I_vGs_fIs[P][Pb] = calloc(54, sizeof ***(ELEMENT->I_vGs_fIs));
			ELEMENT->I_vGs_fIc[P][Pb] = calloc(54, sizeof ***(ELEMENT->I_vGs_fIc));
			ELEMENT->I_vGc_fIs[P][Pb] = calloc(54, sizeof ***(ELEMENT->I_vGc_fIs));
			ELEMENT->I_vGc_fIc[P][Pb] = calloc(54, sizeof ***(ELEMENT->I_vGc_fIc));
			ELEMENT->I_vCs_fIs[P][Pb] = calloc(54, sizeof ***(ELEMENT->I_vCs_fIs));
			ELEMENT->I_vCs_fIc[P][Pb] = calloc(54, sizeof ***(ELEMENT->I_vCs_fIc));
			ELEMENT->I_vCc_fIs[P][Pb] = calloc(54, sizeof ***(ELEMENT->I_vCc_fIs));
			ELEMENT->I_vCc_fIc[P][Pb] = calloc(54, sizeof ***(ELEMENT->I_vCc_fIc));
		}
	}








	// VOLUME Nodes
/*
	ELEMENT->rst_vGs  = malloc(1  * sizeof *(ELEMENT->rst_vGs));  // free
	ELEMENT->rst_vGc  = malloc(NP * sizeof *(ELEMENT->rst_vGc));  // free
	ELEMENT->rst_vCs  = malloc(NP * sizeof *(ELEMENT->rst_vCs));  // free
	ELEMENT->rst_vCc  = malloc(NP * sizeof *(ELEMENT->rst_vCc));  // free
	ELEMENT->rst_vJs  = malloc(NP * sizeof *(ELEMENT->rst_vJs));  // free
	ELEMENT->rst_vJc  = malloc(NP * sizeof *(ELEMENT->rst_vJc));  // free
	ELEMENT->rst_vS   = malloc(NP * sizeof *(ELEMENT->rst_vS));   // free
	ELEMENT->rst_vF   = malloc(NP * sizeof *(ELEMENT->rst_vF));   // free
	ELEMENT->rst_vFrs = malloc(NP * sizeof *(ELEMENT->rst_vFrs)); // free
	ELEMENT->rst_vFrc = malloc(NP * sizeof *(ELEMENT->rst_vFrc)); // free
	ELEMENT->rst_vIs  = malloc(NP * sizeof *(ELEMENT->rst_vIs));  // free
	ELEMENT->rst_vIc  = malloc(NP * sizeof *(ELEMENT->rst_vIc));  // free

	ELEMENT->wvIs = malloc(NP * sizeof *(ELEMENT->wvIs)); // free
	ELEMENT->wvIc = malloc(NP * sizeof *(ELEMENT->wvIc)); // free

	ELEMENT->rst_vGs[0] = NULL;

	for (P = 0; P < NP; P++) {
		ELEMENT->rst_vGc[P]  = NULL;
		ELEMENT->rst_vCs[P]  = NULL;
		ELEMENT->rst_vCc[P]  = NULL;
		ELEMENT->rst_vJs[P]  = NULL;
		ELEMENT->rst_vJc[P]  = NULL;
		ELEMENT->rst_vS[P]   = NULL;
		ELEMENT->rst_vF[P]   = NULL;
		ELEMENT->rst_vFrs[P] = NULL;
		ELEMENT->rst_vFrc[P] = NULL;
		ELEMENT->rst_vIs[P]  = NULL;
		ELEMENT->rst_vIc[P]  = NULL;

		ELEMENT->wvIs[P] = NULL;
		ELEMENT->wvIc[P] = NULL;
	}

	// FACET Nodes
	ELEMENT->rst_fGc  = malloc(NP * sizeof *(ELEMENT->rst_fGc)); // free
	ELEMENT->rst_fIs  = malloc(NP * sizeof *(ELEMENT->rst_fIs)); // free
	ELEMENT->rst_fIc  = malloc(NP * sizeof *(ELEMENT->rst_fIc)); // free

	ELEMENT->wfIs  = malloc(NP * sizeof *(ELEMENT->wfIs)); // free
	ELEMENT->wfIc  = malloc(NP * sizeof *(ELEMENT->wfIc)); // free

	for (P = 0; P < NP; P++) {
		ELEMENT->rst_fGc[P] = NULL;
		ELEMENT->rst_fIs[P] = NULL;
		ELEMENT->rst_fIc[P] = NULL;

		ELEMENT->wfIs[P] = NULL;
		ELEMENT->wfIc[P] = NULL;
	}

	ELEMENT->NfnGc = malloc(NP * sizeof *(ELEMENT->NfnGc)); // free
	ELEMENT->NfnIs = malloc(NP * sizeof *(ELEMENT->NfnIs)); // free
	ELEMENT->NfnIc = malloc(NP * sizeof *(ELEMENT->NfnIc)); // free
*/


	ELEMENT->next = NULL;
	ELEMENT->ELEMENTclass = calloc(2 , sizeof *(ELEMENT->ELEMENTclass)); // free

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
	VOLUME->C_vf     = calloc(6 , sizeof *(VOLUME->C_vf)); // free

	// Initialization
	VOLUME->NvnS = 0;
	VOLUME->What = NULL; // free
	VOLUME->RES  = NULL; // free

	// Solving
	VOLUME->RHS = NULL; // tbd

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
	FACET->VfIn   = 0;
	FACET->VfOut  = 0;
	FACET->BC     = 0;

	FACET->VIn   = NULL; // free (in memory_destructor_V)
	FACET->VOut  = NULL; // free (in memory_destructor_V)
	FACET->VfIn  = 0;
	FACET->VfOut = 0;

	// Geometry
	FACET->curved  = 0;
	FACET->typeInt = 0;

	FACET->n = NULL; // free

	// Solving
	FACET->RHSIn  = NULL; // tbd
	FACET->RHSOut = NULL; // tbd

	FACET->next = NULL;

	return FACET;
}
