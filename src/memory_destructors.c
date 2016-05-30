// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>

#include "database.h"
#include "functions.h"
#include "parameters.h"

/*
 *	Purpose:
 *		Free memory of various structures.
 *
 *	Comments:
 *		ToBeDeleted: Likely need a destructor for elements and element classes which are different.
 *
 *	Notation:
 *
 *	Memory freed:
 *
 *	References:
 *
 */

void memory_destructor_E(struct S_ELEMENT *ELEMENT)
{
	// Initialize DB Parameters
	unsigned int d  = DB.d,
	             NP = DB.NP;

	// Mesh
	free(ELEMENT->Nfve);
	free(ELEMENT->VeCGmsh);
	free(ELEMENT->VeFcon);

	// h-refinement related
	free(ELEMENT->Nfref);
	free(ELEMENT->NfMixed);
	free(ELEMENT->VeF);

	// Normals
	free(ELEMENT->nr);

	// Plotting
	free(ELEMENT->connectivity);
	free(ELEMENT->connect_types);

	// Operators
	free(ELEMENT->NvnGs);
	free(ELEMENT->NvnGc);
	free(ELEMENT->NvnCs);
	free(ELEMENT->NvnCc);
	free(ELEMENT->NvnIs);
	free(ELEMENT->NvnIc);
	free(ELEMENT->NvnS);

	array_free2_ui(NP,ELEMENT->NfnIs);
	array_free2_ui(NP,ELEMENT->NfnIc);

	array_free2_d(NP,ELEMENT->ChiInvS_vS);
	array_free2_d(NP,ELEMENT->ChiS_vP);
	array_free2_d(NP,ELEMENT->ChiS_vIs);
	array_free2_d(NP,ELEMENT->ChiS_vIc);

	array_free2_d(NP,ELEMENT->ICs);
	array_free2_d(NP,ELEMENT->ICc);
	array_free2_d(1 ,ELEMENT->I_vGs_vP);
	array_free2_d(NP,ELEMENT->I_vGs_vGc);
	array_free2_d(NP,ELEMENT->I_vGs_vCs);
	array_free2_d(NP,ELEMENT->I_vGs_vIs);
	array_free2_d(NP,ELEMENT->I_vGs_vIc);
	array_free2_d(NP,ELEMENT->I_vGs_vS);
	array_free2_d(NP,ELEMENT->I_vGc_vP);
	array_free2_d(NP,ELEMENT->I_vGc_vCc);
	array_free2_d(NP,ELEMENT->I_vGc_vIs);
	array_free2_d(NP,ELEMENT->I_vGc_vIc);
	array_free2_d(NP,ELEMENT->I_vGc_vS);
	array_free2_d(NP,ELEMENT->I_vCs_vIs);
	array_free2_d(NP,ELEMENT->I_vCs_vIc);
	array_free2_d(NP,ELEMENT->I_vCc_vIs);
	array_free2_d(NP,ELEMENT->I_vCc_vIc);

	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->ChiS_fIs);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->ChiS_fIc);

	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->I_vGs_fIs);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->I_vGs_fIc);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->I_vGc_fIs);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->I_vGc_fIc);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->I_vCs_fIs);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->I_vCs_fIc);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->I_vCc_fIs);
	array_free4_d(NP,NP,NFREFMAX*NFMAX,ELEMENT->I_vCc_fIc);

	array_free3_d(NP,d,ELEMENT->D_vGs_vCs);
	array_free3_d(NP,d,ELEMENT->D_vGs_vIs);
	array_free3_d(NP,d,ELEMENT->D_vCs_vCs);
	array_free3_d(NP,d,ELEMENT->D_vGc_vCc);
	array_free3_d(NP,d,ELEMENT->D_vGc_vIc);
	array_free3_d(NP,d,ELEMENT->D_vCc_vCc);

	array_free3_d(NP,d,ELEMENT->Ds_Weak);
	array_free3_d(NP,d,ELEMENT->Dc_Weak);

	array_free3_ui(NP,NFORDMAX,ELEMENT->nOrd_fIs);
	array_free3_ui(NP,NFORDMAX,ELEMENT->nOrd_fIc);

	// VOLUME Nodes
/*
	array_free2_d(1 ,ELEMENT->rst_vGs);
	array_free2_d(NP,ELEMENT->rst_vGc);
	array_free2_d(NP,ELEMENT->rst_vCs);
	array_free2_d(NP,ELEMENT->rst_vCc);
	array_free2_d(NP,ELEMENT->rst_vJs);
	array_free2_d(NP,ELEMENT->rst_vJc);
	array_free2_d(NP,ELEMENT->rst_vS);
	array_free2_d(NP,ELEMENT->rst_vF);
	array_free2_d(NP,ELEMENT->rst_vFrs);
	array_free2_d(NP,ELEMENT->rst_vFrc);
	array_free2_d(NP,ELEMENT->rst_vIs);
	array_free2_d(NP,ELEMENT->rst_vIc);

	array_free2_d(NP,ELEMENT->wvIs);
	array_free2_d(NP,ELEMENT->wvIc);

	free(ELEMENT->NvnS);
	free(ELEMENT->NvnF);
	free(ELEMENT->NvnFrs);
	free(ELEMENT->NvnFrc);
	free(ELEMENT->NvnIs);
	free(ELEMENT->NvnIc);

	// FACET Nodes
	array_free3_d(NP,ELEMENT->Nf,ELEMENT->rst_fGc);
	array_free3_d(NP,ELEMENT->Nf,ELEMENT->rst_fIs);
	array_free3_d(NP,ELEMENT->Nf,ELEMENT->rst_fIc);

	array_free2_d(NP,ELEMENT->wfIs);
	array_free2_d(NP,ELEMENT->wfIc);

	free(ELEMENT->NfnGc);
	free(ELEMENT->NfnIs);
	free(ELEMENT->NfnIc);
*/
	free(ELEMENT->ELEMENTclass);
	free(ELEMENT->ELEMENT_FACET);

	free(ELEMENT);
}

void memory_destructor_V(struct S_VOLUME *VOLUME)
{
//	int NP = DB.NP;

	// Structures
	free(VOLUME->XYZ_vC);

	// Geometry
	free(VOLUME->XYZ_S);
	free(VOLUME->XYZ);

	free(VOLUME->detJV_vI);
	free(VOLUME->C_vI);
	array_free2_d(NFMAX,VOLUME->C_vf);

	// Initialization
	free(VOLUME->What);
	free(VOLUME->RES);

	// Solving

	free(VOLUME);
}

void memory_destructor_F(struct S_FACET *FACET)
{
	free(FACET->n);

	free(FACET);
}
