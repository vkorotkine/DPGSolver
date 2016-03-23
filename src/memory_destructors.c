#include <stdlib.h>

#include "database.h"
#include "functions.h"

/*
 *	Purpose:
 *		Free memory of various structures.
 *
 *	Comments:
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
	// Mesh
	free(ELEMENT->Nfve);
	free(ELEMENT->VeC);
	free(ELEMENT->VeE);
	free(ELEMENT->VeF);

	// Operators
	
	// Normals
	free(ELEMENT->nr);

	// VOLUME Nodes
	array_free2_d(1    ,ELEMENT->xir_vGs);
	array_free2_d(DB.NP,ELEMENT->xir_vGc);
	array_free2_d(DB.NP,ELEMENT->xir_vCs);
	array_free2_d(DB.NP,ELEMENT->xir_vCc);
	array_free2_d(DB.NP,ELEMENT->xir_vJs);
	array_free2_d(DB.NP,ELEMENT->xir_vJc);
	array_free2_d(DB.NP,ELEMENT->xir_vS);
	array_free2_d(DB.NP,ELEMENT->xir_vF);
	array_free2_d(DB.NP,ELEMENT->xir_vFrs);
	array_free2_d(DB.NP,ELEMENT->xir_vFrc);
	array_free2_d(DB.NP,ELEMENT->xir_vIs);
	array_free2_d(DB.NP,ELEMENT->xir_vIc);
	array_free2_d(1    ,ELEMENT->xir_vP);

	array_free2_d(DB.NP,ELEMENT->WvIs);
	array_free2_d(DB.NP,ELEMENT->WvIc);

	array_free2_i(1    ,ELEMENT->Con_xir_vP);

	free(ELEMENT->NvnGs);
	free(ELEMENT->NvnGc);
	free(ELEMENT->NvnCs);
	free(ELEMENT->NvnCc);
	free(ELEMENT->NvnJs);
	free(ELEMENT->NvnJc);
	free(ELEMENT->NvnS);
	free(ELEMENT->NvnF);
	free(ELEMENT->NvnFrs);
	free(ELEMENT->NvnFrc);
	free(ELEMENT->NvnIs);
	free(ELEMENT->NvnIc);
	free(ELEMENT->NvnP);

	// FACET Nodes
	array_free3_d(DB.NP,ELEMENT->Nf,ELEMENT->xir_fGc);
	array_free3_d(DB.NP,ELEMENT->Nf,ELEMENT->xir_fIs);
	array_free3_d(DB.NP,ELEMENT->Nf,ELEMENT->xir_fIc);

	array_free2_d(DB.NP,ELEMENT->WfIs);
	array_free2_d(DB.NP,ELEMENT->WfIc);

	free(ELEMENT->NfnGc);
	free(ELEMENT->NfnIs);
	free(ELEMENT->NfnIc);

	free(ELEMENT);
}


