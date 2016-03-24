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
	int NP = DB.NP;

	int i, iMax, P;

	struct S_ELEMENT *ELEMENT;
	ELEMENT = malloc(sizeof *ELEMENT); // free

	// Mesh
	ELEMENT->present = 0;
	ELEMENT->type    = -1;
	ELEMENT->d       = -1;
	ELEMENT->Nve     = 0;
	ELEMENT->Nf      = 0;

	ELEMENT->Nfve = malloc(2    * sizeof *(ELEMENT->Nfve)); // free
	ELEMENT->VeC  = malloc(8    * sizeof *(ELEMENT->VeC));  // free
	ELEMENT->VeE  = malloc(12*2 * sizeof *(ELEMENT->VeE));  // free
	ELEMENT->VeF  = malloc(6*4  * sizeof *(ELEMENT->VeF));  // free

	for (i = 0; i < 2; i++)                 ELEMENT->Nfve[i] = -1;
	for (i = 0; i < 8; i++)                 ELEMENT->VeC[i]  = -1;
	for (i = 0, iMax = 12*2; i < iMax; i++) ELEMENT->VeE[i]  = -1;
	for (i = 0, iMax = 6*4;  i < iMax; i++) ELEMENT->VeF[i]  = -1;

	// Operators

	// Normals
	ELEMENT->nr = malloc(6*3 * sizeof *(ELEMENT->nr)); // free

    for (i = 0, iMax = 6*3; i < iMax; i++) ELEMENT->nr[i] = 0.;

	// VOLUME Nodes
	ELEMENT->xir_vGs  = malloc(1  * sizeof *(ELEMENT->xir_vGs));  // free
	ELEMENT->xir_vGc  = malloc(NP * sizeof *(ELEMENT->xir_vGc));  // free
	ELEMENT->xir_vCs  = malloc(NP * sizeof *(ELEMENT->xir_vCs));  // free
	ELEMENT->xir_vCc  = malloc(NP * sizeof *(ELEMENT->xir_vCc));  // free
	ELEMENT->xir_vJs  = malloc(NP * sizeof *(ELEMENT->xir_vJs));  // free
	ELEMENT->xir_vJc  = malloc(NP * sizeof *(ELEMENT->xir_vJc));  // free
	ELEMENT->xir_vS   = malloc(NP * sizeof *(ELEMENT->xir_vS));   // free
	ELEMENT->xir_vF   = malloc(NP * sizeof *(ELEMENT->xir_vF));   // free
	ELEMENT->xir_vFrs = malloc(NP * sizeof *(ELEMENT->xir_vFrs)); // free
	ELEMENT->xir_vFrc = malloc(NP * sizeof *(ELEMENT->xir_vFrc)); // free
	ELEMENT->xir_vIs  = malloc(NP * sizeof *(ELEMENT->xir_vIs));  // free
	ELEMENT->xir_vIc  = malloc(NP * sizeof *(ELEMENT->xir_vIc));  // free
	ELEMENT->xir_vP   = malloc(1  * sizeof *(ELEMENT->xir_vP));   // free

	ELEMENT->WvIs = malloc(NP * sizeof *(ELEMENT->WvIs)); // free
	ELEMENT->WvIc = malloc(NP * sizeof *(ELEMENT->WvIc)); // free

	ELEMENT->Con_xir_vP = malloc(1 * sizeof *(ELEMENT->Con_xir_vP)); // free

	ELEMENT->xir_vGs[0] = NULL;
	ELEMENT->xir_vP[0]  = NULL;

	ELEMENT->Con_xir_vP[0] = NULL;
	for (P = 0; P < NP; P++) {
		ELEMENT->xir_vGc[P]  = NULL;
		ELEMENT->xir_vCs[P]  = NULL;
		ELEMENT->xir_vCc[P]  = NULL;
		ELEMENT->xir_vJs[P]  = NULL;
		ELEMENT->xir_vJc[P]  = NULL;
		ELEMENT->xir_vS[P]   = NULL;
		ELEMENT->xir_vF[P]   = NULL;
		ELEMENT->xir_vFrs[P] = NULL;
		ELEMENT->xir_vFrc[P] = NULL;
		ELEMENT->xir_vIs[P]  = NULL;
		ELEMENT->xir_vIc[P]  = NULL;

		ELEMENT->WvIs[P] = NULL;
		ELEMENT->WvIc[P] = NULL;
	}

	ELEMENT->NvnGs  = malloc(1  * sizeof *(ELEMENT->NvnGs));  // free
	ELEMENT->NvnGc  = malloc(NP * sizeof *(ELEMENT->NvnGc));  // free
	ELEMENT->NvnCs  = malloc(NP * sizeof *(ELEMENT->NvnCs));  // free
	ELEMENT->NvnCc  = malloc(NP * sizeof *(ELEMENT->NvnCc));  // free
	ELEMENT->NvnJs  = malloc(NP * sizeof *(ELEMENT->NvnJs));  // free
	ELEMENT->NvnJc  = malloc(NP * sizeof *(ELEMENT->NvnJc));  // free
	ELEMENT->NvnS   = malloc(NP * sizeof *(ELEMENT->NvnS));   // free
	ELEMENT->NvnF   = malloc(NP * sizeof *(ELEMENT->NvnF));   // free
	ELEMENT->NvnFrs = malloc(NP * sizeof *(ELEMENT->NvnFrs)); // free
	ELEMENT->NvnFrc = malloc(NP * sizeof *(ELEMENT->NvnFrc)); // free
	ELEMENT->NvnIs  = malloc(NP * sizeof *(ELEMENT->NvnIs));  // free
	ELEMENT->NvnIc  = malloc(NP * sizeof *(ELEMENT->NvnIc));  // free
	ELEMENT->NvnP   = malloc(1  * sizeof *(ELEMENT->NvnP));   // free

	// FACET Nodes
	ELEMENT->xir_fGc  = malloc(NP * sizeof *(ELEMENT->xir_fGc));  // free
	ELEMENT->xir_fIs  = malloc(NP * sizeof *(ELEMENT->xir_fIs));  // free
	ELEMENT->xir_fIc  = malloc(NP * sizeof *(ELEMENT->xir_fIc));  // free

	ELEMENT->WfIs  = malloc(NP * sizeof *(ELEMENT->WfIs));  // free
	ELEMENT->WfIc  = malloc(NP * sizeof *(ELEMENT->WfIc));  // free

	for (P = 0; P < NP; P++) {
		ELEMENT->xir_fGc[P] = NULL;
		ELEMENT->xir_fIs[P] = NULL;
		ELEMENT->xir_fIc[P] = NULL;

		ELEMENT->WfIs[P] = NULL;
		ELEMENT->WfIc[P] = NULL;
	}

	ELEMENT->NfnGc = malloc(NP * sizeof *(ELEMENT->NfnGc)); // free
	ELEMENT->NfnIs = malloc(NP * sizeof *(ELEMENT->NfnIs)); // free
	ELEMENT->NfnIc = malloc(NP * sizeof *(ELEMENT->NfnIc)); // free


	ELEMENT->next = NULL;

	return ELEMENT;
}


struct S_VOLUME *New_VOLUME(void)
{
	int NP = DB.NP;

	struct S_VOLUME *VOLUME;
	VOLUME = malloc(sizeof *VOLUME); // free

	// Structures
	VOLUME->type   = -1;
	VOLUME->Eclass = -1;
	VOLUME->curved = -1;

	// Geometry
	VOLUME->XYZs = malloc(NP * sizeof *(VOLUME->XYZs)); // free

	VOLUME->next = NULL;

	return VOLUME;
}
