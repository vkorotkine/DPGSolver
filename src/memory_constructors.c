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
	ELEMENT->type    = 0;
	ELEMENT->d       = 0;
	ELEMENT->Nve     = 0;
	ELEMENT->Nf      = 0;

	ELEMENT->Nfve = malloc(2    * sizeof *(ELEMENT->Nfve)); // free
	ELEMENT->VeC  = malloc(8    * sizeof *(ELEMENT->VeC));  // free
	ELEMENT->VeE  = malloc(12*2 * sizeof *(ELEMENT->VeE));  // free
	ELEMENT->VeF  = malloc(6*4  * sizeof *(ELEMENT->VeF));  // free

	for (i = 0; i < 2; i++)                 ELEMENT->Nfve[i] = 0;
	for (i = 0; i < 8; i++)                 ELEMENT->VeC[i]  = 0;
	for (i = 0, iMax = 12*2; i < iMax; i++) ELEMENT->VeE[i]  = 0;
	for (i = 0, iMax = 6*4;  i < iMax; i++) ELEMENT->VeF[i]  = 0;

	// Operators

	// Normals
	ELEMENT->nr = malloc(6*3 * sizeof *(ELEMENT->nr)); // free

	for (i = 0, iMax = 6*3; i < iMax; i++) ELEMENT->nr[i] = 0.;

	// VOLUME Nodes
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
	ELEMENT->rst_vP   = malloc(1  * sizeof *(ELEMENT->rst_vP));   // free

	ELEMENT->wvIs = malloc(NP * sizeof *(ELEMENT->wvIs)); // free
	ELEMENT->wvIc = malloc(NP * sizeof *(ELEMENT->wvIc)); // free

	ELEMENT->Con_rst_vP = malloc(1 * sizeof *(ELEMENT->Con_rst_vP)); // free

	ELEMENT->rst_vGs[0] = NULL;
	ELEMENT->rst_vP[0]  = NULL;

	ELEMENT->Con_rst_vP[0] = NULL;
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

	ELEMENT->NvnGs  = malloc(1  * sizeof *(ELEMENT->NvnGs));   // free
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
	ELEMENT->rst_fGc  = malloc(NP * sizeof *(ELEMENT->rst_fGc));  // free
	ELEMENT->rst_fIs  = malloc(NP * sizeof *(ELEMENT->rst_fIs));  // free
	ELEMENT->rst_fIc  = malloc(NP * sizeof *(ELEMENT->rst_fIc));  // free

	ELEMENT->wfIs  = malloc(NP * sizeof *(ELEMENT->wfIs));  // free
	ELEMENT->wfIc  = malloc(NP * sizeof *(ELEMENT->wfIc));  // free

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

	// Operators
	ELEMENT->I_vGs_vGc = malloc(NP * sizeof *(ELEMENT->I_vGs_vGc)); // free




	ELEMENT->next = NULL;

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
	// *XYZc;

	// Geometry
	// *XYZs;

	VOLUME->next    = NULL;
	VOLUME->grpnext = NULL;

	return VOLUME;
}
