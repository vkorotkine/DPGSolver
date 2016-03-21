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
	free(ELEMENT->Nfve);
	free(ELEMENT->VeC);
	free(ELEMENT->VeE);
	free(ELEMENT->VeF);

	free(ELEMENT);
}


