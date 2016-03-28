#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Set up VOLUME and FACET structures.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
*/

int get_element_class(const int type);

void setup_structures()
{
	// Initialize DB Parameters
	int  d       = DB.d,
	     P       = DB.P,
	     NV      = DB.NV,
		 NVC     = DB.NVC,
	     *NE     = DB.NE,
		 *EType  = DB.EType,
		 *VC     = DB.VC,

	     Testing = DB.Testing;

	int  PrintTesting = 0, MPIrank = DB.MPIrank;

	// Standard datatypes
	int i, v,
	    IndE, IndVC,
	    Vs;

	struct S_VOLUME *VOLUME;

	Vs = 0; for (i = 0; i < d; i++) Vs += NE[i];

	VOLUME = New_VOLUME();
	DB.VOLUME = VOLUME;

	for (v = 0, IndE = Vs, IndVC = 0; v < NV; v++) {
		if (v != 0)
			VOLUME = VOLUME->next;

		VOLUME->P      = P;
		VOLUME->type   = EType[IndE];
		VOLUME->Eclass = get_element_class(VOLUME->type);

		if (v == VC[IndVC]) {
			VOLUME->curved = 1;
			IndVC++;
		} else {
			VOLUME->curved = 0;
		}

		if (v != NV-1)
			VOLUME->next = New_VOLUME();

		IndE++;
	}

	if (IndVC > NVC)
		printf("Error: Found too many curved VOLUMEs.\n"), exit(1);

//VOLUME = DB.VOLUME; while(VOLUME != NULL) printf("%d %d\n",VOLUME->type,VOLUME->curved), VOLUME = VOLUME->next;

}

int get_element_class(const int type)
{
	if (type == POINT || type == LINE || type == QUAD || type == HEX)
		return C_TP;
	else if (type == TRI || type == TET)
		return C_SI;
	else if (type == WEDGE)
		return C_WEDGE;
	else if (type == PYR)
		return C_PYR;
	else
		printf("Error: There is not yet an element class associated with the type provided.\n"), exit(1);

}
