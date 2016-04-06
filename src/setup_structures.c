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
		 *EToPrt = DB.EToPrt,
		 *VC     = DB.VC,

	     Testing = DB.Testing;

	int  PrintTesting = 0, MPIrank = DB.MPIrank;

	// Standard datatypes
	int i, v,
	    IndE, IndVC,
	    Vs, firstV, vlocal, NVlocal;

	struct S_VOLUME *VOLUME;


//printf("EToPrt:\n"); array_print_i(1,DB.NE[d],DB.EToPrt,'R');


	Vs = 0; for (i = 0; i < d; i++) Vs += NE[i];

	VOLUME = New_VOLUME();
	DB.VOLUME = VOLUME;

	// Note: only initialize volumes on the current processor.
	for (v = 0, NVlocal = 0; v < NV; v++) {
		if (EToPrt[v] == MPIrank)
			NVlocal++;
	}

	for (v = 0, IndE = Vs, IndVC = 0, firstV = 0, vlocal = 0; v < NV; v++) {
		if (EToPrt[v] == MPIrank) {

			if (firstV != 0)
				VOLUME = VOLUME->next;

			VOLUME->indexl = vlocal;
			VOLUME->indexg = v;
			VOLUME->P      = P;
			VOLUME->type   = EType[IndE];
			VOLUME->Eclass = get_element_class(VOLUME->type);

			if (v == VC[IndVC]) {
				VOLUME->curved = 1;
				IndVC++;
			} else {
				VOLUME->curved = 0;
			}

			if (vlocal != NVlocal-1)
				VOLUME->next = New_VOLUME();

			firstV = 1;
			vlocal++;
		} else {
			// Ensure that appropriate global indices are incremented if necessary
			if (v == VC[IndVC])
				IndVC++;
		}

		IndE++;
	}

	if (IndVC > NVC)
		printf("Error: Found too many curved VOLUMEs.\n"), exit(1);

	// Assign/Overwrite DB parameters 
	DB.NV = NVlocal;

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
