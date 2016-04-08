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

int get_element_class(const unsigned int type)
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

void setup_structures(void)
{
	// Initialize DB Parameters
	unsigned int d       = DB.d,
	             P       = DB.P,
	             NP      = DB.NP,
	             NV      = DB.NV,
	             NVC     = DB.NVC,
	             *NE     = DB.NE,
	             *EType  = DB.EType,
	             *EToPrt = DB.EToPrt,
	             *VC     = DB.VC;
	int          MPIrank = DB.MPIrank;

	int  PrintTesting = 0;

	// Standard datatypes
	unsigned int i, iMax, v,
	             IndE, IndVC, IndVgrp,
	             Vs, firstV, vlocal, NVlocal, NECgrp, NVgrp,
	             uMPIrank;

	struct S_VOLUME *VOLUME, **Vgrp, **Vgrp_tmp;

	// Arbitrary initializations for variables defined in conditionals (to eliminate compiler warnings)
	NECgrp = 0;

	uMPIrank = MPIrank;

	if      (d == 1) NECgrp = 1;
	else if (d == 2) NECgrp = 2;
	else if (d == 3) NECgrp = 4;

	NVgrp = NECgrp*NP*2;
	Vgrp     = malloc(NVgrp * sizeof *Vgrp);     // keep
	Vgrp_tmp = malloc(NVgrp * sizeof *Vgrp_tmp); // free
	for (i = 0, iMax = NVgrp; iMax--; i++) {
		Vgrp[i]     = NULL;
		Vgrp_tmp[i] = NULL;
	}

	Vs = 0; for (i = 0; i < d; i++) Vs += NE[i];

	VOLUME = New_VOLUME();
	DB.VOLUME = VOLUME;

	// Note: only initialize volumes on the current processor.
	for (v = 0, NVlocal = 0; v < NV; v++) {
		if (EToPrt[v] == uMPIrank)
			NVlocal++;
	}

	for (v = 0, IndE = Vs, IndVC = 0, firstV = 0, vlocal = 0; v < NV; v++) {
		if (EToPrt[v] == uMPIrank) {

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

			IndVgrp = (VOLUME->Eclass*NP*2)+(VOLUME->P*2)+(VOLUME->curved);
			if (Vgrp[IndVgrp] == NULL)
				Vgrp[IndVgrp] = VOLUME;
			else
				Vgrp_tmp[IndVgrp]->grpnext = VOLUME;

			Vgrp_tmp[IndVgrp] = VOLUME;

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
	free(Vgrp_tmp);

	if (IndVC > NVC)
		printf("Error: Found too many curved VOLUMEs.\n"), exit(1);


for (i = 0, iMax = NVgrp; iMax--; i++) {
	for (VOLUME = Vgrp[i]; VOLUME != NULL; VOLUME = VOLUME->grpnext) {
		printf("%d %d %d %d\n",i,VOLUME->Eclass,VOLUME->P,VOLUME->curved);
	}
	printf("\t\t%p\n",Vgrp[i]);
}
//exit(1);


	// Assign/Overwrite DB parameters
	DB.NV     = NVlocal;
	DB.NECgrp = NECgrp;

	DB.Vgrp = Vgrp;

//VOLUME = DB.VOLUME; while(VOLUME != NULL) printf("%d %d\n",VOLUME->type,VOLUME->curved), VOLUME = VOLUME->next;

}
