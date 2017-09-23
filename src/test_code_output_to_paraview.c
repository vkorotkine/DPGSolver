// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "test_code_output_to_paraview.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "Test.h"

/*
 *	Purpose:
 *		Provide functions associated with output_to_paraview used for testing purposes.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

char *get_fNameOut(char *output_type)
{
	char string[STRLEN_MIN], *fNameOut;

	fNameOut = malloc(STRLEN_MAX * sizeof *fNameOut); // keep

	strcpy(fNameOut,output_type);
	sprintf(string,"%dD_",DB.d); strcat(fNameOut,string);
								 strcat(fNameOut,DB.MeshType);

	// Choose one of the two options below to clean this up (probably DB and not TestDB) (ToBeDeleted)
	if (DB.Adapt == ADAPT_0) {
		EXIT_UNSUPPORTED;
		sprintf(string,"_ML%d",DB.ML); strcat(fNameOut,string);
		sprintf(string,"P%d_",DB.PGlobal); strcat(fNameOut,string);
	} else {
		sprintf(string,"_ML%d",TestDB.ML); strcat(fNameOut,string);
		sprintf(string,"P%d_",TestDB.PGlobal); strcat(fNameOut,string);
	}

	return fNameOut;
}
