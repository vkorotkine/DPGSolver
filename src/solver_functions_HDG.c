// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "solver_functions_HDG.h"

#include <stdlib.h>
#include <stdio.h>

#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
//#include "S_OpCSR.h" // ToBeModified (Add when using sparse matrices)

#include "element_functions.h"
#include "matrix_functions.h"

struct S_OPERATORS_V *init_mat_ops_VOLUME (struct S_VOLUME const *const VOLUME)
{
	struct S_OPERATORS_V *OPS = malloc(sizeof *OPS); // returned

	unsigned int const P = VOLUME->P;

	struct S_ELEMENT const *const ELEMENT = get_ELEMENT_type(VOLUME->type);

	if (!VOLUME->curved) {
		OPS->D_Weak = (struct S_MATRIX const *const *const)
			mat_constructor2_move('R','D',ELEMENT->NvnS[P],ELEMENT->NvnIs[P],DB.d, ELEMENT->Ds_Weak_VV[P][P][0],NULL);
	} else {
		EXIT_UNSUPPORTED;
	}

	return OPS;
}
