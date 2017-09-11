// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/** \file
 */

#include "multiarray_math.h"

#include "macros.h"

#include "multiarray.h"
#include "matrix.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void transpose_Multiarray_d (struct Multiarray_d* a, const bool mem_only)
{
	if (a->order != 2)
		EXIT_UNSUPPORTED;

	struct Matrix_d* dest = constructor_default_Matrix_d(); // destructed
	set_Matrix_from_Multiarray_d(dest,a,(ptrdiff_t[]){});

	transpose_Matrix_d(dest,mem_only);
	destructor_Matrix_d(dest);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
