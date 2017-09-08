// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/** \file
 */

#include "matrix_math.h"

#include <string.h>
#include <math.h>
#include "mkl.h"

#include "macros.h"

#include "matrix.h"

// Static function declarations ************************************************************************************* //

/// \brief Swap the layout of the \ref Matrix_d\*.
static void swap_layout
	(struct Matrix_d* a ///< The input matrix.
	);

// Interface functions ********************************************************************************************** //

double compute_norm_Matrix_d_row
	(const ptrdiff_t row, const struct Matrix_d*const a, const char*const norm_type)
{
	const double*const data = get_row_Matrix_d(row,a);
	const ptrdiff_t i_max   = a->ext_1;

	double norm = 0.0;
	if (strstr(norm_type,"L2")) {
		for (ptrdiff_t i = 0; i < i_max; ++i)
			norm += data[i]*data[i];
		return sqrt(norm);
	}
	EXIT_UNSUPPORTED;
}

void transpose_Matrix_d (struct Matrix_d* a, const bool mem_only)
{
	mkl_dimatcopy(a->layout,'T',a->ext_0,a->ext_1,1.0,a->data,a->ext_1,a->ext_0);
	if (mem_only) {
		swap_layout(a);
	} else {
		ptrdiff_t tmp = a->ext_0;
		a->ext_0 = a->ext_1;
		a->ext_1 = tmp;
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void swap_layout (struct Matrix_d* a)
{
	a->layout = ( a->layout == 'R' ? 'C' : 'R' );
}
