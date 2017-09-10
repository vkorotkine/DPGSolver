// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/** \file
 */

#include "matrix_math.h"

#include <string.h>
#include <math.h>
#include "mkl.h"

#include "macros.h"
#include "definitions_mkl.h"

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

void mm_d
	(const char trans_a_i, const char trans_b_i, const double alpha, const double beta,
	 const struct const_Matrix_d*const a, const struct const_Matrix_d*const b, struct Matrix_d*const c)
{
	const CBLAS_LAYOUT    layout = ( c->layout == 'R' ? CBRM : CBCM );
	const CBLAS_TRANSPOSE transa = ( (c->layout == a->layout) == (trans_a_i == 'N') ? CBNT : CBT ),
	                      transb = ( (c->layout == b->layout) == (trans_b_i == 'N') ? CBNT : CBT );
	const MKL_INT m = c->ext_0,
	              n = c->ext_1,
	              k = ( trans_a_i == 'N' ? a->ext_1 : a->ext_0 );

	if ((m <= 0) || (n <= 0) || (k <= 0)) {
		EXIT_ERROR("Negative dimension: %d %d %d\n",m,n,k);
	} else if ((m != ( trans_a_i == 'N' ? a->ext_0 : a->ext_1 )) ||
	           (n != ( trans_b_i == 'N' ? b->ext_1 : b->ext_0 )) ||
	           (k != ( trans_b_i == 'N' ? b->ext_0 : b->ext_1 )))
	{
		printf("Invalid matrix dimensions: (m = (%d,%td), n = (%d,%td), k = (%d,%td)\n",
		       m,( trans_a_i == 'N' ? a->ext_0 : a->ext_1 ),
		       n,( trans_b_i == 'N' ? b->ext_1 : b->ext_0 ),
		       k,( trans_b_i == 'N' ? b->ext_0 : b->ext_1 ));
		EXIT_UNSUPPORTED;
	}

	const MKL_INT lda = ( a->layout == 'R' ? a->ext_1 : a->ext_0 ),
	              ldb = ( b->layout == 'R' ? b->ext_1 : b->ext_0 ),
	              ldc = ( c->layout == 'R' ? c->ext_1 : c->ext_0 );

	cblas_dgemm(layout,transa,transb,m,n,k,alpha,a->data,lda,b->data,ldb,beta,c->data,ldc);
}

void mm_NN_d
	(const double alpha, const double beta,
	 const struct const_Matrix_d*const a, const struct const_Matrix_d*const b, struct Matrix_d*const c)
{
	mm_d('N','N',alpha,beta,a,b,c);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void swap_layout (struct Matrix_d* a)
{
	a->layout = ( a->layout == 'R' ? 'C' : 'R' );
}
