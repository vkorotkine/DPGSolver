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
#include "vector.h"

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
	if (a->layout == 'R')
		mkl_dimatcopy(a->layout,'T',a->ext_0,a->ext_1,1.0,a->data,a->ext_1,a->ext_0);
	else if (a->layout == 'C')
		mkl_dimatcopy(a->layout,'T',a->ext_0,a->ext_1,1.0,a->data,a->ext_0,a->ext_1);

	if (mem_only) {
		swap_layout(a);
	} else {
		ptrdiff_t tmp = a->ext_0;
		a->ext_0 = a->ext_1;
		a->ext_1 = tmp;
	}
}

void transpose_const_Matrix_d (const struct const_Matrix_d* a, const bool mem_only)
{
	transpose_Matrix_d((struct Matrix_d*)a,mem_only);
}

void scale_Matrix_d (struct Matrix_d* a, const double val)
{
	ptrdiff_t size = (a->ext_0)*(a->ext_1);
	for (ptrdiff_t i = 0; i < size; ++i)
		a->data[i] *= val;
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
		printf("Invalid matrix dimensions: (m = (%d,%td), n = (%d,%td), k = (%d,%td))\n",
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

void mv_d
	(const char layout_i, const char trans_a_i, const double alpha, const double beta,
	 const struct const_Matrix_d*const a, const struct const_Vector_d*const b, struct Vector_d*const c)
{
	const CBLAS_LAYOUT    layout = ( layout_i == 'R' ? CBRM : CBCM );
	const CBLAS_TRANSPOSE transa = ( (layout_i == a->layout) == (trans_a_i == 'N') ? CBNT : CBT );

	/// \note Unlike the \ref mm_d function, m and n here represent the dimensions of A and not op(A).
	const MKL_INT m = a->ext_0,
	              n = a->ext_1;

	if ((m <= 0) || (n <= 0)) {
		EXIT_ERROR("Negative dimension: %d %d\n",m,n);
	} else if ((m != ( transa == CBNT ? c->ext_0 : b->ext_0 )) ||
	           (n != ( transa == CBNT ? b->ext_0 : c->ext_0 )))
	{
		printf("Invalid matrix/vector dimensions: (m = (%d,%td), n = (%d,%td))\n",
		       m,( transa == CBNT ? b->ext_0 : c->ext_0 ),
		       n,( transa == CBNT ? c->ext_0 : b->ext_0 ));
		EXIT_UNSUPPORTED;
	}

	const MKL_INT lda = ( a->layout == 'R' ? a->ext_1 : a->ext_0 ),
	              ldb = 1,
	              ldc = 1;

	cblas_dgemv(layout,transa,m,n,alpha,a->data,lda,b->data,ldb,beta,c->data,ldc);
}

void scale_Matrix_by_Vector_d
	(const char side, const double alpha, struct Matrix_d*const a, const struct const_Vector_d*const b)
{
	if (alpha != 1.0)
		scale_Matrix_d(a,alpha);

	const ptrdiff_t n_row = a->ext_0,
	                n_col = a->ext_1;

	bool transpose_a = false;
	if (side == 'L') {
		if (b->ext_0 != a->ext_0)
			EXIT_ERROR("Invalid dimensions.");

		if (a->layout == 'C') {
			transpose_a = true;
			transpose_Matrix_d(a,true);
		}

		for (ptrdiff_t row = 0; row < n_row; ++row) {
			const double val = b->data[row];
			double* data_row = get_row_Matrix_d(row,a);
			for (ptrdiff_t col = 0; col < n_col; ++col)
				*data_row++ *= val;
		}
	} else if (side == 'R') {
		if (b->ext_0 != a->ext_1)
			EXIT_ERROR("Invalid dimensions.");

		if (a->layout == 'R') {
			transpose_a = true;
			transpose_Matrix_d(a,true);
		}

		for (ptrdiff_t col = 0; col < n_col; ++col) {
			const double val = b->data[col];
			double* data_col = get_col_Matrix_d(col,a);
			for (ptrdiff_t row = 0; row < n_row; ++row)
				*data_col++ *= val;
		}
	} else {
		EXIT_UNSUPPORTED;
	}
	if (transpose_a)
		transpose_Matrix_d(a,true);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void swap_layout (struct Matrix_d* a)
{
	a->layout = ( a->layout == 'R' ? 'C' : 'R' );
}
