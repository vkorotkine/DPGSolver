/* {{{
This file is part of DPGSolver.

DPGSolver is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or any later version.

DPGSolver is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along with DPGSolver.  If not, see
<http://www.gnu.org/licenses/>.
}}} */
/** \file
 */

#include "complex_matrix_math.h"

#include <assert.h>
#include <stddef.h>
#include <complex.h>
#include "gsl/gsl_permute_matrix_complex_double.h"

#include "macros.h"
#include "definitions_mkl.h"
#include "mkl.h"

#include "complex_matrix.h"
#include "complex_vector.h"
#include "matrix.h"
#include "vector.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void transpose_Matrix_c (struct Matrix_c* a, const bool mem_only)
{
	if (a->layout == 'R')
		mkl_zimatcopy(a->layout,'T',a->ext_0,a->ext_1,1.0,a->data,a->ext_1,a->ext_0);
	else if (a->layout == 'C')
		mkl_zimatcopy(a->layout,'T',a->ext_0,a->ext_1,1.0,a->data,a->ext_0,a->ext_1);

	if (mem_only) {
		swap_layout(&a->layout);
	} else {
		ptrdiff_t tmp = a->ext_0;
		a->ext_0 = a->ext_1;
		a->ext_1 = tmp;
	}
}

void permute_Matrix_c (struct Matrix_c* a, const ptrdiff_t* p)
{
	assert((a->layout == 'R') || a->layout == 'C');
	assert(p != NULL);

	if (a->layout == 'R') {
		const int n_p = a->ext_0;

		size_t perm_st[n_p];
		for (int i = 0; i < n_p; ++i)
			perm_st[i] = p[i];

		const gsl_permutation perm = { .size = n_p, .data = perm_st, };

		transpose_Matrix_c(a,false);
		gsl_matrix_complex A =
			{ .size1 = a->ext_0,
			  .size2 = a->ext_1,
			  .tda   = a->ext_1,
			  .data  = (double*) a->data,
			  .block = NULL,
			  .owner = 0, };

		gsl_permute_matrix_complex(&perm,&A);
		transpose_Matrix_c(a,false);
	} else {
		const int n_p = a->ext_1;

		size_t perm_st[n_p];
		for (int i = 0; i < n_p; ++i)
			perm_st[i] = p[i];

		const gsl_permutation perm = { .size = n_p, .data = perm_st, };

		transpose_Matrix_c(a,true);
		gsl_matrix_complex A =
			{ .size1 = a->ext_0,
			  .size2 = a->ext_1,
			  .tda   = a->ext_1,
			  .data  = (double*) a->data,
			  .block = NULL,
			  .owner = 0, };

		gsl_permute_matrix_complex(&perm,&A);
		transpose_Matrix_c(a,true);
	}
}

void scale_Matrix_c_by_Vector_d
	(const char side, const double alpha, struct Matrix_c*const a, const struct const_Vector_d*const b,
	 const bool invert_diag)
{
	assert(invert_diag == false); // Can be made flexible if necessary.
	assert(alpha == 1.0);

	const ptrdiff_t n_row = a->ext_0,
	                n_col = a->ext_1;

	bool transpose_a = false;
	if (side == 'L') {
		assert(b->ext_0 == a->ext_0);

		if (a->layout == 'C') {
			transpose_a = true;
			transpose_Matrix_c(a,true);
		}

		for (ptrdiff_t row = 0; row < n_row; ++row) {
			const double val = b->data[row];
			double complex* data_row = get_row_Matrix_c(row,a);
			for (ptrdiff_t col = 0; col < n_col; ++col)
				*data_row++ *= val;
		}
	} else if (side == 'R') {
		assert(b->ext_0 == a->ext_1);

		if (a->layout == 'R') {
			transpose_a = true;
			transpose_Matrix_c(a,true);
		}

		for (ptrdiff_t col = 0; col < n_col; ++col) {
			const double val = b->data[col];
			double complex* data_col = get_col_Matrix_c(col,a);
			for (ptrdiff_t row = 0; row < n_row; ++row)
				*data_col++ *= val;
		}
	} else {
		EXIT_UNSUPPORTED;
	}
	if (transpose_a)
		transpose_Matrix_c(a,true);
}

void mm_dcc
	(const char trans_a_i, const char trans_b_i, const double alpha, const double beta,
	 const struct const_Matrix_d*const a, const struct const_Matrix_c*const b, struct Matrix_c*const c)
{
	const struct const_Matrix_c*const a_c = constructor_copy_const_Matrix_c_Matrix_d(a); // destructed
	mm_ccc(trans_a_i,trans_b_i,alpha,beta,a_c,b,c);
	destructor_const_Matrix_c(a_c);
}

void mm_cdc
	(const char trans_a_i, const char trans_b_i, const double alpha, const double beta,
	 const struct const_Matrix_c*const a, const struct const_Matrix_d*const b, struct Matrix_c*const c)
{
	const struct const_Matrix_c*const b_c = constructor_copy_const_Matrix_c_Matrix_d(b); // destructed
	mm_ccc(trans_a_i,trans_b_i,alpha,beta,a,b_c,c);
	destructor_const_Matrix_c(b_c);
}

void mm_ccc
	(const char trans_a_i, const char trans_b_i, const double alpha, const double beta,
	 const struct const_Matrix_c*const a, const struct const_Matrix_c*const b, struct Matrix_c*const c)
{
	const CBLAS_LAYOUT    layout = ( c->layout == 'R' ? CBRM : CBCM );
	const CBLAS_TRANSPOSE transa = ( (c->layout == a->layout) == (trans_a_i == 'N') ? CBNT : CBT ),
	                      transb = ( (c->layout == b->layout) == (trans_b_i == 'N') ? CBNT : CBT );
	const MKL_INT m   = c->ext_0,
	              n   = c->ext_1,
	              k   = ( trans_a_i == 'N' ? a->ext_1 : a->ext_0 ),
	              lda = ( a->layout == 'R' ? a->ext_1 : a->ext_0 ),
	              ldb = ( b->layout == 'R' ? b->ext_1 : b->ext_0 ),
	              ldc = ( c->layout == 'R' ? c->ext_1 : c->ext_0 );

	assert(m == ( trans_a_i == 'N' ? a->ext_0 : a->ext_1 ));
	assert(n == ( trans_b_i == 'N' ? b->ext_1 : b->ext_0 ));
	assert(k == ( trans_b_i == 'N' ? b->ext_0 : b->ext_1 ));

	const double complex alpha_c = alpha,
	                     beta_c  = beta;

	cblas_zgemm(layout,transa,transb,m,n,k,&alpha_c,a->data,lda,b->data,ldb,&beta_c,c->data,ldc);
}

void mm_diag_c
	(const char side, const double alpha, const double beta, const struct const_Matrix_d*const a,
	 const struct const_Vector_c*const b, struct Matrix_c* c, const bool invert_diag)
{
	assert(invert_diag == false); // can be made flexible.
	assert(beta == 1.0);

	assert(a->ext_0 == c->ext_0);
	assert(a->ext_1 == c->ext_1);
	assert(a->layout == c->layout);

	const ptrdiff_t n_row = a->ext_0,
	                n_col = a->ext_1;

	if (side == 'L') {
		assert(b->ext_0 == a->ext_0);

		if (a->layout == 'R') {
			for (ptrdiff_t row = 0; row < n_row; ++row) {
				const double complex val = b->data[row];
				const double* data_a   = get_row_const_Matrix_d(row,a);
				double complex* data_c = get_row_Matrix_c(row,c);
				for (ptrdiff_t col = 0; col < n_col; ++col)
					data_c[col] += alpha*data_a[col]*val;
			}
		} else if (a->layout == 'C') {
			for (ptrdiff_t col = 0; col < n_col; ++col) {
				const double* data_a   = get_col_const_Matrix_d(col,a);
				double complex* data_c = get_col_Matrix_c(col,c);
				for (ptrdiff_t row = 0; row < n_row; ++row) {
					const double complex val = b->data[row];
					data_c[col] += alpha*data_a[col]*val;
				}
			}
		}
	} else if (side == 'R') {
		assert(b->ext_0 == a->ext_1);

		if (a->layout == 'R') {
			for (ptrdiff_t row = 0; row < n_row; ++row) {
				const double* data_a   = get_row_const_Matrix_d(row,a);
				double complex* data_c = get_row_Matrix_c(row,c);
				for (ptrdiff_t col = 0; col < n_col; ++col) {
					const double complex val = b->data[col];
					data_c[col] += alpha*data_a[col]*val;
				}
			}
		} else if (a->layout == 'C') {
			for (ptrdiff_t col = 0; col < n_col; ++col) {
				const double complex val = b->data[col];
				const double* data_a   = get_col_const_Matrix_d(col,a);
				double complex* data_c = get_col_Matrix_c(col,c);
				for (ptrdiff_t row = 0; row < n_row; ++row)
					data_c[row] += alpha*data_a[row]*val;
			}
		}
	} else {
		EXIT_UNSUPPORTED;
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
