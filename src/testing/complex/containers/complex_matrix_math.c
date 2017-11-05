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

#include "macros.h"
#include "definitions_mkl.h"
#include "mkl.h"

#include "complex_matrix.h"
#include "matrix.h"

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

void mm_c
	(const char trans_a_i, const char trans_b_i, const double alpha, const double beta,
	 const struct const_Matrix_d*const a, const struct const_Matrix_c*const b, struct Matrix_c*const c)
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

/// \todo Check whether the conversion to complex is necessary.
	const double complex alpha_c = alpha,
	                     beta_c  = beta;

	const struct const_Matrix_c*const a_c = constructor_copy_const_Matrix_c_Matrix_d(a); // destructed
	cblas_zgemm(layout,transa,transb,m,n,k,&alpha_c,a_c->data,lda,b->data,ldb,&beta_c,c->data,ldc);
	destructor_const_Matrix_c(a_c);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
