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

#include "complex_matrix_constructors.h"

#include <stdlib.h>
#include <assert.h>
#include "definitions_mkl.h"
#include "mkl.h"

#include "complex_matrix.h"
#include "matrix.h"

// Static function declarations ************************************************************************************* //

/** \brief Copy constructor for a \ref Matrix_c\* from a \ref Matrix_d\*.
 *  \return See brief. */
static struct Matrix_c* constructor_copy_Matrix_c_Matrix_d
	(struct Matrix_d* src ///< The source matrix.
	);

// Interface functions ********************************************************************************************** //
// Default constructors ********************************************************************************************* //

// Empty constructors *********************************************************************************************** //

struct Matrix_c* constructor_empty_Matrix_c (const char layout, const ptrdiff_t ext_0, const ptrdiff_t ext_1)
{
	double complex* data = malloc(ext_0*ext_1 * sizeof *data); // keep
	return constructor_move_Matrix_c_c(layout,ext_0,ext_1,true,data);
}

// Copy constructors ************************************************************************************************ //

struct Matrix_c* constructor_copy_Matrix_c (const struct Matrix_c* src)
{
	const ptrdiff_t size = (src->ext_0)*(src->ext_1);
	const double complex*const data_src = src->data;

	double complex* data = malloc(size * sizeof *data); // keep
	for (ptrdiff_t i = 0; i < size; i++)
		data[i] = data_src[i];

	return constructor_move_Matrix_c_c(src->layout,src->ext_0,src->ext_1,true,data);
}

const struct const_Matrix_c* constructor_copy_const_Matrix_c_Matrix_d (const struct const_Matrix_d* src)
{
	return (struct const_Matrix_c*) constructor_copy_Matrix_c_Matrix_d((struct Matrix_d*)src);
}

// Move constructors ************************************************************************************************ //

struct Matrix_c* constructor_move_Matrix_c_c
	(const char layout, const ptrdiff_t ext_0, const ptrdiff_t ext_1, const bool owns_data, double complex*const data)
{
    struct Matrix_c* dest = calloc(1,sizeof *dest); // returned

    dest->layout    = layout;
    dest->ext_0     = ext_0;
    dest->ext_1     = ext_1;
    dest->owns_data = owns_data;
    dest->data      = data;

    return dest;
}

const struct const_Matrix_c* constructor_move_const_Matrix_c_c
	(const char layout, const ptrdiff_t ext_0, const ptrdiff_t ext_1, const bool owns_data,
	 const double complex*const data)
{
	return (const struct const_Matrix_c*)
		constructor_move_Matrix_c_c(layout,ext_0,ext_1,owns_data,(double complex*)data);
}

// Special constructors ********************************************************************************************* //

struct Matrix_c* constructor_sysv_Matrix_c (struct Matrix_c* A_i, struct Matrix_c* B_i)
{
	assert(A_i->layout == B_i->layout); // Can be made flexible in future if necessary.
	assert(A_i->ext_0 == A_i->ext_1);

	// The source matrix is copied as the entries would otherwise be modified while solving the linear system.
/// \todo check if only half of the A matrix needs to be copied.
	struct Matrix_c* A = constructor_copy_Matrix_c(A_i); // destructed
	struct Matrix_c* X = constructor_copy_Matrix_c(B_i); // returned

	const int matrix_layout = ( A->layout == 'R' ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR );
	const lapack_int n      = A->ext_0,
	                 nrhs   = B_i->ext_1;
	double complex* a       = A->data,
	              * x       = X->data;
	const lapack_int lda    = A->ext_0,
	                 ldx    = ( matrix_layout == LAPACK_COL_MAJOR ? n : nrhs );
	lapack_int ipiv[n];

	const lapack_int info = LAPACKE_zsysv(matrix_layout,'U',n,nrhs,a,lda,ipiv,x,ldx);
	assert(info == 0);

	destructor_Matrix_c(A);
	return X;
}

const struct const_Matrix_c* constructor_sysv_const_Matrix_c
	(const struct const_Matrix_c* A_i, const struct const_Matrix_c* B_i)
{
	return (const struct const_Matrix_c*) constructor_sysv_Matrix_c((struct Matrix_c*)A_i,(struct Matrix_c*)B_i);
}

struct Matrix_c* constructor_mm_Matrix_cc
	(const char trans_a_i, const char trans_b_i, const double alpha,
	 const struct const_Matrix_c*const a, const struct const_Matrix_c*const b, const char layout)
{
	const MKL_INT m = ( trans_a_i == 'N' ? a->ext_0 : a->ext_1 ),
	              n = ( trans_b_i == 'N' ? b->ext_1 : b->ext_0 );

	struct Matrix_c* c = constructor_empty_Matrix_c(layout,m,n); // returned

	mm_ccc(trans_a_i,trans_b_i,alpha,0.0,a,b,c);

	return c;
}

const struct const_Matrix_c* constructor_mm_const_Matrix_cc
	(const char trans_a_i, const char trans_b_i, const double alpha,
	 const struct const_Matrix_c*const a, const struct const_Matrix_c*const b, const char layout)
{
	return (const struct const_Matrix_c*) constructor_mm_Matrix_cc(trans_a_i,trans_b_i,alpha,a,b,layout);
}

struct Matrix_c* constructor_mm_diag_Matrix_c_d
	(const double alpha, const struct const_Matrix_c*const a, const struct const_Vector_d*const b, const char side,
	 const bool invert_diag)
{
	struct Matrix_c* c = constructor_copy_Matrix_c((struct Matrix_c*)a);
	scale_Matrix_c_by_Vector_d(side,alpha,c,b,invert_diag);
	return c;
}

const struct const_Matrix_c* constructor_mm_diag_const_Matrix_c_d
	(const double alpha, const struct const_Matrix_c*const a, const struct const_Vector_d*const b, const char side,
	 const bool invert_diag)
{
	return (const struct const_Matrix_c*) constructor_mm_diag_Matrix_c_d(alpha,a,b,side,invert_diag);
}

// Destructors ****************************************************************************************************** //

void destructor_Matrix_c (struct Matrix_c* a)
{
	assert(a != NULL);

	if (a->owns_data)
		free(a->data);
	free(a);
}

void destructor_const_Matrix_c (const struct const_Matrix_c* a)
{
	destructor_Matrix_c((struct Matrix_c*)a);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct Matrix_c* constructor_copy_Matrix_c_Matrix_d (struct Matrix_d* src)
{
	const ptrdiff_t size = (src->ext_0)*(src->ext_1);
	const double*const data_src = src->data;

	double complex* data = calloc(size , sizeof *data); // keep
//	double complex* data = malloc(size * sizeof *data); // keep
	for (ptrdiff_t i = 0; i < size; i++)
		data[i] = data_src[i];

	return constructor_move_Matrix_c_c(src->layout,src->ext_0,src->ext_1,true,data);
}
