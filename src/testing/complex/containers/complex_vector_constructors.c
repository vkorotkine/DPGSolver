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
/// \file

#include "complex_vector_constructors.h"

#include <assert.h>
#include <stdlib.h>
#include "definitions_mkl.h"
#include "mkl.h"

#include "complex_matrix.h"
#include "complex_vector.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //
// Default constructors ********************************************************************************************* //

// Empty constructors *********************************************************************************************** //

struct Vector_c* constructor_empty_Vector_c (const ptrdiff_t ext_0)
{
	double complex* data = malloc(ext_0 * sizeof *data); // keep

	return constructor_move_Vector_c_c(ext_0,true,data);
}

// Zero constructors ************************************************************************************************ //

struct Vector_c* constructor_zero_Vector_c (const ptrdiff_t ext_0)
{
	double complex* data = calloc(ext_0 , sizeof *data); // keep
	return constructor_move_Vector_c_c(ext_0,true,data);
}

// Copy constructors ************************************************************************************************ //

// Move constructors ************************************************************************************************ //

struct Vector_c* constructor_move_Vector_c_c (const ptrdiff_t ext_0, const bool owns_data, double complex*const data)
{
	struct Vector_c* dest = calloc(1,sizeof *dest); // returned

	dest->ext_0     = ext_0;
	dest->owns_data = owns_data;
	dest->data      = data;

	return dest;
}

// Set constructors ************************************************************************************************* //

// Special constructors ********************************************************************************************* //

struct Vector_c* constructor_mv_Vector_c
	(const char trans_a_i, const double alpha, const struct const_Matrix_c*const a,
	 const struct const_Vector_c*const b)
{
	const MKL_INT m = ( trans_a_i == 'N' ? a->ext_0 : a->ext_1 );

	struct Vector_c* c = constructor_empty_Vector_c(m); // returned

	mv_ccc(trans_a_i,alpha,0.0,a,b,c);

	return c;
}

const struct const_Vector_c* constructor_mv_const_Vector_c
	(const char trans_a_i, const double alpha, const struct const_Matrix_c*const a,
	 const struct const_Vector_c*const b)
{
	return (const struct const_Vector_c*) constructor_mv_Vector_c(trans_a_i,alpha,a,b);
}

// Destructors ****************************************************************************************************** //

void destructor_Vector_c (struct Vector_c* a)
{
	assert(a != NULL);

	if (a->owns_data)
		free(a->data);
	free(a);
}

void destructor_const_Vector_c (const struct const_Vector_c* a)
{
	destructor_Vector_c((struct Vector_c*)a);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
