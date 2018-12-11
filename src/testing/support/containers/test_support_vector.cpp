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
/**	\file
 */

#include "test_support_vector.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "test_support.h"
#include "test_support_multiarray.h"

#include "macros.h"
#include "definitions_alloc.h"

#include "file_processing.h"
#include "vector.h"
#include "math_functions.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //
// Difference functions ********************************************************************************************* //

bool diff_Vector_i (const struct Vector_i*const a, const struct Vector_i*const b)
{
	const ptrdiff_t size = a->ext_0;

	if (size != b->ext_0)
		return true;

	for (ptrdiff_t i = 0; i < size; ++i) {
		if (a->data[i] != b->data[i])
			return true;
	}
	return false;
}

bool diff_Vector_d (const struct Vector_d*const a, const struct Vector_d*const b, const double tol)
{
	const ptrdiff_t size = a->ext_0;

	if (size != b->ext_0)
		return true;

	for (ptrdiff_t i = 0; i < size; ++i) {
		if (!equal_d(a->data[i],b->data[i],tol))
			return true;
	}
	return false;
}

bool diff_const_Vector_d (const struct const_Vector_d*const a, const struct const_Vector_d*const b, const double tol)
{
	return diff_Vector_d((const struct Vector_d*const)a,(const struct Vector_d*const)b,tol);
}

bool diff_const_Vector_i (const struct const_Vector_i*const a, const struct const_Vector_i*const b)
{
	return diff_Vector_i((const struct Vector_i*const)a,(const struct Vector_i*const)b);
}

// Printing functions *********************************************************************************************** //

void print_diff_Vector_d (const struct Vector_d*const a, const struct Vector_d*const b, const double tol)
{
	const ptrdiff_t size = a->ext_0;

	if (size != b->ext_0) {
		printf("Note: Attempting to compare Vectors of different size:\n");
		printf("a (%td):\n",a->ext_0);
		print_Vector_d_tol(a,tol);
		printf("b (%td):\n",b->ext_0);
		print_Vector_d_tol(b,tol);
		return;
	}

	const double*const data_a = a->data,
	            *const data_b = b->data;

	double data[size];
	for (ptrdiff_t i = 0; i < size; ++i)
		data[i] = norm_diff_d(1,&data_a[i],&data_b[i],"Inf");

	// Temporarily modify the data pointer to print the vector of differences without constructing a new object.
	double* data_ptr = a->data;

	struct Vector_d* a_tmp = (struct Vector_d*) a;
	a_tmp->data = data;
	print_Vector_d_tol(a_tmp,tol);
	a_tmp->data = data_ptr;
}

void print_diff_const_Vector_d
	(const struct const_Vector_d*const a, const struct const_Vector_d*const b, const double tol)
{
	print_diff_Vector_d((const struct Vector_d*)a,(const struct Vector_d*)b,tol);
}

void print_diff_Vector_i (const struct Vector_i*const a, const struct Vector_i*const b)
{
	const ptrdiff_t size = a->ext_0;

	if (size != b->ext_0) {
		printf("Note: Attempting to compare Vectors of different size:\n");
		print_Vector_i(a);
		print_Vector_i(b);
		return;
	}

	const int*const data_a = a->data,
	         *const data_b = b->data;

	for (ptrdiff_t i = 0; i < size; i++) {
		printf("% 12d ",data_a[i]-data_b[i]);
		if (!((i+1)%8))
			printf("\n");
	}
	printf("\n\n");
}

void print_diff_const_Vector_i
	(const struct const_Vector_i*const a, const struct const_Vector_i*const b)
{
	print_diff_Vector_i((const struct Vector_i*)a,(const struct Vector_i*)b);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

