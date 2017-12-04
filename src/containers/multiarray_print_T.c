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

#include <assert.h>

#include "macros.h"
#include "definitions_tol.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void print_Multiarray_Vector_T (const struct Multiarray_Vector_T*const a)
{
	print_Multiarray_extents(a->order,a->extents);

	const ptrdiff_t size = compute_size(a->order,a->extents);

	for (ptrdiff_t i = 0; i < size; i++) {
		printf("\nIndex: % 3td: ",i);
		if (a->data[i]) {
			printf("\n\n");
			print_Vector_T(a->data[i]);
		} else {
			printf("*** NULL ***");
		}
	}
	printf("\n");
}

void print_const_Multiarray_Vector_T (const struct const_Multiarray_Vector_T*const a)
{
	print_Multiarray_Vector_T((struct Multiarray_Vector_T*)a);
}

void print_Multiarray_T_tol (const struct Multiarray_T*const a, const Real tol)
{
	const int order               = a->order;
	const ptrdiff_t*const extents = a->extents;

	print_Multiarray_extents(order,extents);

	if (order == 1) {
		struct Vector_T* a_V = constructor_move_Vector_T_T(extents[0],false,a->data); // destructed
		print_Vector_T_tol(a_V,tol);
		destructor_Vector_T(a_V);
	} else if (order == 2) {
		struct Matrix_T* a_M = constructor_move_Matrix_T_T(a->layout,extents[0],extents[1],false,a->data); // destructed
		print_Matrix_T_tol(a_M,tol);
		destructor_Matrix_T(a_M);
	} else if (order > 2) {
		struct Matrix_T* a_M = constructor_move_Matrix_T_T(a->layout,extents[0],extents[1],false,NULL); // destructed

		const ptrdiff_t size_tail = compute_size(order-2,&extents[2]);

		ptrdiff_t counter[order-2];
		for (ptrdiff_t i = 0; i < order-2; ++i)
			counter[i] = 0;

		for (ptrdiff_t i = 0; i < size_tail; ++i) {
			print_Multiarray_counter(order,counter);

			const ptrdiff_t ind_M = compute_index_sub_matrix(order,extents,counter);
			a_M->data = &a->data[ind_M];
			print_Matrix_T_tol(a_M,tol);

			increment_counter(order,extents,counter);
		}
		destructor_Matrix_T(a_M);
	} else {
		EXIT_ERROR("Unsupported: %d\n",order);
	}
	printf("\n");
}

void print_const_Multiarray_T_tol (const struct const_Multiarray_T*const a, const Real tol)
{
	print_Multiarray_T_tol((const struct Multiarray_T*)a,tol);
}

void print_Multiarray_Matrix_T_tol (const struct Multiarray_Matrix_T*const a, const Real tol)
{
	print_Multiarray_extents(a->order,a->extents);

	const ptrdiff_t size = compute_size(a->order,a->extents);

	for (ptrdiff_t i = 0; i < size; i++) {
		printf("\nIndex: % 3td: ",i);
		if (a->data[i]) {
			printf("\n\n");
			print_Matrix_T_tol(a->data[i],tol);
		} else {
			printf("*** NULL ***");
		}
	}
	printf("\n");
}

void print_const_Multiarray_Matrix_T_tol (const struct const_Multiarray_Matrix_T*const a, const Real tol)
{
	print_Multiarray_Matrix_T_tol((struct Multiarray_Matrix_T*)a,tol);
}

void print_Multiarray_T (const struct Multiarray_T*const a)
{
	print_Multiarray_T_tol(a,EPS);
}

void print_const_Multiarray_T (const struct const_Multiarray_T*const a)
{
	print_Multiarray_T((const struct Multiarray_T*)a);
}

void print_Multiarray_Matrix_T (const struct Multiarray_Matrix_T*const a)
{
	print_Multiarray_Matrix_T_tol(a,EPS);
}

void print_const_Multiarray_Matrix_T (const struct const_Multiarray_Matrix_T*const a)
{
	print_Multiarray_Matrix_T((const struct Multiarray_Matrix_T*)a);
}

#ifndef TYPE_RC
void fprint_const_Multiarray_Vector_T (FILE* file, const int n_tab, const struct const_Multiarray_Vector_T* a)
{
	const ptrdiff_t size = compute_size(a->order,a->extents);

	for (ptrdiff_t i = 0; i < size; ++i)
		fprint_const_Vector_T(file,n_tab,a->data[i]);
}
#endif
// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
