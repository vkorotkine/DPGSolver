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

#include "complex_multiarray_print.h"

#include <assert.h>

#include "macros.h"
#include "definitions_tol.h"

#include "complex_matrix.h"
#include "complex_multiarray.h"
#include "multiarray.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void print_Multiarray_c (const struct Multiarray_c*const a)
{
	const int order               = a->order;
	const ptrdiff_t*const extents = a->extents;

	print_Multiarray_extents(order,extents);

	if (order == 1) {
		EXIT_ADD_SUPPORT;
	} else if (order == 2) {
		struct Matrix_c* a_M =
			constructor_move_Matrix_c_c(a->layout,extents[0],extents[1],false,a->data); // destructed
		print_Matrix_c(a_M);
		destructor_Matrix_c(a_M);
	} else if (order > 2) {
		struct Matrix_c* a_M =
			constructor_move_Matrix_c_c(a->layout,extents[0],extents[1],false,NULL); // destructed

		const ptrdiff_t size_tail = compute_size(order-2,&extents[2]);

		ptrdiff_t counter[order-2];
		for (ptrdiff_t i = 0; i < order-2; ++i)
			counter[i] = 0;

		for (ptrdiff_t i = 0; i < size_tail; ++i) {
			print_Multiarray_counter(order,counter);

			const ptrdiff_t ind_M = compute_index_sub_matrix(order,extents,counter);
			a_M->data = &a->data[ind_M];
			print_Matrix_c(a_M);

			increment_counter(order,extents,counter);
		}
		destructor_Matrix_c(a_M);
	} else {
		EXIT_ERROR("Unsupported: %d\n",order);
	}
	printf("\n");
}

void print_const_Multiarray_c (const struct const_Multiarray_c*const a)
{
	print_Multiarray_c((const struct Multiarray_c*)a);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
