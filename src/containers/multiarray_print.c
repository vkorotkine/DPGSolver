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

#include "multiarray_print.h"

#include "macros.h"

#include "multiarray.h"
#include "matrix.h"
#include "vector.h"

// Static function declarations ************************************************************************************* //

/// \brief Print the extents of the Multiarray.
static void print_Multiarray_extents
	(const int order,              ///< Defined in \ref Multiarray_d.
	 const ptrdiff_t*const extents ///< Defined in \ref Multiarray_d.
	);

/// \brief Print the counter for the indices of order > 2 when printing sub-Matrices of the Multiarray.
static void print_Multiarray_counter
	(const int order,              ///< Defined in \ref Multiarray_d.
	 const ptrdiff_t*const counter ///< The counter for the indices of order > 2.
	);

///	\brief Increment the counter for the indicies of order > 2 by 1.
static void increment_counter
	(const int order,               ///< Defined in \ref Multiarray_d.
	 const ptrdiff_t*const extents, ///< Defined in \ref Multiarray_d.
	 ptrdiff_t*const counter        ///< The counter for the indices of order > 2.
	);

// Interface functions ********************************************************************************************** //

void print_Multiarray_Vector_i (const struct Multiarray_Vector_i*const a)
{
	print_Multiarray_extents(a->order,a->extents);

	const ptrdiff_t size = compute_size(a->order,a->extents);

	for (ptrdiff_t i = 0; i < size; i++)
		print_Vector_i(a->data[i]);
	printf("\n");
}

void print_const_Multiarray_Vector_i (const struct const_Multiarray_Vector_i*const a)
{
	print_Multiarray_Vector_i((struct Multiarray_Vector_i*)a);
}

void print_Multiarray_d (const struct Multiarray_d*const a, const double tol)
{
	const int order               = a->order;
	const ptrdiff_t*const extents = a->extents;

	print_Multiarray_extents(order,extents);

	if (order == 1) {
		struct Vector_d* a_V = constructor_move_Vector_d_d(extents[0],false,a->data); // destructed
		print_Vector_d(a_V,tol);
		destructor_Vector_d(a_V);
	} else if (order == 2) {
		struct Matrix_d* a_M = constructor_move_Matrix_d_d(a->layout,extents[0],extents[1],false,a->data); // destructed
		print_Matrix_d(a_M,tol);
		destructor_Matrix_d(a_M);
	} else if (order > 2) {
		struct Matrix_d* a_M = constructor_move_Matrix_d_d(a->layout,extents[0],extents[1],false,NULL); // destructed;

		const ptrdiff_t size_tail = compute_size(order-2,&extents[2]);

		ptrdiff_t counter[size_tail];
		for (ptrdiff_t i = 0; i < size_tail; ++i)
			counter[i] = 0;

		for (ptrdiff_t i = 0; i < size_tail; ++i) {
			print_Multiarray_counter(order,counter);

			const ptrdiff_t ind_M = compute_index_sub_matrix(order,extents,counter);
			a_M->data = &a->data[ind_M];
			print_Matrix_d(a_M,tol);

			increment_counter(order,extents,counter);
		}
		destructor_Matrix_d(a_M);
	} else {
		EXIT_UNSUPPORTED;
	}
	printf("\n");
}

void print_const_Multiarray_d (const struct const_Multiarray_d*const a, const double tol)
{
	print_Multiarray_d((const struct Multiarray_d*)a,tol);
}

void print_Multiarray_Matrix_d (const struct Multiarray_Matrix_d*const a, const double tol)
{
	print_Multiarray_extents(a->order,a->extents);

	const ptrdiff_t size = compute_size(a->order,a->extents);

	for (ptrdiff_t i = 0; i < size; i++)
		print_Matrix_d(a->data[i],tol);
	printf("\n");
}

void print_const_Multiarray_Matrix_d (const struct const_Multiarray_Matrix_d*const a, const double tol)
{
	print_Multiarray_Matrix_d((struct Multiarray_Matrix_d*)a,tol);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void print_Multiarray_extents (const int order, const ptrdiff_t*const extents)
{
	printf("Multi-array extents: {");
	for (ptrdiff_t i = 0; i < order; i++)
		printf(" %td,",extents[i]);
	printf(" }\n\n");
}

static void print_Multiarray_counter (const int order, const ptrdiff_t*const counter)
{
	printf("{:,:");
	for (int i = 0; i < order-2; ++i)
		printf(",%td",counter[i]);
	printf("}\n");
}

static void increment_counter (const int order, const ptrdiff_t*const extents, ptrdiff_t*const counter)
{
	const ptrdiff_t*const extents_tail = &extents[2];

	for (int i = 0; i < order-2; ++i) {
		if (counter[i] == extents_tail[i]) {
			counter[i] = 0;
			continue;
		}
		++counter[i];
		return;
	}
	EXIT_ERROR("Counter greater than extents.");
}
