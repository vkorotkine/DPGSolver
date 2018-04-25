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

#include "matrix.h"

#include "definitions_core.h"

#include "vector.h"

// Templated functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "def_templates_matrix.h"
#include "matrix_T.c"
#include "undef_templates_type.h"
#include "undef_templates_matrix.h"

#include "def_templates_type_i.h"
#include "def_templates_matrix_i.h"
#include "matrix_T.c"
#include "undef_templates_type.h"
#include "undef_templates_matrix.h"

// Static function declarations ************************************************************************************* //

/// \brief Container for a DIM-dimensional node with an index.
struct Node {
	double xyz[DIM]; ///< The coordinates.
	int index;       ///< The index of the node.
};

/** \brief Comparison function for std::qsort between `double*` `a` and `b` of length \ref DIM.
 *  \return Standard required return for comparator function of qsort. */
static int cmp_DIM_d
	(const double*const a, ///< Variable 1.
	 const double*const b  ///< Variable 2.
	);

// Interface functions ********************************************************************************************** //

void swap_layout (char*const layout)
{
	*layout = ( *layout == 'R' ? 'C' : 'R' );
}

void swap_layout_and_extents (char*const layout, ptrdiff_t*const ext_0, ptrdiff_t*const ext_1)
{
	swap_layout(layout);

	const ptrdiff_t tmp = *ext_0;
	*ext_0 = *ext_1;
	*ext_1 = tmp;
}

char compute_opposite_layout (const char layout_i)
{
	assert((layout_i == 'R') || (layout_i == 'C'));
	return ( layout_i == 'R' ? 'C' : 'R' );
}

ptrdiff_t compute_index_Matrix
	(const ptrdiff_t i, const ptrdiff_t j, const ptrdiff_t ext_0, const ptrdiff_t ext_1, const char layout)
{
	assert((layout == 'R') || (layout == 'C'));
	return ( layout == 'R' ?  i*ext_1+j : j*ext_0+i );
}

const struct const_Vector_i* row_sort_DIM_Matrix_d (struct Matrix_d*const src, const bool return_indices)
{
	assert(src->layout == 'R');
	assert(src->ext_1 == DIM);

	const ptrdiff_t size = src->ext_0;
	if (!return_indices) {
		qsort(src->data,(size_t)size,DIM*sizeof(src->data[0]),(int (*)(const void *, const void *))cmp_DIM_d);
		return NULL;
	} else {
		struct Node node_data[size];
		for (int i = 0; i < size; ++i) {
			double*const data_src = get_row_Matrix_d(i,src);
			node_data[i].index = i;
			for (int d = 0; d < DIM; ++d)
				node_data[i].xyz[d] = data_src[d];
		}
		qsort(node_data,(size_t)size,sizeof(struct Node),(int (*)(const void *, const void *))cmp_DIM_d);

		int*const ordering_i = malloc((size_t)size * sizeof *ordering_i); // keep
		for (int i = 0; i < size; ++i) {
//			ordering_i[i] = node_data[i].index;
			ordering_i[node_data[i].index] = i;
			double*const data_src = get_row_Matrix_d(i,src);
			for (int d = 0; d < DIM; ++d)
				data_src[d] = node_data[i].xyz[d];
		}
		return constructor_move_const_Vector_i_i(size,true,ordering_i);
	}
	EXIT_ERROR("Should not have made it here.\n");
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static int cmp_DIM_d (const double*const a, const double*const b)
{
	for (int d = DIM-1; d >= 0; --d) {
//	for (int d = 0; d < DIM; ++d) {
		if (a[d] < b[d])
			return -1;
		else if (a[d] > b[d])
			return 1;
	}
	return 0;
}
