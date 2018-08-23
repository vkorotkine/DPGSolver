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

#include "multiarray.h"

#include <string.h>

#include "definitions_alloc.h"

#include "matrix.h"
#include "vector.h"

#include "math_functions.h"

// Static function declarations ************************************************************************************* //

double CMP_TOL = 0.0; ///< Comparison tolerance used for \ref cmp_Vector_tol_T.

// Interface functions ********************************************************************************************** //

#include "def_templates_type_i.h"
#include "multiarray_T.c"
#include "undef_templates_type.h"

#include "def_templates_type_d.h"
#include "multiarray_T.c"
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "multiarray_T.c"
#include "undef_templates_type.h"

ptrdiff_t compute_size (const int order, const ptrdiff_t*const extents)
{

	ptrdiff_t size = 1;
	for (ptrdiff_t i = 0; i < order; i++)
		size *= extents[i];

	return size;
}

bool check_equal_order_extents
	(const int order_1, const int order_2, const ptrdiff_t*const extents_1, const ptrdiff_t*const extents_2)
{
	if (order_1 != order_2)
		return false;

	for (int i = 0; i < order_1; ++i) {
		if (extents_1[i] != extents_2[i])
			return false;
	}
	return true;
}

ptrdiff_t compute_index_sub_vector (const int order, const ptrdiff_t*const extents, const ptrdiff_t*const sub_indices)
{
	return compute_index_sub_container(order,1,extents,sub_indices);
}

ptrdiff_t compute_index_sub_matrix (const int order, const ptrdiff_t*const extents, const ptrdiff_t*const sub_indices)
{
	return compute_index_sub_container(order,2,extents,sub_indices);
}

ptrdiff_t compute_index_sub_container
	(const int order_i, const int order_o, const ptrdiff_t*const extents, const ptrdiff_t*const sub_indices)
{

	const ptrdiff_t*const extents_tail = &extents[order_o];

	ptrdiff_t base = 1;
	for (int i = 0; i < order_o; ++i)
		base *= extents[i];
	ptrdiff_t ind_sub = 0;
	for (int i = 0; i < order_i-order_o; ++i) {
		ind_sub += base*sub_indices[i];
		base *= extents_tail[i];
	}
	return ind_sub;
}

ptrdiff_t compute_index_sub_container_pi
	(const int order_i, const int order_o, const ptrdiff_t*const extents, const int*const sub_indices)
{
	const int n_sub = order_i-order_o;
	ptrdiff_t sub_indices_p[n_sub];
	for (int i = 0; i < n_sub; ++i)
		sub_indices_p[i] = sub_indices[i];

	return compute_index_sub_container(order_i,order_o,extents,sub_indices_p);
}

void check_container_type (FILE* data_file, const char*const container_type)
{
	char line[STRLEN_MAX];
	if (fgets(line,sizeof(line),data_file) != NULL) {};

	char expected_line[STRLEN_MAX];
	strcpy(expected_line,"container ");
	strcat(expected_line,container_type);

	const bool found = ( strstr(line,expected_line) ? true : false );
	if (!found)
		EXIT_ERROR("Reading incorrect container type: %s. (expected: %s)",line,expected_line);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
