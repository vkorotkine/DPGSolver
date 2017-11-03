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

#include "test_support_complex_multiarray.h"

#include <assert.h>

#include "macros.h"

#include "multiarray.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void set_Multiarray_c_Multiarray_d (struct Multiarray_c* a, const struct const_Multiarray_d* b)
{
	const ptrdiff_t size = compute_size(a->order,a->extents);
	assert(size == compute_size(b->order,b->extents));

	for (int i = 0; i < size; ++i)
		a->data[i] = b->data[i];
}

double complex* get_col_Multiarray_c (const ptrdiff_t col, struct Multiarray_c* a)
{
	assert(a->layout == 'C');

	const ptrdiff_t ext_0 = a->extents[0];
	return &a->data[col*ext_0];
}

const double complex* get_col_const_Multiarray_c (const ptrdiff_t col, const struct const_Multiarray_c* a)
{
	return (const double complex*) get_col_Multiarray_c(col,(struct Multiarray_c*)a);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
