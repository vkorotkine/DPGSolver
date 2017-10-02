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

#include "multiarray_math.h"

#include <assert.h>

#include "macros.h"

#include "multiarray.h"
#include "matrix.h"

// Static function declarations ************************************************************************************* //

/// \brief `mutable` version of \ref reinterpret_const_Multiarray_as_Matrix_d.
void reinterpret_Multiarray_as_Matrix_d
	(struct Multiarray_d* a, ///< Defined for \ref reinterpret_const_Multiarray_as_Matrix_d.
	 struct Matrix_d* a_M,   ///< Defined for \ref reinterpret_const_Multiarray_as_Matrix_d.
	 const ptrdiff_t ext_0,  ///< Defined for \ref reinterpret_const_Multiarray_as_Matrix_d.
	 const ptrdiff_t ext_1   ///< Defined for \ref reinterpret_const_Multiarray_as_Matrix_d.
	);

// Interface functions ********************************************************************************************** //

void transpose_Multiarray_d (struct Multiarray_d* a, const bool mem_only)
{
	assert(mem_only == true);
	const int order = a->order;

	const ptrdiff_t ext_0 = a->extents[0],
	                ext_1 = compute_size(order,a->extents)/a->extents[0];

	struct Matrix_d* a_M = constructor_move_Matrix_d_d(a->layout,ext_0,ext_1,false,a->data); // destructed

	transpose_Matrix_d(a_M,mem_only);
	a->layout = a_M->layout;

	destructor_Matrix_d(a_M);
}

void reinterpret_const_Multiarray_as_Matrix_d
	(const struct const_Multiarray_d* a, const struct const_Matrix_d* a_M, const ptrdiff_t ext_0,
	 const ptrdiff_t ext_1)
{
	reinterpret_Multiarray_as_Matrix_d((struct Multiarray_d*)a,(struct Matrix_d*)a_M,ext_0,ext_1);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

void reinterpret_Multiarray_as_Matrix_d
	(struct Multiarray_d* a, struct Matrix_d* a_M, const ptrdiff_t ext_0, const ptrdiff_t ext_1)
{
	assert(compute_size(a->order,a->extents) == ext_0*ext_1);
	a_M->layout    = a->layout;
	a_M->ext_0     = ext_0;
	a_M->ext_1     = ext_1;
	a_M->owns_data = false;
	a_M->data      = a->data;
}
