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

#ifndef DPG__matrix_h__INCLUDED
#define DPG__matrix_h__INCLUDED
/** \file
 *  \brief Provides real Matrix_\* containers and related functions.
 *
 *  Potentially relevant comments may be found in \ref multiarray.h.
 *
 *  Matrices are 2D Multiarrays.
 */

#include <stddef.h>

#include "matrix_constructors.h"
#include "matrix_math.h"
#include "matrix_print.h"

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "def_templates_matrix.h"
#include "matrix_T.h"
#include "undef_templates_type.h"
#include "undef_templates_matrix.h"

#include "def_templates_type_i.h"
#include "def_templates_matrix_i.h"
#include "matrix_T.h"
#include "undef_templates_type.h"
#include "undef_templates_matrix.h"

/// \brief Swap the layout.
void swap_layout
	(char*const layout ///< Pointer to the layout variable.
	);

/// \brief Swap the layout and the extents.
void swap_layout_and_extents
	(char*const layout,     ///< Pointer to the layout variable.
	 ptrdiff_t*const ext_0, ///< Pointer to the ext_0 variable.
	 ptrdiff_t*const ext_1  ///< Pointer to the ext_1 variable.
	);

/** \brief Compute the opposite layout.
 *  \return See brief. */
char compute_opposite_layout
	(const char layout_i ///< The input layout.
	);

/** \brief See return.
 *  \return The index of a Matrix corresponding to the given row/column input. */
ptrdiff_t compute_index_Matrix
	(const ptrdiff_t i,     ///< The index in the `ext_0` direction.
	 const ptrdiff_t j,     ///< The index in the `ext_1` direction.
	 const ptrdiff_t ext_0, ///< Defined in \ref Matrix_T.
	 const ptrdiff_t ext_1, ///< Defined in \ref Matrix_T.
	 const char layout      ///< Defined in \ref Matrix_T.
	);

/** \brief Sort the input 'Row'-major \ref Matrix_T\* having ext_1 == DIM such that the output is ordered by fastest
 *         increasing 0th column, 1st column, etc.
 *  \return Optionally return indices or `NULL`.
 *
 *  In the case of indices being returned, the values are such that: sorted[ind[i]] == unsorted[i].
 */
const struct const_Vector_i* row_sort_DIM_Matrix_d
	(struct Matrix_d*const src, ///< The Matrix to be sorted.
	 const bool return_indices  ///< Flag for whether the sorted indices should be returned.
	);

#endif // DPG__matrix_h__INCLUDED
