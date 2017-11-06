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

#include "complex_matrix.h"

#include <assert.h>

#include "macros.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

double complex* get_row_Matrix_c (const ptrdiff_t row, const struct Matrix_c* a)
{
	assert(a->layout == 'R');
	return &a->data[row*(a->ext_1)];
}

double complex* get_col_Matrix_c (const ptrdiff_t col, const struct Matrix_c* a)
{
	assert(a->layout == 'C');
	return &a->data[col*(a->ext_0)];
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
