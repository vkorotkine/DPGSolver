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

// Templated functions ********************************************************************************************** //

#include "def_templates_type_dc.h"
#include "def_templates_matrix_c.h"
#include "matrix_T.c"
#include "undef_templates_type.h"
#include "undef_templates_matrix.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void set_block_Matrix_c_d
	(struct Matrix_c* a, const struct const_Matrix_d* a_sub, const ptrdiff_t row0, const ptrdiff_t col0,
	 const char set_type)
{
	const struct const_Matrix_c*const a_sub_c = constructor_copy_const_Matrix_c_Matrix_d(a_sub); // destructed
	set_block_Matrix_c(a,a_sub_c,row0,col0,set_type);
	destructor_const_Matrix_c(a_sub_c);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
