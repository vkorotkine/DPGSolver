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

#ifndef DPG__complex_matrix_h__INCLUDED
#define DPG__complex_matrix_h__INCLUDED
/** \file
 *  \brief Provides `complex` Matrix_\* containers and related functions.
 *
 *  Potentially relevant comments may be found in \ref matrix.h.
 */

#include <stddef.h>
#include <stdbool.h>
#include <complex.h>

#include "complex_matrix_constructors.h"
#include "complex_matrix_math.h"
#include "complex_matrix_print.h"

// Interface functions ********************************************************************************************** //

// double complex
#include "def_templates_type_dc.h"
#include "def_templates_matrix_c.h"
#include "matrix_T.h"
#include "undef_templates_type.h"
#include "undef_templates_matrix.h"

/// \brief `complex` version of \ref set_block_Matrix_T with a real input.
void set_block_Matrix_c_d
	(struct Matrix_c* a,                 ///< See brief.
	 const struct const_Matrix_d* a_sub, ///< See brief.
	 const ptrdiff_t row0,               ///< See brief.
	 const ptrdiff_t col0,               ///< See brief.
	 const char set_type                 ///< See brief.
	);

#endif // DPG__complex_matrix_h__INCLUDED
