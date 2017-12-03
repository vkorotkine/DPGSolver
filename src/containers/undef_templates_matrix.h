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
 *  \brief Undefine macro definitions for c-style templated relating to templates_matrix_\*.h.
 */

#include "undef_templates_matrix_constructors.h"
#include "undef_templates_matrix_math.h"
#include "undef_templates_matrix_print.h"

#undef Matrix_T
#undef Matrix_R
#undef const_Matrix_T
#undef const_Matrix_R

#undef Matrix_CSR_T
#undef const_Matrix_CSR_T

#undef set_value_fptr_T
#undef set_value_insert_T
#undef set_value_add_T

#undef get_row_Matrix_T
#undef get_row_const_Matrix_T
#undef get_col_Matrix_T
#undef get_col_const_Matrix_T
#undef get_slice_Matrix_T
#undef get_slice_const_Matrix_T
#undef get_val_Matrix_T
#undef get_val_const_Matrix_T

#undef set_row_Matrix_T
#undef set_col_Matrix_T
#undef set_to_value_Matrix_T
#undef set_col_to_val_Matrix_T

#undef set_block_Matrix_T
#undef set_block_Matrix_T_R
