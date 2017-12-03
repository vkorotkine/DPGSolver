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
 *  \brief Provides the macro definitions used for c-style templating related to the `double complex` matrix
 *         containers/functions.
 */

#include "def_templates_matrix_constructors_c.h"
#include "def_templates_matrix_math_c.h"

///\{ \name Data types
#define Matrix_T       Matrix_c
#define Matrix_R       Matrix_d
#define const_Matrix_T const_Matrix_c
#define const_Matrix_R const_Matrix_d

#define Matrix_CSR_T       Matrix_CSR_c
#define const_Matrix_CSR_T const_Matrix_CSR_c
///\}

///\{ \name Function names
#define set_value_fptr_T   set_value_fptr_c
#define set_value_insert_T set_value_insert_c
#define set_value_add_T    set_value_add_c

#define get_row_Matrix_T         get_row_Matrix_c
#define get_row_const_Matrix_T   get_row_const_Matrix_c
#define get_col_Matrix_T         get_col_Matrix_c
#define get_col_const_Matrix_T   get_col_const_Matrix_c
#define get_slice_Matrix_T       get_slice_Matrix_c
#define get_slice_const_Matrix_T get_slice_const_Matrix_c
#define get_val_Matrix_T         get_val_Matrix_c
#define get_val_const_Matrix_T   get_val_const_Matrix_c

#define set_row_Matrix_T        set_row_Matrix_c
#define set_col_Matrix_T        set_col_Matrix_c
#define set_to_value_Matrix_T   set_to_value_Matrix_c
#define set_col_to_val_Matrix_T set_col_to_val_Matrix_c

#define set_block_Matrix_T   set_block_Matrix_c
#define set_block_Matrix_T_R set_block_Matrix_c_d
///\}
