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
 *  \brief Provides the macro definitions used for c-style templating related to the `int` matrix
 *         containers/functions.
 */

#include "def_templates_matrix_constructors_i.h"
#include "def_templates_matrix_print_i.h"

///\{ \name Data types
#define Matrix_T       Matrix_i
#define Matrix_R       Matrix_d
#define const_Matrix_T const_Matrix_i
#define const_Matrix_R const_Matrix_d

#define Matrix_CSR_T       Matrix_CSR_i
#define const_Matrix_CSR_T const_Matrix_CSR_i
///\}

///\{ \name Function names
#define set_value_fptr_T   set_value_fptr_i
#define set_value_insert_T set_value_insert_i
#define set_value_add_T    set_value_add_i

#define get_row_Matrix_T         get_row_Matrix_i
#define get_row_const_Matrix_T   get_row_const_Matrix_i
#define get_col_Matrix_T         get_col_Matrix_i
#define get_col_const_Matrix_T   get_col_const_Matrix_i
#define get_row_Matrix_R         get_row_Matrix_i
#define get_row_const_Matrix_R   get_row_const_Matrix_i
#define get_col_Matrix_R         get_col_Matrix_i
#define get_col_const_Matrix_R   get_col_const_Matrix_i
#define get_slice_Matrix_T       get_slice_Matrix_i
#define get_slice_const_Matrix_T get_slice_const_Matrix_i
#define get_val_Matrix_T         get_val_Matrix_i
#define get_val_const_Matrix_T   get_val_const_Matrix_i

#define set_row_Matrix_T        set_row_Matrix_i
#define set_col_Matrix_T        set_col_Matrix_i
#define set_to_value_Matrix_T   set_to_value_Matrix_i
#define set_col_to_val_Matrix_T set_col_to_val_Matrix_i

#define set_block_Matrix_T   set_block_Matrix_i
#define set_block_Matrix_T_R set_block_Matrix_i
///\}
