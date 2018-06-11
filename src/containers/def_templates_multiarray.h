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
 *  \brief Provides the macro definitions used for c-style templating related to the multiarray containers/functions.
 */

#include "def_templates_multiarray_constructors.h"
#include "def_templates_multiarray_math.h"
#include "def_templates_multiarray_print.h"

#if defined TYPE_RC
#if TYPE_RC == TYPE_REAL

///\{ \name Data types
#define Multiarray_T       Multiarray_d
#define const_Multiarray_T const_Multiarray_d

#define Multiarray_Vector_T       Multiarray_Vector_d
#define const_Multiarray_Vector_T const_Multiarray_Vector_d
#define Multiarray_Matrix_T       Multiarray_Matrix_d
#define const_Multiarray_Matrix_T const_Multiarray_Matrix_d
///\}

///\{ \name Function names
#define Vector_T_indexed Vector_d_indexed

#define constructor_move_Vector_T_indexed constructor_move_Vector_d_indexed
#define destructor_Vector_T_indexed       destructor_Vector_d_indexed
#define cmp_Vector_T_indexed              cmp_Vector_d_indexed
#define reorder_Multiarray_Vector_T       reorder_Multiarray_Vector_d
#define compute_total_entries             compute_total_entries_d

#define get_row_Multiarray_T          get_row_Multiarray_d
#define get_row_const_Multiarray_T    get_row_const_Multiarray_d
#define get_col_Multiarray_T          get_col_Multiarray_d
#define get_col_const_Multiarray_T    get_col_const_Multiarray_d
#define remove_col_Multiarray_T       remove_col_Multiarray_d
#define set_to_value_Multiarray_T     set_to_value_Multiarray_d
#define set_Multiarray_Vector_T_T     set_Multiarray_Vector_d_d
#define set_Multiarray_T              set_Multiarray_d
#define set_Multiarray_T_Multiarray_R set_Multiarray_d_Multiarray_d

#define sort_Multiarray_Vector_T               sort_Multiarray_Vector_d
#define collapse_Multiarray_Vector_T           collapse_Multiarray_Vector_d
#define resize_Multiarray_T                    resize_Multiarray_d
#define get_const_Multiarray_Vector_T          get_const_Multiarray_Vector_d
#define interpret_Multiarray_as_Vector_T       interpret_Multiarray_as_Vector_d
#define interpret_const_Multiarray_as_Vector_T interpret_const_Multiarray_as_Vector_d
#define interpret_Multiarray_as_Matrix_T       interpret_Multiarray_as_Matrix_d
#define interpret_const_Multiarray_as_Matrix_T interpret_const_Multiarray_as_Matrix_d
#define interpret_Multiarray_as_slice_T        interpret_Multiarray_as_slice_d
#define interpret_const_Multiarray_as_slice_T  interpret_const_Multiarray_as_slice_d
#define interpret_Multiarray_slice_as_Vector_T interpret_Multiarray_slice_as_Vector_d
#define copy_into_Multiarray_T                 copy_into_Multiarray_d
#define copy_into_Multiarray_T_from_R          copy_into_Multiarray_d_from_d
#define push_back_Multiarray_T                 push_back_Multiarray_d
#define make_unique_row_Multiarray_T           make_unique_row_Multiarray_d
#define update_rows_Multiarray_T               update_rows_Multiarray_d
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Data types
#define Multiarray_T       Multiarray_c
#define const_Multiarray_T const_Multiarray_c

#define Multiarray_Vector_T       Multiarray_Vector_c
#define const_Multiarray_Vector_T const_Multiarray_Vector_c
#define Multiarray_Matrix_T       Multiarray_Matrix_c
#define const_Multiarray_Matrix_T const_Multiarray_Matrix_c
///\}

///\{ \name Function names
#define Vector_T_indexed Vector_c_indexed

#define constructor_move_Vector_T_indexed constructor_move_Vector_c_indexed
#define destructor_Vector_T_indexed       destructor_Vector_c_indexed
#define cmp_Vector_T_indexed              cmp_Vector_c_indexed
#define reorder_Multiarray_Vector_T       reorder_Multiarray_Vector_c
#define compute_total_entries             compute_total_entries_c

#define get_row_Multiarray_T          get_row_Multiarray_c
#define get_row_const_Multiarray_T    get_row_const_Multiarray_c
#define get_col_Multiarray_T          get_col_Multiarray_c
#define get_col_const_Multiarray_T    get_col_const_Multiarray_c
#define remove_col_Multiarray_T       remove_col_Multiarray_c
#define set_to_value_Multiarray_T     set_to_value_Multiarray_c
#define set_Multiarray_Vector_T_T     set_Multiarray_Vector_c_c
#define set_Multiarray_T              set_Multiarray_c
#define set_Multiarray_T_Multiarray_R set_Multiarray_c_Multiarray_d

#define sort_Multiarray_Vector_T               sort_Multiarray_Vector_c
#define collapse_Multiarray_Vector_T           collapse_Multiarray_Vector_c
#define resize_Multiarray_T                    resize_Multiarray_c
#define get_const_Multiarray_Vector_T          get_const_Multiarray_Vector_c
#define interpret_Multiarray_as_Vector_T       interpret_Multiarray_as_Vector_c
#define interpret_const_Multiarray_as_Vector_T interpret_const_Multiarray_as_Vector_c
#define interpret_Multiarray_as_Matrix_T       interpret_Multiarray_as_Matrix_c
#define interpret_const_Multiarray_as_Matrix_T interpret_const_Multiarray_as_Matrix_c
#define interpret_Multiarray_as_slice_T        interpret_Multiarray_as_slice_c
#define interpret_const_Multiarray_as_slice_T  interpret_const_Multiarray_as_slice_c
#define interpret_Multiarray_slice_as_Vector_T interpret_Multiarray_slice_as_Vector_c
#define copy_into_Multiarray_T                 copy_into_Multiarray_c
#define copy_into_Multiarray_T_from_R          copy_into_Multiarray_c_from_d
#define push_back_Multiarray_T                 push_back_Multiarray_c
#define make_unique_row_Multiarray_T           make_unique_row_Multiarray_c
#define update_rows_Multiarray_T               update_rows_Multiarray_c
///\}

#endif

///\{ \name Real Data types/Function names
#define Multiarray_R       Multiarray_d
#define const_Multiarray_R const_Multiarray_d

#define get_row_Multiarray_R          get_row_Multiarray_d
#define get_row_const_Multiarray_R    get_row_const_Multiarray_d
#define get_col_Multiarray_R          get_col_Multiarray_d
#define get_col_const_Multiarray_R    get_col_const_Multiarray_d
#define set_to_value_Multiarray_R     set_to_value_Multiarray_d
#define set_Multiarray_R              set_Multiarray_d

#define resize_Multiarray_R                    resize_Multiarray_d
#define interpret_const_Multiarray_as_Vector_R interpret_const_Multiarray_as_Vector_d
#define interpret_Multiarray_as_Matrix_R       interpret_Multiarray_as_Matrix_d
#define interpret_const_Multiarray_as_Matrix_R interpret_const_Multiarray_as_Matrix_d
#define update_rows_Multiarray_R               update_rows_Multiarray_d
///\}



#elif defined TYPE_I
#if TYPE_I == TYPE_II

///\{ \name Data types
#define Multiarray_T       Multiarray_i
#define Multiarray_R       Multiarray_d
#define const_Multiarray_T const_Multiarray_i
#define const_Multiarray_R const_Multiarray_d

#define Multiarray_Vector_T       Multiarray_Vector_i
#define const_Multiarray_Vector_T const_Multiarray_Vector_i
#define Multiarray_Matrix_T       Multiarray_Matrix_i
#define const_Multiarray_Matrix_T const_Multiarray_Matrix_i
///\}

///\{ \name Function names
#define Vector_T_indexed Vector_i_indexed

#define constructor_move_Vector_T_indexed constructor_move_Vector_i_indexed
#define destructor_Vector_T_indexed       destructor_Vector_i_indexed
#define cmp_Vector_T_indexed              cmp_Vector_i_indexed
#define reorder_Multiarray_Vector_T       reorder_Multiarray_Vector_i
#define compute_total_entries             compute_total_entries_i

#define get_row_Multiarray_T          get_row_Multiarray_i
#define get_row_const_Multiarray_T    get_row_const_Multiarray_i
#define get_col_Multiarray_T          get_col_Multiarray_i
#define get_col_const_Multiarray_T    get_col_const_Multiarray_i
#define set_to_value_Multiarray_T     set_to_value_Multiarray_i
#define set_Multiarray_Vector_T_T     set_Multiarray_Vector_i_i
#define set_Multiarray_T              set_Multiarray_i
#define set_Multiarray_T_Multiarray_R set_Multiarray_i_Multiarray_d

#define sort_Multiarray_Vector_T               sort_Multiarray_Vector_i
#define collapse_Multiarray_Vector_T           collapse_Multiarray_Vector_i
#define resize_Multiarray_T                    resize_Multiarray_i
#define get_const_Multiarray_Vector_T          get_const_Multiarray_Vector_i
#define interpret_const_Multiarray_as_Vector_T interpret_const_Multiarray_as_Vector_i
#define interpret_const_Multiarray_as_Vector_R interpret_const_Multiarray_as_Vector_i
#define interpret_Multiarray_as_Matrix_T       interpret_Multiarray_as_Matrix_i
///\}

#endif

#endif
