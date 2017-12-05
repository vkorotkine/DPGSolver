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
 *  \brief Provides the macro definitions used for c-style templating related to the `double` vector constructor
 *         functions.
 */

///\{ \name Function names
#define constructor_default_Vector_T       constructor_default_Vector_d
#define constructor_default_const_Vector_T constructor_default_const_Vector_d
#define constructor_default_Vector_T_2     constructor_default_Vector_d_2

#define constructor_empty_Vector_T constructor_empty_Vector_d

#define constructor_zero_Vector_T constructor_zero_Vector_d

#define constructor_copy_Vector_T         constructor_copy_Vector_d
#define constructor_copy_Vector_T_T       constructor_copy_Vector_d_d
#define constructor_copy_const_Vector_T_T constructor_copy_const_Vector_d_d

#define constructor_move_Vector_T_T                constructor_move_Vector_d_d
#define constructor_move_const_Vector_T_T          constructor_move_const_Vector_d_d
#define constructor_move_const_Vector_Matrix_row_T constructor_move_const_Vector_Matrix_row_d
#define constructor_move_Vector_T_Matrix_T         constructor_move_Vector_d_Matrix_d
#define const_constructor_move_Vector_T            const_constructor_move_Vector_d
#define const_constructor_move_const_Vector_T      const_constructor_move_const_Vector_d

#define constructor_set_Vector_T_Multiarray_T       constructor_set_Vector_d_Multiarray_d
#define constructor_set_const_Vector_T_Multiarray_T constructor_set_const_Vector_d_Multiarray_d

#define constructor_inverse_Vector_T                  constructor_inverse_Vector_d
#define constructor_inverse_const_Vector_T            constructor_inverse_const_Vector_d
#define constructor_dot_mult_const_Vector_T           constructor_dot_mult_const_Vector_d
#define constructor_sum_Vector_T_const_Matrix_T       constructor_sum_Vector_d_const_Matrix_d
#define constructor_sum_const_Vector_T_const_Matrix_T constructor_sum_const_Vector_d_const_Matrix_d
#define constructor_mv_Vector_T                       constructor_mv_Vector_d
#define constructor_mv_const_Vector_T                 constructor_mv_const_Vector_d
#define constructor_sgesv_Vector_T                    constructor_sgesv_Vector_d
#define constructor_sgesv_const_Vector_T              constructor_sgesv_const_Vector_d
#define set_Vector_from_Matrix_T                      set_Vector_from_Matrix_d
#define set_const_Vector_from_Matrix_T                set_const_Vector_from_Matrix_d
#define set_Vector_from_Multiarray_T                  set_Vector_from_Multiarray_d
#define set_const_Vector_from_Multiarray_T            set_const_Vector_from_Multiarray_d

#define destructor_Vector_T       destructor_Vector_d
#define destructor_const_Vector_T destructor_const_Vector_d
///\}

///\{ \name Function aliases
#define constructor_dot_mult_const_Vector_R           constructor_dot_mult_const_Vector_d

#define destructor_const_Vector_R destructor_const_Vector_d
///\}
