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
 *  \brief Provides the macro definitions used for c-style templating related to the vector constructor functions.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Function names
#define constructor_default_Vector_T       constructor_default_Vector_d
#define constructor_default_const_Vector_T constructor_default_const_Vector_d
#define constructor_default_Vector_T_2     constructor_default_Vector_d_2

#define constructor_empty_Vector_T constructor_empty_Vector_d

#define constructor_zero_Vector_T constructor_zero_Vector_d

#define constructor_copy_Vector_T                constructor_copy_Vector_d
#define constructor_copy_const_Vector_T          constructor_copy_const_Vector_d
#define constructor_copy_Vector_T_T              constructor_copy_Vector_d_d
#define constructor_copy_const_Vector_T_T        constructor_copy_const_Vector_d_d
#define constructor_copy_Vector_T_Vector_R       constructor_copy_Vector_d_Vector_d
#define constructor_copy_const_Vector_T_Vector_R constructor_copy_const_Vector_d_Vector_d

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
#define constructor_repeated_const_Vector_T           constructor_repeated_const_Vector_d
#define constructor_sum_Vectors_Vector_T              constructor_sum_Vectors_Vector_d
#define constructor_sum_Vectors_const_Vector_T        constructor_sum_Vectors_const_Vector_d
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

#define constructor_file_name_Vector_T       constructor_file_name_Vector_d
#define constructor_file_name_const_Vector_T constructor_file_name_const_Vector_d
#define constructor_file_Vector_T            constructor_file_Vector_d
#define constructor_file_const_Vector_T      constructor_file_const_Vector_d

#define destructor_Vector_T       destructor_Vector_d
#define destructor_const_Vector_T destructor_const_Vector_d
#define destructor_conditional_Vector_T       destructor_conditional_Vector_d
#define destructor_conditional_const_Vector_T destructor_conditional_const_Vector_d
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function names
#define constructor_default_Vector_T       constructor_default_Vector_c
#define constructor_default_const_Vector_T constructor_default_const_Vector_c
#define constructor_default_Vector_T_2     constructor_default_Vector_c_2

#define constructor_empty_Vector_T constructor_empty_Vector_c

#define constructor_zero_Vector_T constructor_zero_Vector_c

#define constructor_copy_Vector_T                constructor_copy_Vector_c
#define constructor_copy_const_Vector_T          constructor_copy_const_Vector_c
#define constructor_copy_Vector_T_T              constructor_copy_Vector_c_c
#define constructor_copy_const_Vector_T_T        constructor_copy_const_Vector_c_c
#define constructor_copy_Vector_T_Vector_R       constructor_copy_Vector_c_Vector_d
#define constructor_copy_const_Vector_T_Vector_R constructor_copy_const_Vector_c_Vector_d

#define constructor_move_Vector_T_T                constructor_move_Vector_c_c
#define constructor_move_const_Vector_T_T          constructor_move_const_Vector_c_c
#define constructor_move_const_Vector_Matrix_row_T constructor_move_const_Vector_Matrix_row_c
#define constructor_move_Vector_T_Matrix_T         constructor_move_Vector_c_Matrix_c
#define const_constructor_move_Vector_T            const_constructor_move_Vector_c
#define const_constructor_move_const_Vector_T      const_constructor_move_const_Vector_c

#define constructor_set_Vector_T_Multiarray_T       constructor_set_Vector_c_Multiarray_c
#define constructor_set_const_Vector_T_Multiarray_T constructor_set_const_Vector_c_Multiarray_c

#define constructor_inverse_Vector_T                  constructor_inverse_Vector_c
#define constructor_inverse_const_Vector_T            constructor_inverse_const_Vector_c
#define constructor_dot_mult_const_Vector_T           constructor_dot_mult_const_Vector_c
#define constructor_repeated_const_Vector_T           constructor_repeated_const_Vector_c
#define constructor_sum_Vectors_Vector_T              constructor_sum_Vectors_Vector_c
#define constructor_sum_Vectors_const_Vector_T        constructor_sum_Vectors_const_Vector_c
#define constructor_sum_Vector_T_const_Matrix_T       constructor_sum_Vector_c_const_Matrix_c
#define constructor_sum_const_Vector_T_const_Matrix_T constructor_sum_const_Vector_c_const_Matrix_c
#define constructor_mv_Vector_T                       constructor_mv_Vector_c
#define constructor_mv_const_Vector_T                 constructor_mv_const_Vector_c
#define constructor_sgesv_Vector_T                    constructor_sgesv_Vector_c
#define constructor_sgesv_const_Vector_T              constructor_sgesv_const_Vector_c
#define set_Vector_from_Matrix_T                      set_Vector_from_Matrix_c
#define set_const_Vector_from_Matrix_T                set_const_Vector_from_Matrix_c
#define set_Vector_from_Multiarray_T                  set_Vector_from_Multiarray_c
#define set_const_Vector_from_Multiarray_T            set_const_Vector_from_Multiarray_c

#define destructor_Vector_T       destructor_Vector_c
#define destructor_const_Vector_T destructor_const_Vector_c
#define destructor_conditional_Vector_T       destructor_conditional_Vector_c
#define destructor_conditional_const_Vector_T destructor_conditional_const_Vector_c
///\}

#endif

///\{ \name Real Data types/Function names
#define constructor_copy_const_Vector_R_R constructor_copy_const_Vector_d_d

#define constructor_inverse_const_Vector_R            constructor_inverse_const_Vector_d
#define constructor_dot_mult_const_Vector_R           constructor_dot_mult_const_Vector_d
#define constructor_repeated_const_Vector_R           constructor_repeated_const_Vector_d

#define destructor_const_Vector_R destructor_const_Vector_d
///\}
