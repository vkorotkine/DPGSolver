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
 *  \brief Provides the macro definitions used for c-style templating related to the `int` vector constructor
 *         functions.
 */

///\{ \name Function names
#define constructor_default_Vector_T       constructor_default_Vector_i
#define constructor_default_const_Vector_T constructor_default_const_Vector_i
#define constructor_default_Vector_T_2     constructor_default_Vector_i_2

#define constructor_empty_Vector_T constructor_empty_Vector_i

#define constructor_zero_Vector_T constructor_zero_Vector_i

#define constructor_copy_Vector_T         constructor_copy_Vector_i
#define constructor_copy_const_Vector_T   constructor_copy_const_Vector_i
#define constructor_copy_Vector_T_T       constructor_copy_Vector_i_i
#define constructor_copy_const_Vector_T_T constructor_copy_const_Vector_i_i

#define constructor_move_Vector_T_T                constructor_move_Vector_i_i
#define constructor_move_const_Vector_T_T          constructor_move_const_Vector_i_i
#define constructor_move_const_Vector_Matrix_row_T constructor_move_const_Vector_Matrix_row_i
#define constructor_move_Vector_T_Matrix_T         constructor_move_Vector_i_Matrix_i
#define const_constructor_move_Vector_T            const_constructor_move_Vector_i
#define const_constructor_move_const_Vector_T      const_constructor_move_const_Vector_i

#define constructor_file_name_Vector_T       constructor_file_name_Vector_i
#define constructor_file_name_const_Vector_T constructor_file_name_const_Vector_i
#define constructor_file_Vector_T            constructor_file_Vector_i
#define constructor_file_const_Vector_T      constructor_file_const_Vector_i

#define destructor_Vector_T                   destructor_Vector_i
#define destructor_const_Vector_T             destructor_const_Vector_i
#define destructor_conditional_Vector_T       destructor_conditional_Vector_i
#define destructor_conditional_const_Vector_T destructor_conditional_const_Vector_i
///\}
