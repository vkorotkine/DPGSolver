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
 *  \brief Provides the macro definitions used for c-style templating related to the `int` matrix constructor
 *         functions.
 */

///\{ \name Function names
#define constructor_default_Matrix_T       constructor_default_Matrix_i
#define constructor_default_const_Matrix_T constructor_default_const_Matrix_i

#define constructor_empty_Matrix_T constructor_empty_Matrix_i

#define constructor_copy_Matrix_T_T             constructor_copy_Matrix_i_i
#define constructor_copy_const_Matrix_T_T       constructor_copy_const_Matrix_i_i
#define constructor_copy_Matrix_T               constructor_copy_Matrix_i
#define constructor_copy_const_Matrix_T         constructor_copy_const_Matrix_i
#define constructor_copy_extract_const_Matrix_T constructor_copy_extract_const_Matrix_i
#define const_constructor_copy_Matrix_T         const_constructor_copy_Matrix_i

#define constructor_move_Matrix_T_T           constructor_move_Matrix_i_i
#define constructor_move_const_Matrix_T_T     constructor_move_const_Matrix_i_i
#define const_constructor_move_Matrix_T       const_constructor_move_Matrix_i
#define const_constructor_move_const_Matrix_T const_constructor_move_const_Matrix_i

#define destructor_Matrix_T       destructor_Matrix_i
#define destructor_const_Matrix_T destructor_const_Matrix_i
///\}
