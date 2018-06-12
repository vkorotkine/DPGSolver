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
 *  \brief Provides the macro definitions used for c-style templating related to the operator functions.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Function names
#define constructor_mm_NN1_Operator_Multiarray_T       constructor_mm_NN1_Operator_Multiarray_d
#define constructor_mm_NN1_Operator_const_Multiarray_T constructor_mm_NN1_Operator_const_Multiarray_d
#define constructor_mm_NN1_Operator_Multiarray_T_Multiarray_R constructor_mm_NN1_Operator_Multiarray_d_Multiarray_d
#define constructor_mm_NN1_Operator_const_Multiarray_T_Multiarray_R constructor_mm_NN1_Operator_const_Multiarray_d_Multiarray_d

#define set_ops_tp_n_rows_cols_T set_ops_tp_n_rows_cols

#define mm_NNC_Operator_Multiarray_T  mm_NNC_Operator_Multiarray_d
#define mm_NN1C_Operator_Multiarray_T mm_NN1C_Operator_Multiarray_d
#define mm_NN1_Operator_Multiarray_T  mm_NN1_Operator_Multiarray_d
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function names
#define constructor_mm_NN1_Operator_Multiarray_T       constructor_mm_NN1_Operator_Multiarray_c
#define constructor_mm_NN1_Operator_const_Multiarray_T constructor_mm_NN1_Operator_const_Multiarray_c
#define constructor_mm_NN1_Operator_Multiarray_T_Multiarray_R constructor_mm_NN1_Operator_Multiarray_c_Multiarray_d
#define constructor_mm_NN1_Operator_const_Multiarray_T_Multiarray_R constructor_mm_NN1_Operator_const_Multiarray_c_Multiarray_d

#define set_ops_tp_n_rows_cols_T set_ops_tp_n_rows_cols_c

#define mm_NNC_Operator_Multiarray_T  mm_NNC_Operator_Multiarray_c
#define mm_NN1C_Operator_Multiarray_T mm_NN1C_Operator_Multiarray_c
#define mm_NN1_Operator_Multiarray_T  mm_NN1_Operator_Multiarray_c
///\}

#endif


///\{ \name Real Data types/Function names
#define constructor_mm_NN1_Operator_Multiarray_R       constructor_mm_NN1_Operator_Multiarray_d
#define constructor_mm_NN1_Operator_const_Multiarray_R constructor_mm_NN1_Operator_const_Multiarray_d

#define mm_NN1C_Operator_Multiarray_R mm_NN1C_Operator_Multiarray_d
///\}
