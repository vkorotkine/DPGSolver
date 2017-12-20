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
 *  \brief Provides the macro definitions used for c-style templating related to the matrix constructor functions.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Function names
#define constructor_default_Matrix_T       constructor_default_Matrix_d
#define constructor_default_const_Matrix_T constructor_default_const_Matrix_d

#define constructor_empty_Matrix_T constructor_empty_Matrix_d
#define constructor_zero_Matrix_T  constructor_zero_Matrix_d

#define constructor_copy_Matrix_T_T              constructor_copy_Matrix_d_d
#define constructor_copy_const_Matrix_T_T        constructor_copy_const_Matrix_d_d
#define constructor_copy_Matrix_T                constructor_copy_Matrix_d
#define constructor_copy_const_Matrix_T          constructor_copy_const_Matrix_d
#define constructor_copy_extract_const_Matrix_T  constructor_copy_extract_const_Matrix_d
#define const_constructor_copy_Matrix_T          const_constructor_copy_Matrix_d
#define constructor_copy_Matrix_T_Matrix_R       constructor_copy_Matrix_d_Matrix_d
#define constructor_copy_const_Matrix_T_Matrix_R constructor_copy_const_Matrix_d_Matrix_d

#define constructor_move_Matrix_T_T           constructor_move_Matrix_d_d
#define constructor_move_const_Matrix_T_T     constructor_move_const_Matrix_d_d
#define const_constructor_move_Matrix_T       const_constructor_move_Matrix_d
#define const_constructor_move_const_Matrix_T const_constructor_move_const_Matrix_d

#define constructor_sub_block_Matrix_T            constructor_sub_block_Matrix_d
#define constructor_subset_const_Matrix_T         constructor_subset_const_Matrix_d
#define constructor_copy_transpose_Matrix_T       constructor_copy_transpose_Matrix_d
#define constructor_block_diagonal_const_Matrix_T constructor_block_diagonal_const_Matrix_d
#define constructor_diagonal_Matrix_T_T           constructor_diagonal_Matrix_d_d
#define constructor_identity_Matrix_T             constructor_identity_Matrix_d
#define constructor_identity_const_Matrix_T       constructor_identity_const_Matrix_d
#define constructor_inverse_Matrix_T              constructor_inverse_Matrix_d
#define constructor_inverse_const_Matrix_T        constructor_inverse_const_Matrix_d
#define constructor_sgesv_Matrix_T                constructor_sgesv_Matrix_d
#define constructor_sgesv_const_Matrix_T          constructor_sgesv_const_Matrix_d
#define constructor_sysv_Matrix_T                 constructor_sysv_Matrix_d
#define constructor_sysv_const_Matrix_T           constructor_sysv_const_Matrix_d
#define constructor_mm_Matrix_T                   constructor_mm_Matrix_d
#define constructor_mm_const_Matrix_T             constructor_mm_const_Matrix_d
#define constructor_mm_NN1R_Matrix_T              constructor_mm_NN1R_Matrix_d
#define constructor_mm_NN1R_const_Matrix_T        constructor_mm_NN1R_const_Matrix_d
#define constructor_mm_NN1C_Matrix_T              constructor_mm_NN1C_Matrix_d
#define constructor_mm_NN1C_const_Matrix_T        constructor_mm_NN1C_const_Matrix_d
#define constructor_mm_diag_Matrix_T_R            constructor_mm_diag_Matrix_d
#define constructor_mm_diag_const_Matrix_T_R      constructor_mm_diag_const_Matrix_d
#define set_Matrix_from_Multiarray_T              set_Matrix_from_Multiarray_d
#define set_const_Matrix_from_Multiarray_T        set_const_Matrix_from_Multiarray_d
#define set_Matrix_from_Multiarray_Matrix_T       set_Matrix_from_Multiarray_Matrix_d
#define set_const_Matrix_from_Multiarray_Matrix_T set_const_Matrix_from_Multiarray_Matrix_d

#define destructor_Matrix_T                   destructor_Matrix_d
#define destructor_const_Matrix_T             destructor_const_Matrix_d
#define destructor_conditional_Matrix_T       destructor_conditional_Matrix_d
#define destructor_conditional_const_Matrix_T destructor_conditional_const_Matrix_d
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function names
#define constructor_default_Matrix_T       constructor_default_Matrix_c
#define constructor_default_const_Matrix_T constructor_default_const_Matrix_c

#define constructor_empty_Matrix_T constructor_empty_Matrix_c
#define constructor_zero_Matrix_T  constructor_zero_Matrix_c

#define constructor_copy_Matrix_T_T             constructor_copy_Matrix_c_c
#define constructor_copy_const_Matrix_T_T       constructor_copy_const_Matrix_c_c
#define constructor_copy_Matrix_T               constructor_copy_Matrix_c
#define constructor_copy_const_Matrix_T         constructor_copy_const_Matrix_c
#define constructor_copy_extract_const_Matrix_T constructor_copy_extract_const_Matrix_c
#define const_constructor_copy_Matrix_T         const_constructor_copy_Matrix_c
#define constructor_copy_Matrix_T_Matrix_R       constructor_copy_Matrix_c_Matrix_d
#define constructor_copy_const_Matrix_T_Matrix_R constructor_copy_const_Matrix_c_Matrix_d

#define constructor_move_Matrix_T_T           constructor_move_Matrix_c_c
#define constructor_move_const_Matrix_T_T     constructor_move_const_Matrix_c_c
#define const_constructor_move_Matrix_T       const_constructor_move_Matrix_c
#define const_constructor_move_const_Matrix_T const_constructor_move_const_Matrix_c

#define constructor_sub_block_Matrix_T            constructor_sub_block_Matrix_c
#define constructor_subset_const_Matrix_T         constructor_subset_const_Matrix_c
#define constructor_copy_transpose_Matrix_T       constructor_copy_transpose_Matrix_c
#define constructor_block_diagonal_const_Matrix_T constructor_block_diagonal_const_Matrix_c
#define constructor_diagonal_Matrix_T_T           constructor_diagonal_Matrix_c_c
#define constructor_identity_Matrix_T             constructor_identity_Matrix_c
#define constructor_identity_const_Matrix_T       constructor_identity_const_Matrix_c
#define constructor_inverse_Matrix_T              constructor_inverse_Matrix_c
#define constructor_inverse_const_Matrix_T        constructor_inverse_const_Matrix_c
#define constructor_sgesv_Matrix_T                constructor_sgesv_Matrix_c
#define constructor_sgesv_const_Matrix_T          constructor_sgesv_const_Matrix_c
#define constructor_sysv_Matrix_T                 constructor_sysv_Matrix_c
#define constructor_sysv_const_Matrix_T           constructor_sysv_const_Matrix_c
#define constructor_mm_Matrix_T                   constructor_mm_Matrix_c
#define constructor_mm_const_Matrix_T             constructor_mm_const_Matrix_c
#define constructor_mm_NN1R_Matrix_T              constructor_mm_NN1R_Matrix_c
#define constructor_mm_NN1R_const_Matrix_T        constructor_mm_NN1R_const_Matrix_c
#define constructor_mm_NN1C_Matrix_T              constructor_mm_NN1C_Matrix_c
#define constructor_mm_NN1C_const_Matrix_T        constructor_mm_NN1C_const_Matrix_c
#define constructor_mm_diag_Matrix_T_R            constructor_mm_diag_Matrix_c_d
#define constructor_mm_diag_const_Matrix_T_R      constructor_mm_diag_const_Matrix_c_d
#define set_Matrix_from_Multiarray_T              set_Matrix_from_Multiarray_c
#define set_const_Matrix_from_Multiarray_T        set_const_Matrix_from_Multiarray_c
#define set_Matrix_from_Multiarray_Matrix_T       set_Matrix_from_Multiarray_Matrix_c
#define set_const_Matrix_from_Multiarray_Matrix_T set_const_Matrix_from_Multiarray_Matrix_c

#define destructor_Matrix_T                   destructor_Matrix_c
#define destructor_const_Matrix_T             destructor_const_Matrix_c
#define destructor_conditional_Matrix_T       destructor_conditional_Matrix_c
#define destructor_conditional_const_Matrix_T destructor_conditional_const_Matrix_c
///\}

#endif

///\{ \name Real Data types/Function names
#define constructor_copy_Matrix_R                constructor_copy_Matrix_d
#define constructor_copy_const_Matrix_R          constructor_copy_const_Matrix_d

#define constructor_inverse_const_Matrix_R        constructor_inverse_const_Matrix_d
#define constructor_mm_const_Matrix_R             constructor_mm_const_Matrix_d
#define constructor_mm_diag_const_Matrix_R        constructor_mm_diag_const_Matrix_d

#define destructor_Matrix_R       destructor_Matrix_d
#define destructor_const_Matrix_R destructor_const_Matrix_d
///\}
