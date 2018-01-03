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
 *  \brief Provides the macro definitions used for c-style templating related to the matrix math functions.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name External container/function names
#define gsl_matrix_T gsl_matrix

#define gsl_permute_matrix_T gsl_permute_matrix
#define mkl_Timatcopy        mkl_dimatcopy
#define LAPACKE_Tsysv        LAPACKE_dsysv
#define LAPACKE_Tsgesv       LAPACKE_dsgesv
///\}

///\{ \name Function names
#define compute_norm_Matrix_T_row  compute_norm_Matrix_d_row
#define transpose_Matrix_T         transpose_Matrix_d
#define invert_sub_block_Matrix_T  invert_sub_block_Matrix_d
#define scale_Matrix_T             scale_Matrix_d
#define permute_Matrix_T           permute_Matrix_d
#define permute_Matrix_T_V         permute_Matrix_d_V
#define mm_T                       mm_d
#define mm_RTT                     mm_d
#define mm_TRT                     mm_d
#define mv_T                       mv_d
#define scale_Matrix_T_by_Vector_R scale_Matrix_d_by_Vector_d
#define mm_diag_T                  mm_diag_d
#define reinterpret_const_Matrix_T reinterpret_const_Matrix_d

#define update_layout_Multiarray_Matrix_T update_layout_Multiarray_Matrix_d
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name External container/function names
#define gsl_matrix_T gsl_matrix_complex

#define gsl_permute_matrix_T gsl_permute_matrix_complex
#define mkl_Timatcopy        mkl_zimatcopy
#define LAPACKE_Tsysv        LAPACKE_zsysv
#define LAPACKE_Tsgesv       LAPACKE_zcgesv
///\}

///\{ \name Function names
#define compute_norm_Matrix_T_row  compute_norm_Matrix_c_row
#define transpose_Matrix_T         transpose_Matrix_c
#define invert_sub_block_Matrix_T  invert_sub_block_Matrix_c
#define scale_Matrix_T             scale_Matrix_c
#define permute_Matrix_T           permute_Matrix_c
#define permute_Matrix_T_V         permute_Matrix_c_V
#define mm_T                       mm_c
#define mm_RTT                     mm_dcc
#define mm_TRT                     mm_cdc
#define mv_T                       mv_c
#define scale_Matrix_T_by_Vector_R scale_Matrix_c_by_Vector_d
#define mm_diag_T                  mm_diag_c
#define reinterpret_const_Matrix_T reinterpret_const_Matrix_c

#define update_layout_Multiarray_Matrix_T update_layout_Multiarray_Matrix_c
///\}

#endif

///\{ \name Real Data types/Function names
#define permute_Matrix_R_V         permute_Matrix_d_V
#define scale_Matrix_R_by_Vector_R scale_Matrix_d_by_Vector_d
#define mm_R                       mm_d
///\}
