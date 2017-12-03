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
 *  \brief Provides the macro definitions used for c-style templating related to the `double` matrix math
 *         functions.
 */

///\{ \name Function names
#define compute_norm_Matrix_T_row  compute_norm_Matrix_d_row
#define transpose_Matrix_T         transpose_Matrix_d
#define invert_sub_block_Matrix_T  invert_sub_block_Matrix_d
#define scale_Matrix_T             scale_Matrix_d
#define permute_Matrix_T           permute_Matrix_d
#define permute_Matrix_T_V         permute_Matrix_d_V
#define mm_RTT                     mm_d
#define mm_T                       mm_d
#define mv_T                       mv_d
#define scale_Matrix_by_Vector_T   scale_Matrix_by_Vector_d
#define mm_diag_T                  mm_diag_d
#define reinterpret_const_Matrix_T reinterpret_const_Matrix_d
///\}
