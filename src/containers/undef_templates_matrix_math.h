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
 *  \brief Undefine macro definitions for c-style templating relating to def_templates_matrix_math_\*.h.
 */



#undef gsl_matrix_T

#undef gsl_permute_matrix_T
#undef mkl_Timatcopy
#undef LAPACKE_Tsysv
#undef LAPACKE_Tsgesv



#undef compute_norm_Matrix_T_row
#undef transpose_Matrix_T
#undef invert_sub_block_Matrix_T
#undef scale_Matrix_T
#undef permute_Matrix_T
#undef permute_Matrix_T_V
#undef mm_T
#undef mm_RTT
#undef mm_TRT
#undef mv_T
#undef scale_Matrix_T_by_Vector_R
#undef mm_diag_T
#undef reinterpret_const_Matrix_T
