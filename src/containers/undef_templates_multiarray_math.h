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
 *  \brief Undefine macro definitions for c-style templated relating to templates_multiarray_math_\*.h.
 */

#undef reinterpret_Multiarray_as_Matrix_T
#undef reinterpret_Matrix_as_Multiarray_T

#undef transpose_Multiarray_T
#undef scale_Multiarray_T
#undef normalize_Multiarray_T
#undef permute_Multiarray_T
#undef permute_Multiarray_T_V
#undef scale_Multiarray_T_by_Vector_R
#undef subtract_in_place_Multiarray_T
#undef mm_NNC_Multiarray_T
#undef mm_NN1C_Multiarray_T
#undef mm_NN1C_overwrite_Multiarray_T
#undef reinterpret_const_Multiarray_as_Matrix_T
#undef reinterpret_const_Matrix_as_Multiarray_T
#undef compute_extents_mm_MMa
