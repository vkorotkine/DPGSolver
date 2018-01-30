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
 *  \brief Provides the macro definitions used for c-style templating related to the multiarray math functions.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Function names
#define reinterpret_Multiarray_as_Matrix_T reinterpret_Multiarray_as_Matrix_d
#define reinterpret_Matrix_as_Multiarray_T reinterpret_Matrix_as_Multiarray_d

#define transpose_Multiarray_T                   transpose_Multiarray_d
#define scale_Multiarray_T                       scale_Multiarray_d
#define normalize_Multiarray_T                   normalize_Multiarray_d
#define permute_Multiarray_T                     permute_Multiarray_d
#define permute_Multiarray_T_V                   permute_Multiarray_d_V
#define scale_Multiarray_T_by_Vector_R           scale_Multiarray_by_Vector_d
#define add_in_place_Multiarray_T                add_in_place_Multiarray_d
#define multiply_in_place_Multiarray_TR          multiply_in_place_Multiarray_dd
#define subtract_in_place_Multiarray_T           subtract_in_place_Multiarray_d
#define mm_NNC_Multiarray_T                      mm_NNC_Multiarray_d
#define mm_NN1C_Multiarray_T                     mm_NN1C_Multiarray_d
#define mm_NN1C_overwrite_Multiarray_T           mm_NN1C_overwrite_Multiarray_d
#define reinterpret_const_Multiarray_as_Matrix_T reinterpret_const_Multiarray_as_Matrix_d
#define reinterpret_const_Matrix_as_Multiarray_T reinterpret_const_Matrix_as_Multiarray_d
#define compute_extents_mm_MMa                   compute_extents_mm_MMa
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function names
#define reinterpret_Multiarray_as_Matrix_T reinterpret_Multiarray_as_Matrix_c
#define reinterpret_Matrix_as_Multiarray_T reinterpret_Matrix_as_Multiarray_c

#define transpose_Multiarray_T                   transpose_Multiarray_c
#define scale_Multiarray_T                       scale_Multiarray_c
#define normalize_Multiarray_T                   normalize_Multiarray_c
#define permute_Multiarray_T                     permute_Multiarray_c
#define permute_Multiarray_T_V                   permute_Multiarray_c_V
#define scale_Multiarray_T_by_Vector_R           scale_Multiarray_c_by_Vector_d
#define add_in_place_Multiarray_T                add_in_place_Multiarray_c
#define multiply_in_place_Multiarray_TR          multiply_in_place_Multiarray_cd
#define subtract_in_place_Multiarray_T           subtract_in_place_Multiarray_c
#define mm_NNC_Multiarray_T                      mm_NNC_Multiarray_c
#define mm_NN1C_Multiarray_T                     mm_NN1C_Multiarray_c
#define mm_NN1C_overwrite_Multiarray_T           mm_NN1C_overwrite_Multiarray_c
#define reinterpret_const_Multiarray_as_Matrix_T reinterpret_const_Multiarray_as_Matrix_c
#define reinterpret_const_Matrix_as_Multiarray_T reinterpret_const_Matrix_as_Multiarray_c
#define compute_extents_mm_MMa                   compute_extents_mm_MMa
///\}

#endif

///\{ \name Function names
#define transpose_Multiarray_R transpose_Multiarray_d
#define normalize_Multiarray_R                   normalize_Multiarray_d
#define permute_Multiarray_R                     permute_Multiarray_d
///\}
