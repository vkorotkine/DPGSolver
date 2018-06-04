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
 *  \brief Provides the macro definitions used for c-style templating related to the face rlhs computing functions.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Function pointers
#define compute_rlhs_f_fptr_T compute_rlhs_f_fptr
///\}

///\{ \name Function names
#define get_operator__tw0_vt_fc_T               get_operator__tw0_vt_fc
#define get_operator__cv0_vs_fc_T               get_operator__cv0_vs_fc
#define get_operator__cv0_vr_fc_T               get_operator__cv0_vr_fc
#define permute_Matrix_T_fc                     permute_Matrix_d_fc
#define permute_Matrix_R_fc                     permute_Matrix_d_fc
#define get_operator__nc_fc_T                   get_operator__nc_fc
#define constructor_Numerical_Flux_Input_data_T constructor_Numerical_Flux_Input_data
#define destructor_Numerical_Flux_Input_data_T  destructor_Numerical_Flux_Input_data
#define constructor_lhs_f_1_T                   constructor_lhs_f_1
#define constructor_lhs_p_f_2_T                 constructor_lhs_p_f_2
#define add_to_flux_imbalance_face_nf_w_T       add_to_flux_imbalance_face_nf_w
#define compute_rhs_f_dg_like_T                 compute_rhs_f_dg_like
#define permute_Multiarray_T_fc                 permute_Multiarray_d_fc
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function pointers
#define compute_rlhs_f_fptr_T compute_rlhs_f_fptr_c
///\}

///\{ \name Function names
#define get_operator__tw0_vt_fc_T               get_operator__tw0_vt_fc_c
#define get_operator__cv0_vs_fc_T               get_operator__cv0_vs_fc_c
#define get_operator__cv0_vr_fc_T               get_operator__cv0_vr_fc_c
#define permute_Matrix_T_fc                     permute_Matrix_c_fc
#define permute_Matrix_R_fc                     permute_Matrix_d_fc_c
#define get_operator__nc_fc_T                   get_operator__nc_fc_c
#define constructor_Numerical_Flux_Input_data_T constructor_Numerical_Flux_Input_data_c
#define destructor_Numerical_Flux_Input_data_T  destructor_Numerical_Flux_Input_data_c
#define constructor_lhs_f_1_T                   constructor_lhs_f_1_c
#define constructor_lhs_p_f_2_T                 constructor_lhs_p_f_2_c
#define add_to_flux_imbalance_face_nf_w_T       add_to_flux_imbalance_face_nf_w_c
#define compute_rhs_f_dg_like_T                 compute_rhs_f_dg_like_c
#define permute_Multiarray_T_fc                 permute_Multiarray_c_fc
///\}

#endif
