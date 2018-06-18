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
 *  \brief Provides the macro definitions used for c-style templating related to the rlhs computing functions for the opg
 *         faces.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Data types.
#define Lhs_Operators_OPG_T Lhs_Operators_OPG
///\}

///\{ \name Function pointers
#define compute_rlhs_opg_f_fptr_T compute_rlhs_opg_f_fptr
///\}

///\{ \name Function names
#define compute_face_rlhs_opg_T compute_face_rlhs_opg
#define update_coef_nf_f_opg_T update_coef_nf_f_opg
#define constructor_Lhs_Operators_OPG_T constructor_Lhs_Operators_OPG
#define destructor_Lhs_Operators_OPG_T destructor_Lhs_Operators_OPG
///\}

///\{ \name Static names
#define S_Params_T S_Params_d
#define Num_Flux_T Num_Flux_d
#define set_s_params_T set_s_params_d
#define constructor_Numerical_Flux_OPG_T constructor_Numerical_Flux_OPG_d
#define constructor_Flux_Ref_OPG_T constructor_Flux_Ref_OPG_d
#define scale_by_Jacobian_e_T scale_by_Jacobian_e_d
#define constructor_Numerical_Flux_Input_data_opg_T constructor_Numerical_Flux_Input_data_opg_d
#define constructor_jump_test_s_fc_T constructor_jump_test_s_fc_d
#define compute_rhs_f_opg_dg_like_T compute_rhs_f_opg_dg_like_d
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Data types.
#define Lhs_Operators_OPG_T Lhs_Operators_OPG_c
///\}

///\{ \name Function pointers
#define compute_rlhs_opg_f_fptr_T compute_rlhs_opg_f_fptr_c
///\}

///\{ \name Function names
#define compute_face_rlhs_opg_T compute_face_rlhs_opg_c
#define update_coef_nf_f_opg_T update_coef_nf_f_opg_c
#define constructor_Lhs_Operators_OPG_T constructor_Lhs_Operators_OPG_c
#define destructor_Lhs_Operators_OPG_T destructor_Lhs_Operators_OPG_c
///\}

///\{ \name Static names
#define S_Params_T S_Params_T_c
#define Num_Flux_T Num_Flux_T_c
#define set_s_params_T set_s_params_T_c
#define constructor_Numerical_Flux_OPG_T constructor_Numerical_Flux_OPG_T_c
#define constructor_Flux_Ref_OPG_T constructor_Flux_Ref_OPG_T_c
#define scale_by_Jacobian_e_T scale_by_Jacobian_e_T_c
#define constructor_Numerical_Flux_Input_data_opg_T constructor_Numerical_Flux_Input_data_opg_T_c
#define constructor_jump_test_s_fc_T constructor_jump_test_s_fc_c
#define compute_rhs_f_opg_dg_like_T compute_rhs_f_opg_dg_like_c
///\}

#endif
