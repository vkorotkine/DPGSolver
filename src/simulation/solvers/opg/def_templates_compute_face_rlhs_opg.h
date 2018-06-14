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

///\{ \name Function pointers
#define compute_rlhs_opg_f_fptr_T compute_rlhs_opg_f_fptr
///\}

///\{ \name Function names
#define compute_face_rlhs_opg_T compute_face_rlhs_opg
#define update_coef_nf_f_opg_T update_coef_nf_f_opg
///\}

///\{ \name Static names
#define S_Params_T S_Params_T
#define Num_Flux_T Num_Flux_T
#define set_s_params_T set_s_params_T
#define constructor_Numerical_Flux_OPG_T constructor_Numerical_Flux_OPG_T
#define constructor_Flux_OPG_T constructor_Flux_OPG_T
#define scale_by_Jacobian_e_T scale_by_Jacobian_e_T
#define constructor_Numerical_Flux_Input_data_opg_T constructor_Numerical_Flux_Input_data_opg_T
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function pointers
#define compute_rlhs_opg_f_fptr_T compute_rlhs_opg_f_fptr_c
///\}

///\{ \name Function names
#define compute_face_rlhs_opg_T compute_face_rlhs_opg_c
#define update_coef_nf_f_opg_T update_coef_nf_f_opg_c
///\}

///\{ \name Static names
#define S_Params_T S_Params_T_c
#define Num_Flux_T Num_Flux_T_c
#define set_s_params_T set_s_params_T_c
#define constructor_Numerical_Flux_OPG_T constructor_Numerical_Flux_OPG_T_c
#define constructor_Flux_OPG_T constructor_Flux_OPG_T_c
#define scale_by_Jacobian_e_T scale_by_Jacobian_e_T_c
#define constructor_Numerical_Flux_Input_data_opg_T constructor_Numerical_Flux_Input_data_opg_T_c
///\}

#endif
