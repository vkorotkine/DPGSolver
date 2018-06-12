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
 *  \brief Provides the macro definitions used for c-style templating related to the rlhs computing functions for the dg
 *         faces.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Function names
#define compute_face_rlhs_dg_T                     compute_face_rlhs_dg
#define compute_flux_imbalances_faces_dg_T         compute_flux_imbalances_faces_dg
#define constructor_Numerical_Flux_Input_data_dg_T constructor_Numerical_Flux_Input_data_dg
///\}

#define scale_by_Jacobian_fptr_T scale_by_Jacobian_fptr_T
#define S_Params_T S_Params_T
#define Num_Flux_T Num_Flux_T
#define set_s_params_T set_s_params_T
#define add_to_flux_imbalance add_to_flux_imbalance
#define constructor_Boundary_Value_Input_g_face_fcl constructor_Boundary_Value_Input_g_face_fcl
#define constructor_Boundary_Value_g_face_fcl constructor_Boundary_Value_g_face_fcl
#define scale_by_Jacobian_e_T scale_by_Jacobian_e_T
#define constructor_partial_grad_fc_interp constructor_partial_grad_fc_interp
#define compute_scaling_weak_gradient compute_scaling_weak_gradient

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function names
#define compute_face_rlhs_dg_T                     compute_face_rlhs_dg_c
#define compute_flux_imbalances_faces_dg_T         compute_flux_imbalances_faces_dg_c
#define constructor_Numerical_Flux_Input_data_dg_T constructor_Numerical_Flux_Input_data_dg_c
///\}

#define scale_by_Jacobian_fptr_T scale_by_Jacobian_fptr_T_c
#define S_Params_T S_Params_T_c
#define Num_Flux_T Num_Flux_T_c
#define set_s_params_T set_s_params_T_c
#define add_to_flux_imbalance add_to_flux_imbalance_c
#define constructor_Boundary_Value_Input_g_face_fcl constructor_Boundary_Value_Input_g_face_fcl_c
#define constructor_Boundary_Value_g_face_fcl constructor_Boundary_Value_g_face_fcl_c
#define scale_by_Jacobian_e_T scale_by_Jacobian_e_T_c
#define constructor_partial_grad_fc_interp constructor_partial_grad_fc_interp_c
#define compute_scaling_weak_gradient compute_scaling_weak_gradient_c

#endif
