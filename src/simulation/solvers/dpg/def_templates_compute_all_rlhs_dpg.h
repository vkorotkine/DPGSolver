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
 *  \brief Provides the macro definitions used for c-style templating related to the dpg rlhs computing functions.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Function names
#define compute_all_rlhs_dpg_T                compute_all_rlhs_dpg
#define constructor_lhs_l_internal_face_dpg_T constructor_lhs_l_internal_face_dpg
#define compute_n_dof_nf_T                    compute_n_dof_nf
#define constructor_petsc_idxm_dpg_T          constructor_petsc_idxm_dpg
#define add_to_rlhs__face_T                   add_to_rlhs__face
#define compute_flux_imbalances_faces_dpg_T   compute_flux_imbalances_faces_dpg
#define add_to_rlhs__face_boundary_T          add_to_rlhs__face_boundary
///\}

///\{ \name Static names
#define constructor_norm_DPG_fptr constructor_norm_DPG_fptr
#define compute_rlhs_dpg_fptr compute_rlhs_dpg_fptr
#define S_Params_DPG S_Params_DPG
#define Norm_DPG Norm_DPG
#define set_s_params_dpg set_s_params_dpg
#define set_idxm set_idxm
#define add_to_rlhs__face_internal add_to_rlhs__face_internal
#define add_to_flux_imbalance add_to_flux_imbalance
#define constructor_Numerical_Flux_dpg constructor_Numerical_Flux_dpg
#define scale_by_Jacobian scale_by_Jacobian
#define increment_rhs_boundary_face increment_rhs_boundary_face
#define increment_lhs_boundary_face increment_lhs_boundary_face
#define compute_rlhs_1 compute_rlhs_1
#define constructor_norm_DPG__h0 constructor_norm_DPG__h0
#define constructor_norm_DPG__h1 constructor_norm_DPG__h1
#define constructor_norm_DPG__h1_upwind constructor_norm_DPG__h1_upwind
#define set_exact_normal_flux set_exact_normal_flux
#define constructor_rhs_v_1 constructor_rhs_v_1
#define increment_rhs_source increment_rhs_source
#define constructor_rhs_opt_neg constructor_rhs_opt_neg
#define destructor_Norm_DPG destructor_Norm_DPG
#define get_operator__ones_coef_vt get_operator__ones_coef_vt
#define constructor_l_mult_M constructor_l_mult_M
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function names
#define compute_all_rlhs_dpg_T                compute_all_rhs_dpg_c
#define constructor_lhs_l_internal_face_dpg_T constructor_lhs_l_internal_face_dpg_c
#define compute_n_dof_nf_T                    compute_n_dof_nf_c
#define constructor_petsc_idxm_dpg_T          constructor_petsc_idxm_dpg_c
#define add_to_rlhs__face_T                   add_to_rlhs__face_c
#define compute_flux_imbalances_faces_dpg_T   compute_flux_imbalances_faces_dpg_c
#define add_to_rlhs__face_boundary_T          add_to_rlhs__face_boundary_c
///\}

///\{ \name Static names
#define constructor_norm_DPG_fptr constructor_norm_DPG_fptr_c
#define compute_rlhs_dpg_fptr compute_rlhs_dpg_fptr_c
#define S_Params_DPG S_Params_DPG_c
#define Norm_DPG Norm_DPG_c
#define set_s_params_dpg set_s_params_dpg_c
#define set_idxm set_idxm_c
#define add_to_rlhs__face_internal add_to_rlhs__face_internal_c
#define add_to_flux_imbalance add_to_flux_imbalance_c
#define constructor_Numerical_Flux_dpg constructor_Numerical_Flux_dpg_c
#define scale_by_Jacobian scale_by_Jacobian_c
#define increment_rhs_boundary_face increment_rhs_boundary_face_c
#define increment_lhs_boundary_face increment_lhs_boundary_face_c
#define compute_rlhs_1 compute_rlhs_1_c
#define constructor_norm_DPG__h0 constructor_norm_DPG__h0_c
#define constructor_norm_DPG__h1 constructor_norm_DPG__h1_c
#define constructor_norm_DPG__h1_upwind constructor_norm_DPG__h1_upwind_c
#define set_exact_normal_flux set_exact_normal_flux_c
#define constructor_rhs_v_1 constructor_rhs_v_1_c
#define increment_rhs_source increment_rhs_source_c
#define constructor_rhs_opt_neg constructor_rhs_opt_neg_c
#define destructor_Norm_DPG destructor_Norm_DPG_c
#define get_operator__ones_coef_vt get_operator__ones_coef_vt_c
#define constructor_l_mult_M constructor_l_mult_M_c
///\}

#endif
