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
 *  \brief Provides the macro definitions used for c-style templating related to the gradient coefficient computing
 *         functions for the dg scheme.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Function names
#define compute_grad_coef_dg_T compute_grad_coef_dg
///\}

#define compute_grad_coef_volumes compute_grad_coef_volumes
#define compute_grad_coef_faces compute_grad_coef_faces
#define get_operator__cv1_vs_vc get_operator__cv1_vs_vc
#define constructor_grad_xyz_p constructor_grad_xyz_p
#define assert_numerical_solution_is_central assert_numerical_solution_is_central
#define constructor_diff_s_num_s constructor_diff_s_num_s
#define constructor_jdet_n_diff_fc constructor_jdet_n_diff_fc
#define destructor_jdet_n_diff_fc destructor_jdet_n_diff_fc
#define compute_g_coef_f_i compute_g_coef_f_i
#define compute_d_g_coef_f__d_s_coef_i compute_d_g_coef_f__d_s_coef_i
#define compute_g_coef_f_i_using_lin compute_g_coef_f_i_using_lin
#define compute_g_coef_related_boundary compute_g_coef_related_boundary
#define add_face_grad_coef_f_to_volumes add_face_grad_coef_f_to_volumes
#define constructor_m_inv_tw0_vt_fc constructor_m_inv_tw0_vt_fc
#define get_jn_fc_V get_jn_fc_V
#define get_sol_scale get_sol_scale

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function names
#define compute_grad_coef_dg_T compute_grad_coef_dg_c
///\}

#define compute_grad_coef_volumes compute_grad_coef_volumes_c
#define compute_grad_coef_faces compute_grad_coef_faces_c
#define get_operator__cv1_vs_vc get_operator__cv1_vs_vc_c
#define constructor_grad_xyz_p constructor_grad_xyz_p_c
#define assert_numerical_solution_is_central assert_numerical_solution_is_central_c
#define constructor_diff_s_num_s constructor_diff_s_num_s_c
#define constructor_jdet_n_diff_fc constructor_jdet_n_diff_fc_c
#define destructor_jdet_n_diff_fc destructor_jdet_n_diff_fc_c
#define compute_g_coef_f_i compute_g_coef_f_i_c
#define compute_d_g_coef_f__d_s_coef_i compute_d_g_coef_f__d_s_coef_i_c
#define compute_g_coef_f_i_using_lin compute_g_coef_f_i_using_lin_c
#define compute_g_coef_related_boundary compute_g_coef_related_boundary_c
#define add_face_grad_coef_f_to_volumes add_face_grad_coef_f_to_volumes_c
#define constructor_m_inv_tw0_vt_fc constructor_m_inv_tw0_vt_fc_c
#define get_jn_fc_V get_jn_fc_V_c
#define get_sol_scale get_sol_scale_c

#endif
