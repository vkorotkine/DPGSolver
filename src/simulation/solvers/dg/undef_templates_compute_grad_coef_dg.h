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
 *  \brief Undefine macro definitions for c-style templating relating to gradient coefficient computing functions for
 *         the dg scheme.
 */

///\{ \name Function names
#undef compute_grad_coef_dg_T
///\}

#undef compute_grad_coef_volumes
#undef compute_grad_coef_faces
#undef get_operator__cv1_vs_vc
#undef constructor_grad_xyz_p
#undef assert_numerical_solution_is_central
#undef constructor_diff_s_num_s
#undef constructor_jdet_n_diff_fc
#undef destructor_jdet_n_diff_fc
#undef compute_g_coef_f_i
#undef compute_d_g_coef_f__d_s_coef_i
#undef compute_g_coef_f_i_using_lin
#undef compute_g_coef_related_boundary
#undef add_face_grad_coef_f_to_volumes
#undef constructor_m_inv_tw0_vt_fc
#undef get_jn_fc_V
#undef get_sol_scale
