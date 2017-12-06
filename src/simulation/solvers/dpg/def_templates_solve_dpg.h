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

///\{ \name Data types
#define S_Params_DPG_T S_Params_DPG
///\}

///\{ \name Function pointers
#define constructor_norm_op_fptr_T constructor_norm_op_fptr
#define compute_rlhs_fptr_T        compute_rlhs_fptr
///\}

///\{ \name Function names
#define get_operator__cvt1_vt_vc__rlhs_T      get_operator__cvt1_vt_vc__rlhs
#define constructor_lhs_l_internal_face_dpg_T constructor_lhs_l_internal_face_dpg
#define compute_n_dof_nf_T                    compute_n_dof_nf
#define constructor_petsc_idxm_dpg_T                    constructor_petsc_idxm_dpg


#define set_s_params_dpg_T               set_s_params_dpg
#define get_operator__cv0_ff_fc_T        get_operator__cv0_ff_fc
#define set_idxm_T                       set_idxm
#define compute_rlhs_1_T                 compute_rlhs_1
#define constructor_norm_op__h1_upwind_T constructor_norm_op__h1_upwind
#define constructor_rhs_v_1_T            constructor_rhs_v_1
#define increment_and_add_dof_rlhs_f_1_T increment_and_add_dof_rlhs_f_1
#define increment_rhs_source_T           increment_rhs_source
#define increment_rlhs_internal_face_T   increment_rlhs_internal_face
#define increment_rlhs_boundary_face_T   increment_rlhs_boundary_face
#define scale_by_Jacobian_T              scale_by_Jacobian
#define increment_rhs_boundary_face_T    increment_rhs_boundary_face
#define increment_lhs_boundary_face_T    increment_lhs_boundary_face
///\}
