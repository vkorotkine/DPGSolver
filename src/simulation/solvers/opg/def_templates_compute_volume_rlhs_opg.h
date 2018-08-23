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
 *         volumes.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Function names
#define compute_volume_rlhs_opg_T compute_volume_rlhs_opg
#define constructor_operator__test_s_coef_to_sol_coef_T constructor_operator__test_s_coef_to_sol_coef_d
#define update_coef_s_v_opg_T update_coef_s_v_opg_d
#define constructor_test_diff_op_1v_opg_T constructor_test_diff_op_1v_opg_d
///\}

///\{ \name Static names
#define S_Params_T S_Params
#define set_s_params_T set_s_params_T
#define constructor_Flux_Ref_vol_opg_T constructor_Flux_Ref_vol_opg_d
#define constructor_sol_vs_T constructor_sol_vs_d
#define constructor_grad_vs_T constructor_grad_vs_d
#define constructor_xyz_vs_T constructor_xyz_vs_d
#define compute_rhs_v_opg_T compute_rhs_v_opg_d
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function names
#define compute_volume_rlhs_opg_T compute_volume_rlhs_opg_c
#define constructor_operator__test_s_coef_to_sol_coef_T constructor_operator__test_s_coef_to_sol_coef_c
#define update_coef_s_v_opg_T update_coef_s_v_opg_c
#define constructor_test_diff_op_1v_opg_T constructor_test_diff_op_1v_opg_c
///\}

///\{ \name Static names
#define S_Params_T S_Params_T_c
#define set_s_params_T set_s_params_T_c
#define constructor_Flux_Ref_vol_opg_T constructor_Flux_Ref_vol_opg_c
#define constructor_sol_vs_T constructor_sol_vs_c
#define constructor_grad_vs_T constructor_grad_vs_c
#define constructor_xyz_vs_T constructor_xyz_vs_c
#define compute_rhs_v_opg_T compute_rhs_v_opg_c
///\}

#endif
