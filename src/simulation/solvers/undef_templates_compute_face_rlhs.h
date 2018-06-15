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
 *  \brief Undefine macro definitions for c-style templating relating to face rlhs computing functions.
 */

///\{ \name Function pointers
#undef compute_rlhs_f_fptr_T
///\}

#undef get_operator__tw0_vt_fc_T
#undef get_operator__cv0_vs_fc_T
#undef get_operator__cv0_vr_fc_T
#undef get_operator__cv0_ff_fc_T
#undef permute_Matrix_T_fc
#undef permute_Matrix_R_fc
#undef get_operator__nc_fc_T
#undef constructor_Numerical_Flux_Input_data_T
#undef destructor_Numerical_Flux_Input_data_T
#undef constructor_lhs_f_1_T
#undef constructor_lhs_p_f_2_T
#undef add_to_flux_imbalance_face_nf_w_T
#undef compute_rhs_f_dg_like_T
#undef permute_Multiarray_T_fc
#undef constructor_Flux_Input_data_f_T
#undef destructor_Flux_Input_data_f_T

#undef finalize_face_rhs_dg_like_T
