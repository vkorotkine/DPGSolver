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
 *  \brief Provides the macro definitions used for c-style templating related to the volume rlhs computing functions.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Data types
#define S_Params_Volume_Structor_T S_Params_Volume_Structor
#define Flux_Ref_T                 Flux_Ref
///\}

///\{ \name Function pointers
#define constructor_sol_vc_fptr_T constructor_sol_vc_fptr
#define destructor_sol_vc_fptr_T  destructor_sol_vc_fptr
#define compute_rlhs_v_fptr_T     compute_rlhs_v_fptr
///\}

///\{ \name Function names
#define set_S_Params_Volume_Structor_T set_S_Params_Volume_Structor
#define constructor_Flux_Ref_vol_T     constructor_Flux_Ref_vol
#define destructor_Flux_Ref_T          destructor_Flux_Ref
#define compute_rhs_v_dg_like_T        compute_rhs_v_dg_like
#define constructor_lhs_v_1_T          constructor_lhs_v_1
#define constructor_lhs_p_v_2_T        constructor_lhs_p_v_2
#define get_operator__cv0_vs_vc_T      get_operator__cv0_vs_vc
#define get_operator__cv0_vt_vc_T      get_operator__cv0_vt_vc
#define get_operator__cv0_vr_vc_T      get_operator__cv0_vr_vc
#define get_operator__tw1_vt_vc_T      get_operator__tw1_vt_vc
#define get_operator__cv1_vt_vc_T      get_operator__cv1_vt_vc
///\}

///\{ \name Static names
#define constructor_NULL constructor_NULL
#define constructor_sol_vc_interp constructor_sol_vc_interp
#define constructor_grad_vc_interp constructor_grad_vc_interp
#define constructor_sol_vc_col constructor_sol_vc_col
#define constructor_grad_vc_col constructor_grad_vc_col
#define constructor_xyz_vc constructor_xyz_vc
#define destructor_NULL destructor_NULL
#define destructor_sol_vc_interp destructor_sol_vc_interp
#define destructor_sol_vc_col destructor_sol_vc_col
#define constructor_Flux_Ref constructor_Flux_Ref
#define constructor_flux_ref_T constructor_flux_ref_T
#define get_operator__cv0_vg_vc_T get_operator__cv0_vg_vc_T
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Data types
#define S_Params_Volume_Structor_T S_Params_Volume_Structor_c
#define Flux_Ref_T                 Flux_Ref_c
///\}

///\{ \name Function pointers
#define constructor_sol_vc_fptr_T constructor_sol_vc_fptr_c
#define destructor_sol_vc_fptr_T  destructor_sol_vc_fptr_c
#define compute_rlhs_v_fptr_T     compute_rlhs_v_fptr_c
///\}

///\{ \name Function names
#define set_S_Params_Volume_Structor_T set_S_Params_Volume_Structor_c
#define constructor_Flux_Ref_vol_T     constructor_Flux_Ref_vol_c
#define destructor_Flux_Ref_T          destructor_Flux_Ref_c
#define compute_rhs_v_dg_like_T        compute_rhs_v_dg_like_c
#define constructor_lhs_v_1_T          constructor_lhs_v_1_c
#define constructor_lhs_p_v_2_T        constructor_lhs_p_v_2_c
#define get_operator__cv0_vs_vc_T      get_operator__cv0_vs_vc_c
#define get_operator__cv0_vt_vc_T      get_operator__cv0_vt_vc_c
#define get_operator__cv0_vr_vc_T      get_operator__cv0_vr_vc_c
#define get_operator__tw1_vt_vc_T      get_operator__tw1_vt_vc_c
#define get_operator__cv1_vt_vc_T      get_operator__cv1_vt_vc_c
///\}

///\{ \name Static names
#define constructor_NULL constructor_NULL_c
#define constructor_sol_vc_interp constructor_sol_vc_interp_c
#define constructor_grad_vc_interp constructor_grad_vc_interp_c
#define constructor_sol_vc_col constructor_sol_vc_col_c
#define constructor_grad_vc_col constructor_grad_vc_col_c
#define constructor_xyz_vc constructor_xyz_vc_c
#define destructor_NULL destructor_NULL_c
#define destructor_sol_vc_interp destructor_sol_vc_interp_c
#define destructor_sol_vc_col destructor_sol_vc_col_c
#define constructor_Flux_Ref constructor_Flux_Ref_c
#define constructor_flux_ref_T constructor_flux_ref_T_c
#define get_operator__cv0_vg_vc_T get_operator__cv0_vg_vc_T_c
///\}

#endif
