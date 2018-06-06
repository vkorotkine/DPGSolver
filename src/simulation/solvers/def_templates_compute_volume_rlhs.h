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

#endif
