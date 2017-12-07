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

///\{ \name Data types
#define Num_Flux_T Num_Flux
#define S_Params_T S_Params
///\}

///\{ \name Function pointers
#define scale_by_Jacobian_fptr_T scale_by_Jacobian_fptr
#define compute_rlhs_fptr_T      compute_rlhs_fptr
///\}

///\{ \name Function names
#define compute_face_rlhs_dg_T compute_face_rlhs_dg

#define set_s_params_T        set_s_params
#define scale_by_Jacobian_e_T scale_by_Jacobian_e
#define compute_rhs_f_dg_T    compute_rhs_f_dg

#define finalize_face_rhs_dg_T finalize_face_rhs_dg
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Data types
#define Num_Flux_T Num_Flux_c
#define S_Params_T S_Params_c
///\}

///\{ \name Function pointers
#define scale_by_Jacobian_fptr_T scale_by_Jacobian_fptr_c
#define compute_rlhs_fptr_T      compute_rlhs_fptr_c
///\}

///\{ \name Function names
#define compute_face_rlhs_dg_T compute_face_rlhs_dg_c

#define set_s_params_T        set_s_params_c
#define scale_by_Jacobian_e_T scale_by_Jacobian_e_c
#define compute_rhs_f_dg_T    compute_rhs_f_dg_c

#define finalize_face_rhs_dg_T finalize_face_rhs_dg_c
///\}

#endif
