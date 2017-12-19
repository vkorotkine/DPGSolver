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
#define get_operator__cvt1_vt_vc__rlhs_T      get_operator__cvt1_vt_vc__rlhs
#define constructor_lhs_l_internal_face_dpg_T constructor_lhs_l_internal_face_dpg
#define compute_n_dof_nf_T                    compute_n_dof_nf
#define constructor_petsc_idxm_dpg_T          constructor_petsc_idxm_dpg
#define add_to_rlhs__face_T                   add_to_rlhs__face
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function names
#define compute_all_rlhs_dpg_T                compute_all_rhs_dpg_c
#define get_operator__cvt1_vt_vc__rlhs_T      get_operator__cvt1_vt_vc__rlhs_c
#define constructor_lhs_l_internal_face_dpg_T constructor_lhs_l_internal_face_dpg_c
#define compute_n_dof_nf_T                    compute_n_dof_nf_c
#define constructor_petsc_idxm_dpg_T          constructor_petsc_idxm_dpg_c
#define add_to_rlhs__face_T                   add_to_rlhs__face_c
///\}

#endif
