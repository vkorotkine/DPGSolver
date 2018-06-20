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
 *  \brief Undefine macro definitions for c-style templating relating to rlhs computing functions for the dg faces.
 */

///\{ \name Function names
#undef compute_face_rlhs_dg_T
#undef compute_flux_imbalances_faces_dg_T
#undef constructor_Numerical_Flux_Input_data_dg_T
///\}

#undef S_Params_T
#undef Num_Flux_T
#undef set_s_params_T
#undef add_to_flux_imbalance
#undef constructor_Boundary_Value_Input_g_face_fcl
#undef constructor_Boundary_Value_g_face_fcl
#undef constructor_partial_grad_fc_interp
#undef compute_scaling_weak_gradient
