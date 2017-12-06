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

///\{ \name Data types
#undef Num_Flux_T
#undef S_Params_T
///\}

///\{ \name Function pointers
#undef scale_by_Jacobian_fptr_T
#undef compute_rlhs_fptr_T
///\}

///\{ \name Function names
#undef compute_face_rlhs_dg_T

#undef set_s_params_T
#undef scale_by_Jacobian_e_T
#undef compute_rhs_f_dg_T

#undef finalize_face_rhs_dg_T
///\}
