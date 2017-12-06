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
 *  \brief Undefine macro definitions for c-style templating relating to volume rlhs computing functions.
 */

///\{ \name Data types
#undef S_Params_Volume_Structor_T
#undef Flux_Ref_T
///\}

///\{ \name Function pointers
#undef constructor_sol_vc_fptr_T
#undef destructor_sol_vc_fptr_T
///\}

///\{ \name Function names
#undef set_S_Params_Volume_Structor_T
#undef constructor_Flux_Ref_vol_T
#undef destructor_Flux_Ref_T
#undef constructor_lhs_v_1_T
#undef get_operator__cv0_vs_vc_T
#undef get_operator__tw1_vt_vc_T

#undef constructor_sol_vc_interp_T
#undef constructor_sol_vc_col_T
#undef destructor_sol_vc_interp_T
#undef destructor_sol_vc_col_T
#undef constructor_Flux_Ref_T
#undef constructor_flux_ref_T
///\}
