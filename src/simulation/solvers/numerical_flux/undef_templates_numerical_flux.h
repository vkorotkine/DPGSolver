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
 *  \brief Undefine macro definitions for c-style templating relating to numerical flux containers/functions.
 */

///\{ \name Data types
#undef Numerical_Flux_Input_T
#undef mutable_Numerical_Flux_T
#undef Numerical_Flux_T

#undef Neigh_Info_NF_T
#undef m_Neigh_Info_NF_T
///\}

///\{ \name Function pointers
#undef compute_Numerical_Flux_fptr_T
///\}

///\{ \name Function names
#undef constructor_Numerical_Flux_Input_T
#undef destructor_Numerical_Flux_Input_T

#undef constructor_Numerical_Flux_T
#undef destructor_Numerical_Flux_T

#undef compute_Numerical_Flux_1_T
#undef compute_Numerical_Flux_2_T
#undef compute_Numerical_Flux_12_T
///\}

///\{ \name Function names (pde specific)
#undef compute_Numerical_Flux_T_advection_upwind
#undef compute_Numerical_Flux_T_advection_upwind_jacobian

#undef compute_Numerical_Flux_T_diffusion_br2
#undef compute_Numerical_Flux_T_diffusion_br2_jacobian

#undef compute_Numerical_Flux_T_euler_lax_friedrichs
#undef compute_Numerical_Flux_T_euler_lax_friedrichs_jacobian
#undef compute_Numerical_Flux_T_euler_roe_pike
#undef compute_Numerical_Flux_T_euler_roe_pike_jacobian
///\}
