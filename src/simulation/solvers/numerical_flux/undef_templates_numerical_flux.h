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

#undef compute_Numerical_Flux_T_diffusion_central
#undef compute_Numerical_Flux_T_diffusion_central_jacobian

#undef compute_Numerical_Flux_T_euler_lax_friedrichs
#undef compute_Numerical_Flux_T_euler_lax_friedrichs_jacobian
#undef compute_Numerical_Flux_T_euler_roe_pike
#undef compute_Numerical_Flux_T_euler_roe_pike_jacobian

#undef compute_Numerical_Flux_T_navier_stokes_central
#undef compute_Numerical_Flux_T_navier_stokes_central_jacobian
///\}

#undef combine_num_flux_boundary_T
#undef combine_num_flux_boundary_dnnf_ds_T
#undef combine_num_flux_boundary_dnnf_dg_g_T
#undef combine_num_flux_boundary_dnnf_dg_s_T
#undef combine_num_flux_boundary_T
#undef Fluxes_LR
#undef constructor_Fluxes_LR
#undef destructor_Fluxes_LR
#undef set_Numerical_Flux_member
#undef set_provided_Numerical_Flux_members
#undef set_provided_Numerical_Flux_jacobian_members
#undef compute_Numerical_Flux_T_central
#undef compute_Numerical_Flux_T_central_jacobian
#undef set_Numerical_Flux_Energy_member
#undef min_abs_T
#undef max_abs_T
