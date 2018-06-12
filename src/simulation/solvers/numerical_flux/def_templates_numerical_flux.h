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
 *  \brief Provides the macro definitions used for c-style templating related to the numerical flux functions.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Data types
#define Numerical_Flux_Input_T   Numerical_Flux_Input
#define mutable_Numerical_Flux_T mutable_Numerical_Flux
#define Numerical_Flux_T         Numerical_Flux

#define Neigh_Info_NF_T   Neigh_Info_NF
#define m_Neigh_Info_NF_T m_Neigh_Info_NF
///\}

///\{ \name Function pointers
#define compute_Numerical_Flux_fptr_T compute_Numerical_Flux_fptr
///\}

///\{ \name Function names
#define constructor_Numerical_Flux_Input_T     constructor_Numerical_Flux_Input
#define destructor_Numerical_Flux_Input_T      destructor_Numerical_Flux_Input

#define constructor_Numerical_Flux_T constructor_Numerical_Flux
#define destructor_Numerical_Flux_T  destructor_Numerical_Flux

#define compute_Numerical_Flux_1_T  compute_Numerical_Flux_1
#define compute_Numerical_Flux_2_T  compute_Numerical_Flux_2
#define compute_Numerical_Flux_12_T compute_Numerical_Flux_12
///\}

///\{ \name Function names (pde specific)
#define compute_Numerical_Flux_T_advection_upwind          compute_Numerical_Flux_advection_upwind
#define compute_Numerical_Flux_T_advection_upwind_jacobian compute_Numerical_Flux_advection_upwind_jacobian

#define compute_Numerical_Flux_T_diffusion_central          compute_Numerical_Flux_diffusion_central
#define compute_Numerical_Flux_T_diffusion_central_jacobian compute_Numerical_Flux_diffusion_central_jacobian

#define compute_Numerical_Flux_T_euler_lax_friedrichs          compute_Numerical_Flux_euler_lax_friedrichs
#define compute_Numerical_Flux_T_euler_lax_friedrichs_jacobian compute_Numerical_Flux_euler_lax_friedrichs_jacobian
#define compute_Numerical_Flux_T_euler_roe_pike                compute_Numerical_Flux_euler_roe_pike
#define compute_Numerical_Flux_T_euler_roe_pike_jacobian       compute_Numerical_Flux_euler_roe_pike_jacobian

#define compute_Numerical_Flux_T_navier_stokes_central          compute_Numerical_Flux_navier_stokes_central
#define compute_Numerical_Flux_T_navier_stokes_central_jacobian compute_Numerical_Flux_navier_stokes_central_jacobian
///\}

#define combine_num_flux_boundary_T combine_num_flux_boundary_T
#define combine_num_flux_boundary_dnnf_ds_T combine_num_flux_boundary_dnnf_ds_T
#define combine_num_flux_boundary_dnnf_dg_g_T combine_num_flux_boundary_dnnf_dg_g_T
#define combine_num_flux_boundary_dnnf_dg_s_T combine_num_flux_boundary_dnnf_dg_s_T
#define combine_num_flux_boundary_T combine_num_flux_boundary_T
#define Fluxes_LR Fluxes_LR
#define constructor_Fluxes_LR constructor_Fluxes_LR
#define destructor_Fluxes_LR destructor_Fluxes_LR
#define set_Numerical_Flux_member set_Numerical_Flux_member
#define set_provided_Numerical_Flux_members set_provided_Numerical_Flux_members
#define set_provided_Numerical_Flux_jacobian_members set_provided_Numerical_Flux_jacobian_members
#define compute_Numerical_Flux_T_central compute_Numerical_Flux_T_central
#define compute_Numerical_Flux_T_central_jacobian compute_Numerical_Flux_T_central_jacobian
#define set_Numerical_Flux_Energy_member set_Numerical_Flux_Energy_member
#define min_abs_real_T min_abs_real_T
#define max_abs_real_T max_abs_real_T

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Data types
#define Numerical_Flux_Input_T   Numerical_Flux_Input_c
#define mutable_Numerical_Flux_T mutable_Numerical_Flux_c
#define Numerical_Flux_T         Numerical_Flux_c

#define Neigh_Info_NF_T   Neigh_Info_NF_c
#define m_Neigh_Info_NF_T m_Neigh_Info_NF_c
///\}

///\{ \name Function pointers
#define compute_Numerical_Flux_fptr_T compute_Numerical_Flux_fptr_c
///\}

///\{ \name Function names
#define constructor_Numerical_Flux_Input_T     constructor_Numerical_Flux_Input_c
#define destructor_Numerical_Flux_Input_T      destructor_Numerical_Flux_Input_c

#define constructor_Numerical_Flux_T constructor_Numerical_Flux_c
#define destructor_Numerical_Flux_T  destructor_Numerical_Flux_c

#define compute_Numerical_Flux_1_T  compute_Numerical_Flux_1_c
#define compute_Numerical_Flux_2_T  compute_Numerical_Flux_2_c
#define compute_Numerical_Flux_12_T compute_Numerical_Flux_12_c
///\}

///\{ \name Function names (pde specific)
#define compute_Numerical_Flux_T_advection_upwind          compute_Numerical_Flux_c_advection_upwind
#define compute_Numerical_Flux_T_advection_upwind_jacobian compute_Numerical_Flux_c_advection_upwind_jacobian

#define compute_Numerical_Flux_T_diffusion_central          compute_Numerical_Flux_c_diffusion_central
#define compute_Numerical_Flux_T_diffusion_central_jacobian compute_Numerical_Flux_c_diffusion_central_jacobian

#define compute_Numerical_Flux_T_euler_lax_friedrichs          compute_Numerical_Flux_c_euler_lax_friedrichs
#define compute_Numerical_Flux_T_euler_lax_friedrichs_jacobian compute_Numerical_Flux_c_euler_lax_friedrichs_jacobian
#define compute_Numerical_Flux_T_euler_roe_pike                compute_Numerical_Flux_c_euler_roe_pike
#define compute_Numerical_Flux_T_euler_roe_pike_jacobian       compute_Numerical_Flux_c_euler_roe_pike_jacobian

#define compute_Numerical_Flux_T_navier_stokes_central          compute_Numerical_Flux_c_navier_stokes_central
#define compute_Numerical_Flux_T_navier_stokes_central_jacobian compute_Numerical_Flux_c_navier_stokes_central_jacobian
///\}

#define combine_num_flux_boundary_T combine_num_flux_boundary_T_c
#define combine_num_flux_boundary_dnnf_ds_T combine_num_flux_boundary_dnnf_ds_T_c
#define combine_num_flux_boundary_dnnf_dg_g_T combine_num_flux_boundary_dnnf_dg_g_T_c
#define combine_num_flux_boundary_dnnf_dg_s_T combine_num_flux_boundary_dnnf_dg_s_T_c
#define combine_num_flux_boundary_T combine_num_flux_boundary_T_c
#define Fluxes_LR Fluxes_LR_c
#define constructor_Fluxes_LR constructor_Fluxes_LR_c
#define destructor_Fluxes_LR destructor_Fluxes_LR_c
#define set_Numerical_Flux_member set_Numerical_Flux_member_c
#define set_provided_Numerical_Flux_members set_provided_Numerical_Flux_members_c
#define set_provided_Numerical_Flux_jacobian_members set_provided_Numerical_Flux_jacobian_members_c
#define compute_Numerical_Flux_T_central compute_Numerical_Flux_T_central_c
#define compute_Numerical_Flux_T_central_jacobian compute_Numerical_Flux_T_central_jacobian_c
#define set_Numerical_Flux_Energy_member set_Numerical_Flux_Energy_member_c
#define min_abs_real_T min_abs_real_T_c
#define max_abs_real_T max_abs_real_T_c

#endif
