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
 *  \brief Provides the macro definitions used for c-style templating related to the real flux functions.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Data types
#define Flux_Input_T   Flux_Input
#define Flux_T         Flux
#define mutable_Flux_T mutable_Flux
///\}

///\{ \name Function pointers
#define compute_Flux_fptr_T compute_Flux_fptr
///\}

///\{ \name Function names
#define constructor_Flux_Input_T constructor_Flux_Input
#define constructor_Flux_Input_T_e constructor_Flux_Input_e
#define destructor_Flux_Input_T  destructor_Flux_Input
#define constructor_Flux_T       constructor_Flux
#define destructor_Flux_T        destructor_Flux
#define destructor_conditional_Flux_T destructor_conditional_Flux
#define compute_Flux_1_T         compute_Flux_1
#define compute_Flux_2_T         compute_Flux_2
#define compute_Flux_12_T        compute_Flux_12
#define increment_pointers_T     increment_pointers
///\}

///\{ \name Function names (pde_specific)
#define compute_Flux_T_advection     compute_Flux_advection
#define compute_Flux_T_diffusion     compute_Flux_diffusion
#define compute_Flux_T_euler         compute_Flux_euler
#define compute_Flux_T_navier_stokes compute_Flux_navier_stokes

#define compute_V2_from_uvw_T    compute_V2_from_uvw
#define compute_V2_from_rhouvw_T compute_V2_from_rhouvw
///\}

///\{ \name Static names
#define Flux_Data_Advection Flux_Data_Advection
#define Flux_Data_Diffusion Flux_Data_Diffusion
#define Flux_Data_Euler Flux_Data_Euler
#define Flux_Data_Navier_Stokes Flux_Data_Navier_Stokes
#define Partials_Scalar Partials_Scalar
#define Partials_Vector Partials_Vector
#define Partials_Tensor Partials_Tensor
#define compute_Flux_Advection_fptr compute_Flux_Advection_fptr
#define compute_Flux_diffusion_fptr compute_Flux_diffusion_fptr
#define compute_Flux_Euler_fptr compute_Flux_Euler_fptr
#define compute_Flux_Navier_Stokes_fptr compute_Flux_Navier_Stokes_fptr
#define compute_dmu_ds_fptr compute_dmu_ds_fptr
#define get_compute_Flux_Advection_fptr get_compute_Flux_Advection_fptr
#define compute_Flux_Advection_100 compute_Flux_Advection_100
#define compute_Flux_Advection_110 compute_Flux_Advection_110
#define compute_Flux_advection_0 compute_Flux_advection_0
#define compute_Flux_advection_1 compute_Flux_advection_1
#define get_compute_Flux_diffusion_fptr get_compute_Flux_diffusion_fptr
#define compute_Flux_diffusion_100 compute_Flux_diffusion_100
#define compute_Flux_diffusion_101 compute_Flux_diffusion_101
#define compute_Flux_diffusion_0 compute_Flux_diffusion_0
#define compute_Flux_diffusion_1 compute_Flux_diffusion_1
#define get_compute_Flux_Euler_fptr get_compute_Flux_Euler_fptr
#define compute_Flux_Euler_100 compute_Flux_Euler_100
#define compute_Flux_Euler_110 compute_Flux_Euler_110
#define compute_Flux_Euler_111 compute_Flux_Euler_111
#define compute_Flux_euler_0 compute_Flux_euler_0
#define compute_Flux_euler_1 compute_Flux_euler_1
#define compute_Flux_euler_2 compute_Flux_euler_2
#define get_compute_Flux_Navier_Stokes_fptr get_compute_Flux_Navier_Stokes_fptr
#define check_if_mu_is_const check_if_mu_is_const
#define get_compute_dmu_ds_fptr get_compute_dmu_ds_fptr
#define compute_Pr compute_Pr
#define compute_mu_p compute_mu_p
#define compute_uvw_p compute_uvw_p
#define compute_duvw_p compute_duvw_p
#define compute_tau_p compute_tau_p
#define compute_dTs_p compute_dTs_p
#define compute_Flux_Navier_Stokes_100 compute_Flux_Navier_Stokes_100
#define compute_Flux_Navier_Stokes_111 compute_Flux_Navier_Stokes_111
#define compute_dmu_ds_constant compute_dmu_ds_constant
#define compute_dmu_ds_sutherland compute_dmu_ds_sutherland
#define compute_uvw compute_uvw
#define compute_duvw_ds compute_duvw_ds
#define compute_duvw compute_duvw
#define compute_dduvw_ds compute_dduvw_ds
#define compute_dduvw_dg compute_dduvw_dg
#define compute_tau compute_tau
#define compute_dtau_ds compute_dtau_ds
#define compute_dtau_dg compute_dtau_dg
#define compute_dTs compute_dTs
#define compute_ddTs_ds compute_ddTs_ds
#define compute_ddTs_dg compute_ddTs_dg
#define compute_Flux_navier_stokes_0 compute_Flux_navier_stokes_0
#define compute_Flux_navier_stokes_1s compute_Flux_navier_stokes_1s
#define compute_Flux_navier_stokes_1g compute_Flux_navier_stokes_1g
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Data types
#define Flux_Input_T   Flux_Input_c
#define Flux_T         Flux_c
#define mutable_Flux_T mutable_Flux_c
///\}

///\{ \name Function pointers
#define compute_Flux_fptr_T compute_Flux_c_fptr
///\}

///\{ \name Function names
#define constructor_Flux_Input_T constructor_Flux_Input_c
#define constructor_Flux_Input_T_e constructor_Flux_Input_c_e
#define destructor_Flux_Input_T  destructor_Flux_Input_c
#define constructor_Flux_T       constructor_Flux_c
#define destructor_Flux_T        destructor_Flux_c
#define destructor_conditional_Flux_T destructor_conditional_Flux_c
#define compute_Flux_1_T         compute_Flux_c_1
#define compute_Flux_2_T         compute_Flux_c_2
#define compute_Flux_12_T        compute_Flux_c_12
#define increment_pointers_T     increment_pointers_c
///\}

///\{ \name Function names (pde_specific)
#define compute_Flux_T_advection     compute_Flux_c_advection
#define compute_Flux_T_diffusion     compute_Flux_c_diffusion
#define compute_Flux_T_euler         compute_Flux_c_euler
#define compute_Flux_T_navier_stokes compute_Flux_c_navier_stokes

#define compute_V2_from_uvw_T    compute_V2_from_uvw_c
#define compute_V2_from_rhouvw_T compute_V2_from_rhouvw_c
///\}

///\{ \name Static names
#define Flux_Data_Advection Flux_Data_Advection_c
#define Flux_Data_Diffusion Flux_Data_Diffusion_c
#define Flux_Data_Euler Flux_Data_Euler_c
#define Flux_Data_Navier_Stokes Flux_Data_Navier_Stokes_c
#define Partials_Scalar Partials_Scalar_c
#define Partials_Vector Partials_Vector_c
#define Partials_Tensor Partials_Tensor_c
#define compute_Flux_Advection_fptr compute_Flux_Advection_fptr_c
#define compute_Flux_diffusion_fptr compute_Flux_diffusion_fptr_c
#define compute_Flux_Euler_fptr compute_Flux_Euler_fptr_c
#define compute_Flux_Navier_Stokes_fptr compute_Flux_Navier_Stokes_fptr_c
#define compute_dmu_ds_fptr compute_dmu_ds_fptr_c
#define get_compute_Flux_Advection_fptr get_compute_Flux_Advection_fptr_c
#define compute_Flux_Advection_100 compute_Flux_Advection_100_c
#define compute_Flux_Advection_110 compute_Flux_Advection_110_c
#define compute_Flux_advection_0 compute_Flux_advection_0_c
#define compute_Flux_advection_1 compute_Flux_advection_1_c
#define get_compute_Flux_diffusion_fptr get_compute_Flux_diffusion_fptr_c
#define compute_Flux_diffusion_100 compute_Flux_diffusion_100_c
#define compute_Flux_diffusion_101 compute_Flux_diffusion_101_c
#define compute_Flux_diffusion_0 compute_Flux_diffusion_0_c
#define compute_Flux_diffusion_1 compute_Flux_diffusion_1_c
#define get_compute_Flux_Euler_fptr get_compute_Flux_Euler_fptr_c
#define compute_Flux_Euler_100 compute_Flux_Euler_100_c
#define compute_Flux_Euler_110 compute_Flux_Euler_110_c
#define compute_Flux_Euler_111 compute_Flux_Euler_111_c
#define compute_Flux_euler_0 compute_Flux_euler_0_c
#define compute_Flux_euler_1 compute_Flux_euler_1_c
#define compute_Flux_euler_2 compute_Flux_euler_2_c
#define get_compute_Flux_Navier_Stokes_fptr get_compute_Flux_Navier_Stokes_fptr_c
#define check_if_mu_is_const check_if_mu_is_const_c
#define get_compute_dmu_ds_fptr get_compute_dmu_ds_fptr_c
#define compute_Pr compute_Pr_c
#define compute_mu_p compute_mu_p_c
#define compute_uvw_p compute_uvw_p_c
#define compute_duvw_p compute_duvw_p_c
#define compute_tau_p compute_tau_p_c
#define compute_dTs_p compute_dTs_p_c
#define compute_Flux_Navier_Stokes_100 compute_Flux_Navier_Stokes_100_c
#define compute_Flux_Navier_Stokes_111 compute_Flux_Navier_Stokes_111_c
#define compute_dmu_ds_constant compute_dmu_ds_constant_c
#define compute_dmu_ds_sutherland compute_dmu_ds_sutherland_c
#define compute_uvw compute_uvw_c
#define compute_duvw_ds compute_duvw_ds_c
#define compute_duvw compute_duvw_c
#define compute_dduvw_ds compute_dduvw_ds_c
#define compute_dduvw_dg compute_dduvw_dg_c
#define compute_tau compute_tau_c
#define compute_dtau_ds compute_dtau_ds_c
#define compute_dtau_dg compute_dtau_dg_c
#define compute_dTs compute_dTs_c
#define compute_ddTs_ds compute_ddTs_ds_c
#define compute_ddTs_dg compute_ddTs_dg_c
#define compute_Flux_navier_stokes_0 compute_Flux_navier_stokes_0_c
#define compute_Flux_navier_stokes_1s compute_Flux_navier_stokes_1s_c
#define compute_Flux_navier_stokes_1g compute_Flux_navier_stokes_1g_c
///\}

#endif
