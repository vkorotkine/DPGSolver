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
 *  \brief Undefine macro definitions for c-style templating relating to flux containers/functions.
 */

///\{ \name Data types
#undef Flux_Input_T
#undef Flux_T
#undef mutable_Flux_T
///\}

///\{ \name Function pointers
#undef compute_Flux_fptr_T
///\}

///\{ \name Function names
#undef constructor_Flux_Input_T
#undef constructor_Flux_Input_T_e
#undef destructor_Flux_Input_T
#undef constructor_Flux_T
#undef destructor_Flux_T
#undef destructor_conditional_Flux_T
#undef compute_Flux_1_T
#undef compute_Flux_2_T
#undef compute_Flux_12_T
#undef increment_pointers_T
///\}

///\{ \name Function names (pde specific)
#undef compute_Flux_T_advection
#undef compute_Flux_T_diffusion
#undef compute_Flux_T_euler
#undef compute_Flux_T_navier_stokes

#undef compute_V2_from_uvw_T
#undef compute_V2_from_rhouvw_T
///\}

#undef Flux_Data_Advection
#undef Flux_Data_Diffusion
#undef Flux_Data_Euler
#undef Flux_Data_Navier_Stokes
#undef Partials_Scalar
#undef Partials_Vector
#undef Partials_Tensor
#undef compute_Flux_Advection_fptr
#undef compute_Flux_diffusion_fptr
#undef compute_Flux_Euler_fptr
#undef compute_Flux_Navier_Stokes_fptr
#undef compute_dmu_ds_fptr
#undef get_compute_Flux_Advection_fptr
#undef compute_Flux_Advection_100
#undef compute_Flux_Advection_110
#undef compute_Flux_advection_0
#undef compute_Flux_advection_1
#undef get_compute_Flux_diffusion_fptr
#undef compute_Flux_diffusion_100
#undef compute_Flux_diffusion_101
#undef compute_Flux_diffusion_0
#undef compute_Flux_diffusion_1
#undef get_compute_Flux_Euler_fptr
#undef compute_Flux_Euler_100
#undef compute_Flux_Euler_110
#undef compute_Flux_Euler_111
#undef compute_Flux_euler_0
#undef compute_Flux_euler_1
#undef compute_Flux_euler_2
#undef get_compute_Flux_Navier_Stokes_fptr
#undef check_if_mu_is_const
#undef get_compute_dmu_ds_fptr
#undef compute_Pr
#undef compute_mu_p
#undef compute_uvw_p
#undef compute_duvw_p
#undef compute_tau_p
#undef compute_dTs_p
#undef compute_Flux_Navier_Stokes_100
#undef compute_Flux_Navier_Stokes_111
#undef compute_dmu_ds_constant
#undef compute_dmu_ds_sutherland
#undef compute_uvw
#undef compute_duvw_ds
#undef compute_duvw
#undef compute_dduvw_ds
#undef compute_dduvw_dg
#undef compute_tau
#undef compute_dtau_ds
#undef compute_dtau_dg
#undef compute_dTs
#undef compute_ddTs_ds
#undef compute_ddTs_dg
#undef compute_Flux_navier_stokes_0
#undef compute_Flux_navier_stokes_1s
#undef compute_Flux_navier_stokes_1g
