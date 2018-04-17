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
 *  \brief Provides the macro definitions used for c-style templating related to the boundary
 *         containers/functions.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Data types
#define Boundary_Value_Input_T   Boundary_Value_Input
#define Boundary_Value_T         Boundary_Value
#define mutable_Boundary_Value_T mutable_Boundary_Value
///\}

///\{ \name Function pointers
#define constructor_Boundary_Value_Input_face_fptr_T constructor_Boundary_Value_Input_face_fptr
#define constructor_Boundary_Value_fptr_T            constructor_Boundary_Value_fptr

#define constructor_s_fc_interp_T constructor_s_fc_interp
///\}

///\{ \name Function names
#define constructor_Boundary_Value_Input_face_s_fcl_interp_T  constructor_Boundary_Value_Input_face_s_fcl_interp
#define constructor_Boundary_Value_Input_face_sg_fcl_interp_T constructor_Boundary_Value_Input_face_sg_fcl_interp
#define destructor_Boundary_Value_Input_T                     destructor_Boundary_Value_Input
#define constructor_Boundary_Value_s_fcl_interp_T             constructor_Boundary_Value_s_fcl_interp
#define destructor_Boundary_Value_T                           destructor_Boundary_Value
#define constructor_Boundary_Value_T_grad_from_internal       constructor_Boundary_Value_grad_from_internal
///\}

///\{ \name Function names (pde specific)
#define constructor_Boundary_Value_T_advection_inflow   constructor_Boundary_Value_advection_inflow
#define constructor_Boundary_Value_T_advection_outflow  constructor_Boundary_Value_advection_outflow
#define constructor_Boundary_Value_T_advection_slipwall constructor_Boundary_Value_advection_slipwall
#define constructor_Boundary_Value_T_advection_upwind   constructor_Boundary_Value_advection_upwind

#define constructor_Boundary_Value_T_diffusion_dirichlet constructor_Boundary_Value_diffusion_dirichlet
#define constructor_Boundary_Value_T_diffusion_neumann   constructor_Boundary_Value_diffusion_neumann

#define constructor_Boundary_Value_T_euler_riemann            constructor_Boundary_Value_euler_riemann
#define constructor_Boundary_Value_T_euler_slipwall           constructor_Boundary_Value_euler_slipwall
#define constructor_Boundary_Value_T_euler_supersonic_inflow  constructor_Boundary_Value_euler_supersonic_inflow
#define constructor_Boundary_Value_T_euler_supersonic_outflow constructor_Boundary_Value_euler_supersonic_outflow
#define constructor_Boundary_Value_T_euler_back_pressure      constructor_Boundary_Value_euler_back_pressure
#define constructor_Boundary_Value_T_euler_total_tp           constructor_Boundary_Value_euler_total_tp

#define constructor_Boundary_Value_T_navier_stokes_no_slip_all_rotating   constructor_Boundary_Value_navier_stokes_no_slip_all_rotating
#define constructor_Boundary_Value_T_navier_stokes_no_slip_flux_adiabatic constructor_Boundary_Value_navier_stokes_no_slip_flux_adiabatic
#define constructor_Boundary_Value_T_navier_stokes_no_slip_flux_diabatic  constructor_Boundary_Value_navier_stokes_no_slip_flux_diabatic
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Data types
#define Boundary_Value_Input_T   Boundary_Value_Input_c
#define Boundary_Value_T         Boundary_Value_c
#define mutable_Boundary_Value_T mutable_Boundary_Value_c
///\}

///\{ \name Function pointers
#define constructor_Boundary_Value_Input_face_fptr_T constructor_Boundary_Value_Input_c_face_fptr
#define constructor_Boundary_Value_fptr_T            constructor_Boundary_Value_c_fptr

#define constructor_s_fc_interp_T constructor_s_fc_interp_c
///\}

///\{ \name Function names
#define constructor_Boundary_Value_Input_face_s_fcl_interp_T  constructor_Boundary_Value_Input_c_face_s_fcl_interp
#define constructor_Boundary_Value_Input_face_sg_fcl_interp_T constructor_Boundary_Value_Input_c_face_sg_fcl_interp
#define destructor_Boundary_Value_Input_T                     destructor_Boundary_Value_Input_c
#define constructor_Boundary_Value_s_fcl_interp_T             constructor_Boundary_Value_c_s_fcl_interp
#define destructor_Boundary_Value_T                           destructor_Boundary_Value_c
#define constructor_Boundary_Value_T_grad_from_internal       constructor_Boundary_Value_c_grad_from_internal
///\}

///\{ \name Function names
#define constructor_Boundary_Value_T_advection_inflow   constructor_Boundary_Value_c_advection_inflow
#define constructor_Boundary_Value_T_advection_outflow  constructor_Boundary_Value_c_advection_outflow
#define constructor_Boundary_Value_T_advection_slipwall constructor_Boundary_Value_c_advection_slipwall
#define constructor_Boundary_Value_T_advection_upwind   constructor_Boundary_Value_c_advection_upwind

#define constructor_Boundary_Value_T_diffusion_dirichlet constructor_Boundary_Value_c_diffusion_dirichlet
#define constructor_Boundary_Value_T_diffusion_neumann   constructor_Boundary_Value_c_diffusion_neumann

#define constructor_Boundary_Value_T_euler_riemann            constructor_Boundary_Value_c_euler_riemann
#define constructor_Boundary_Value_T_euler_slipwall           constructor_Boundary_Value_c_euler_slipwall
#define constructor_Boundary_Value_T_euler_supersonic_inflow  constructor_Boundary_Value_c_euler_supersonic_inflow
#define constructor_Boundary_Value_T_euler_supersonic_outflow constructor_Boundary_Value_c_euler_supersonic_outflow
#define constructor_Boundary_Value_T_euler_back_pressure      constructor_Boundary_Value_c_euler_back_pressure
#define constructor_Boundary_Value_T_euler_total_tp           constructor_Boundary_Value_c_euler_total_tp

#define constructor_Boundary_Value_T_navier_stokes_no_slip_all_rotating   constructor_Boundary_Value_c_navier_stokes_no_slip_all_rotating
#define constructor_Boundary_Value_T_navier_stokes_no_slip_flux_adiabatic constructor_Boundary_Value_c_navier_stokes_no_slip_flux_adiabatic
#define constructor_Boundary_Value_T_navier_stokes_no_slip_flux_diabatic  constructor_Boundary_Value_c_navier_stokes_no_slip_flux_diabatic
///\}

#endif
