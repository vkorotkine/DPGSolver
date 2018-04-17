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
 *  \brief Undefine macro definitions for c-style templating relating to boundary containers/functions.
 */

///\{ \name Data types
#undef Boundary_Value_Input_T
#undef Boundary_Value_T
#undef mutable_Boundary_Value_T
///\}

///\{ \name Function pointers
#undef constructor_Boundary_Value_Input_face_fptr_T
#undef constructor_Boundary_Value_fptr_T

#undef constructor_s_fc_interp_T
///\}

///\{ \name Function names (general)
#undef constructor_Boundary_Value_Input_face_s_fcl_interp_T
#undef constructor_Boundary_Value_Input_face_sg_fcl_interp_T
#undef destructor_Boundary_Value_Input_T
#undef constructor_Boundary_Value_s_fcl_interp_T
#undef destructor_Boundary_Value_T
#undef constructor_Boundary_Value_T_grad_from_internal
///\}

///\{ \name Function names (pde specific)
#undef constructor_Boundary_Value_T_advection_inflow
#undef constructor_Boundary_Value_T_advection_outflow
#undef constructor_Boundary_Value_T_advection_slipwall
#undef constructor_Boundary_Value_T_advection_upwind

#undef constructor_Boundary_Value_T_diffusion_dirichlet
#undef constructor_Boundary_Value_T_diffusion_neumann

#undef constructor_Boundary_Value_T_euler_riemann
#undef constructor_Boundary_Value_T_euler_slipwall
#undef constructor_Boundary_Value_T_euler_supersonic_inflow
#undef constructor_Boundary_Value_T_euler_supersonic_outflow
#undef constructor_Boundary_Value_T_euler_back_pressure
#undef constructor_Boundary_Value_T_euler_total_tp

#undef constructor_Boundary_Value_T_navier_stokes_no_slip_all_rotating
#undef constructor_Boundary_Value_T_navier_stokes_no_slip_flux_adiabatic
#undef constructor_Boundary_Value_T_navier_stokes_no_slip_flux_diabatic
///\}
