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

#ifndef DPG__templates_numerical_flux_h__INCLUDED
#define DPG__templates_numerical_flux_h__INCLUDED
/** \file
 *  \brief Provides the macro definitions used for c-style templating related to the real numerical flux functions.
 */

#define Numerical_Flux_Input_T   Numerical_Flux_Input   ///< Numerical_Flux_Input parameter.
#define Numerical_Flux_Input_R   Numerical_Flux_Input   ///< Real Numerical_Flux_Input parameter.
#define mutable_Numerical_Flux_T mutable_Numerical_Flux ///< mutable_Numerical_Flux parameter.

#define compute_Numerical_Flux_T_advection_upwind          compute_Numerical_Flux_advection_upwind          ///< Standard.
#define compute_Numerical_Flux_T_advection_upwind_jacobian compute_Numerical_Flux_advection_upwind_jacobian ///< Standard.

#define compute_Numerical_Flux_T_euler_lax_friedrichs          compute_Numerical_Flux_euler_lax_friedrichs          ///< Standard.
#define compute_Numerical_Flux_T_euler_lax_friedrichs_jacobian compute_Numerical_Flux_euler_lax_friedrichs_jacobian ///< Standard.
#define compute_Numerical_Flux_T_euler_roe_pike                compute_Numerical_Flux_euler_roe_pike                ///< Standard.
#define compute_Numerical_Flux_T_euler_roe_pike_jacobian       compute_Numerical_Flux_euler_roe_pike_jacobian       ///< Standard.

#endif // DPG__templates_numerical_flux_h__INCLUDED
