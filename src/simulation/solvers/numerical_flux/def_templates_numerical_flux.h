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
#define Numerical_Flux_Input_R   Numerical_Flux_Input
#define mutable_Numerical_Flux_T mutable_Numerical_Flux
#define Numerical_Flux_T         Numerical_Flux
///\}

///\{ \name Function pointers
#define compute_Numerical_Flux_fptr_T compute_Numerical_Flux_fptr
///\}

///\{ \name Function names
#define constructor_Numerical_Flux_Input_T     constructor_Numerical_Flux_Input
#define destructor_Numerical_Flux_Input_T      destructor_Numerical_Flux_Input
#define destructor_Numerical_Flux_Input_data_T destructor_Numerical_Flux_Input_data

#define constructor_Numerical_Flux_T constructor_Numerical_Flux
#define destructor_Numerical_Flux_T  destructor_Numerical_Flux

#define compute_Numerical_Flux_1_T compute_Numerical_Flux_1
#define compute_Numerical_Flux_12_T compute_Numerical_Flux_12
///\}

///\{ \name Function names (pde specific)
#define compute_Numerical_Flux_T_advection_upwind          compute_Numerical_Flux_advection_upwind
#define compute_Numerical_Flux_T_advection_upwind_jacobian compute_Numerical_Flux_advection_upwind_jacobian

#define compute_Numerical_Flux_T_euler_lax_friedrichs          compute_Numerical_Flux_euler_lax_friedrichs
#define compute_Numerical_Flux_T_euler_lax_friedrichs_jacobian compute_Numerical_Flux_euler_lax_friedrichs_jacobian
#define compute_Numerical_Flux_T_euler_roe_pike                compute_Numerical_Flux_euler_roe_pike
#define compute_Numerical_Flux_T_euler_roe_pike_jacobian       compute_Numerical_Flux_euler_roe_pike_jacobian
///\}

#elif TYPE_RC == TYPE_COMPLEX

#endif
