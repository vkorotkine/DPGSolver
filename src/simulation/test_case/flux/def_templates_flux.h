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
#define destructor_Flux_Input_T  destructor_Flux_Input
#define constructor_Flux_T       constructor_Flux
#define destructor_Flux_T        destructor_Flux
#define compute_Flux_1_T         compute_Flux_1
#define compute_Flux_12_T        compute_Flux_12
///\}

///\{ \name Function names (pde_specific)
#define compute_Flux_T_advection          compute_Flux_advection
#define compute_Flux_T_advection_jacobian compute_Flux_advection_jacobian

#define compute_Flux_T_euler          compute_Flux_euler
#define compute_Flux_T_euler_jacobian compute_Flux_euler_jacobian
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
#define destructor_Flux_Input_T  destructor_Flux_Input_c
#define constructor_Flux_T       constructor_Flux_c
#define destructor_Flux_T        destructor_Flux_c
#define compute_Flux_1_T         compute_Flux_c_1
#define compute_Flux_12_T        compute_Flux_c_12
///\}

///\{ \name Function names (pde_specific)
#define compute_Flux_T_advection          compute_Flux_c_advection
#define compute_Flux_T_advection_jacobian compute_Flux_c_advection_jacobian

#define compute_Flux_T_euler          compute_Flux_c_euler
#define compute_Flux_T_euler_jacobian compute_Flux_c_euler_jacobian
///\}

#endif
