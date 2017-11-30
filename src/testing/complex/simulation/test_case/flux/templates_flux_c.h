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

#ifndef DPG__templates_flux_c_h__INCLUDED
#define DPG__templates_flux_c_h__INCLUDED
/** \file
 *  \brief Provides the macro definitions used for c-style templating related to the real flux functions.
 */

#define Flux_Input_T   Flux_Input_c   ///< Flux_Input parameter.
#define Flux_Input_R   Flux_Input     ///< Real Flux_Input parameter.
#define mutable_Flux_T mutable_Flux_c ///< mutable_Flux parameter.

#define compute_Flux_T_advection          compute_Flux_c_advection          ///< Standard.
#define compute_Flux_T_advection_jacobian compute_Flux_c_advection_jacobian ///< Standard.

#define compute_Flux_T_euler          compute_Flux_c_euler          ///< Standard.
#define compute_Flux_T_euler_jacobian compute_Flux_c_euler_jacobian ///< Standard.

#endif // DPG__templates_flux_c_h__INCLUDED
