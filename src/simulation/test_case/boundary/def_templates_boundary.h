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
 *  \brief Provides the macro definitions used for c-style templating related to the real boundary functions.
 */

///\{ \name Data types
#define Boundary_Value_T Boundary_Value
#define Boundary_Value_Input_T Boundary_Value_Input
#define Boundary_Value_Input_R Boundary_Value_Input
///\}

///\{ \name Function names
#define constructor_Boundary_Value_T_advection_inflow constructor_Boundary_Value_advection_inflow
#define constructor_Boundary_Value_T_advection_outflow constructor_Boundary_Value_advection_outflow

#define constructor_Boundary_Value_T_euler_riemann constructor_Boundary_Value_euler_riemann
#define constructor_Boundary_Value_T_euler_slipwall constructor_Boundary_Value_euler_slipwall
///\}
