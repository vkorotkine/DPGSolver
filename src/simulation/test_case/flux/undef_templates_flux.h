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
///\}
