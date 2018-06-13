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
 *  \brief Provides the macro definitions used for c-style templating related to the solution functions for the
 *         navier-stokes equations.
 */

#include "taylor_couette/def_templates_solution_taylor_couette.h"

#if TYPE_RC == TYPE_REAL

///\{ \name Function pointers
#define compute_mu_fptr_T compute_mu_fptr_d
///\}

///\{ \name Function names
#define set_function_pointers_solution_navier_stokes_T set_function_pointers_solution_navier_stokes
#define convert_variables_gradients_T                  convert_variables_gradients
#define get_compute_mu_fptr_T                          get_compute_mu_fptr_d
#define set_viscosity_type_T                           set_viscosity_type_d
#define compute_mu_constant_T                          compute_mu_constant_d
#define compute_mu_sutherland_T                        compute_mu_sutherland_d
///\}

///\{ \name Static names
#define set_function_pointers_num_flux set_function_pointers_num_flux
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function pointers
#define compute_mu_fptr_T compute_mu_fptr_c
///\}

///\{ \name Function names
#define set_function_pointers_solution_navier_stokes_T set_function_pointers_solution_navier_stokes_c
#define convert_variables_gradients_T                  convert_variables_gradients_c
#define get_compute_mu_fptr_T                          get_compute_mu_fptr_c
#define set_viscosity_type_T                           set_viscosity_type_c
#define compute_mu_constant_T                          compute_mu_constant_c
#define compute_mu_sutherland_T                        compute_mu_sutherland_c
///\}

///\{ \name Static names
#define set_function_pointers_num_flux set_function_pointers_num_flux_c
///\}

#endif
