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
 *         linear advection equation.
 */

#include "def_templates_solution_advection_default.h"
#include "peterson/def_templates_solution_advection_peterson.h"
#include "free_stream_advection/def_templates_solution_free_stream_advection.h"
#include "vortex_advection/def_templates_solution_vortex_advection.h"
#include "hyperbolic_tan/def_templates_solution_hyperbolic_tan.h"

#if TYPE_RC == TYPE_REAL

///\{ \name Function names
#define set_function_pointers_solution_advection_T set_function_pointers_solution_advection
///\}

#define set_function_pointers_num_flux_T set_function_pointers_num_flux_T

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function names
#define set_function_pointers_solution_advection_T set_function_pointers_solution_advection_c
///\}

#define set_function_pointers_num_flux_T set_function_pointers_num_flux_T_c

#endif
