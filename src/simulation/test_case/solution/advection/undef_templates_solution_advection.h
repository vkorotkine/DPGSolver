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
 *  \brief Undefine macro definitions for c-style templated containers/functions relating to solution for the
 *         linear advection equation.
 */

#include "undef_templates_solution_advection_default.h"
#include "peterson/undef_templates_solution_advection_peterson.h"
#include "free_stream_advection/undef_templates_solution_free_stream_advection.h"
#include "vortex_advection/undef_templates_solution_vortex_advection.h"
#include "hyperbolic_tan/undef_templates_solution_hyperbolic_tan.h"

#undef compute_b_adv_fptr_T
#undef Sol_Data__Advection_T
#undef set_function_pointers_solution_advection_T
#undef get_sol_data_advection_T
#undef read_data_advection_T
#undef compute_b_adv_constant_T
#undef compute_b_adv_vortex_T

#undef set_function_pointers_num_flux_T
