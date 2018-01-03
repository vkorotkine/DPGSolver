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
 *         euler equations.
 */

#include "periodic_vortex/undef_templates_solution_periodic_vortex.h"
#include "supersonic_vortex/undef_templates_solution_supersonic_vortex.h"
#include "free_stream/undef_templates_solution_free_stream.h"

#undef set_function_pointers_solution_euler_T
#undef convert_variables_T
