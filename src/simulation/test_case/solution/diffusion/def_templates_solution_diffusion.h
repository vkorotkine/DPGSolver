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
 *         diffusion equation.
 */

#include "default_steady/def_templates_solution_diffusion_default_steady.h"

#if TYPE_RC == TYPE_REAL

///\{ \name Function names
#define set_function_pointers_solution_diffusion_T set_function_pointers_solution_diffusion
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function names
#define set_function_pointers_solution_diffusion_T set_function_pointers_solution_diffusion_c
///\}

#endif
