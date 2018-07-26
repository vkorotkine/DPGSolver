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
 *  \brief Undefine macro definitions for c-style templating relating to \ref OPG_Solver_Face_T containers/functions.
 */

#include "undef_templates_face_solver.h"

///\{ \name Data types
#undef OPG_Solver_Face_T
///\}

#undef constructor_rlhs_f_b_test_penalty_T

///\{ \name Function names
#undef constructor_derived_OPG_Solver_Face_T
#undef destructor_derived_OPG_Solver_Face_T

#undef get_operator__cv1_vt_fc_T
///\}

#undef constructor_inverse_mass_face_T
#undef set_function_pointers_penalty_T
#undef set_function_pointers_penalty_boundary_T
#undef set_function_pointers_penalty_boundary_advection_T
