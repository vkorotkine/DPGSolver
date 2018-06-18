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
 *  \brief Provides the macro definitions used for c-style templating related to the \ref OPG_Solver_Face_T
 *         containers/functions.
 */

#include "def_templates_face_solver.h"

#if TYPE_RC == TYPE_REAL

///\{ \name Data types
#define OPG_Solver_Face_T OPG_Solver_Face
///\}

///\{ \name Function pointers
#define constructor_rlhs_f_b_test_penalty_T constructor_rlhs_f_b_test_penalty_d
///\}

///\{ \name Function names
#define constructor_derived_OPG_Solver_Face_T constructor_derived_OPG_Solver_Face
#define destructor_derived_OPG_Solver_Face_T  destructor_derived_OPG_Solver_Face

#define get_operator__cv0_vt_fc_T get_operator__cv0_vt_fc_d
#define get_operator__cv1_vt_fc_T get_operator__cv1_vt_fc_d
///\}

///\{ \name Static names
#define constructor_inverse_mass_face_T constructor_inverse_mass_face_d
#define constructor_mass_face_T constructor_mass_face_d
#define set_function_pointers_penalty_T set_function_pointers_penalty_d
#define set_function_pointers_penalty_boundary_T set_function_pointers_penalty_boundary_d
#define set_function_pointers_penalty_boundary_advection_T set_function_pointers_penalty_boundary_advection_d
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Data types
#define OPG_Solver_Face_T OPG_Solver_Face_c
///\}

///\{ \name Function pointers
#define constructor_rlhs_f_b_test_penalty_T constructor_rlhs_f_b_test_penalty_c
///\}

///\{ \name Function names
#define constructor_derived_OPG_Solver_Face_T constructor_derived_OPG_Solver_Face_c
#define destructor_derived_OPG_Solver_Face_T  destructor_derived_OPG_Solver_Face_c

#define get_operator__cv0_vt_fc_T get_operator__cv0_vt_fc_c
#define get_operator__cv1_vt_fc_T get_operator__cv1_vt_fc_c
///\}

///\{ \name Static names
#define constructor_inverse_mass_face_T constructor_inverse_mass_face_T_c
#define constructor_mass_face_T constructor_mass_face_T_c
#define set_function_pointers_penalty_T set_function_pointers_penalty_c
#define set_function_pointers_penalty_boundary_T set_function_pointers_penalty_boundary_c
#define set_function_pointers_penalty_boundary_advection_T set_function_pointers_penalty_boundary_advection_c
///\}

#endif
