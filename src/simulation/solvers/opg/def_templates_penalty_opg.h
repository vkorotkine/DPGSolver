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
 *  \brief Provides macro definitions for instantiations relating to \ref def_templates_penalty_opg_T.h and
 *         \ref def_templates_penalty_opg_T.c.
 */

#include "def_templates_penalty_opg_advection.h"

#if TYPE_RC == TYPE_REAL

///\{ \name Function names
#define reset_penalty_indicators_opg_T reset_penalty_indicators_opg_d
#define constructor_rlhs_f_test_penalty_unsupported_T constructor_rlhs_f_test_penalty_unsupported_d
#define constructor_rlhs_f_test_penalty_do_nothing_T constructor_rlhs_f_test_penalty_do_nothing_d
#define read_bc_test_s_type_T read_bc_test_s_type_d
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function names
#define reset_penalty_indicators_opg_T reset_penalty_indicators_opg_c
#define constructor_rlhs_f_test_penalty_unsupported_T constructor_rlhs_f_test_penalty_unsupported_c
#define constructor_rlhs_f_test_penalty_do_nothing_T constructor_rlhs_f_test_penalty_do_nothing_c
#define read_bc_test_s_type_T read_bc_test_s_type_c
///\}

#endif
