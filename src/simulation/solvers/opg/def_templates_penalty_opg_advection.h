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
 *  \brief Provides macro definitions for instantiations relating to \ref def_templates_penalty_opg_advection_T.h and
 *         \ref def_templates_penalty_opg_advection_T.c.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Function names
#define constructor_rlhs_f_test_penalty_advection_upwind_T constructor_rlhs_f_test_penalty_advection_upwind_d
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Function names
#define constructor_rlhs_f_test_penalty_advection_upwind_T constructor_rlhs_f_test_penalty_advection_upwind_c
///\}

#endif
