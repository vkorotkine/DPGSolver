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

#ifndef DPG__definitions_opg_h__INCLUDED
#define DPG__definitions_opg_h__INCLUDED
/** \file
 *  \brief Provides the definitions relating to the 'o'ptimal 'p'etrov-'g'alerkin method.
 */

///\{ \name Availalble options for \ref OPG_Solver_Face_T::bc_test_s in relation to \ref Solver_Volume_T::test_s_coef.
#define BC_TEST_S_DO_NOTHING 100 ///< Do not impose any constraints on the boundary.
#define BC_TEST_S_SPEC_CONST 101 ///< Set the boundary values to a "spec"ified "const"ant.
///\}

#endif // DPG__definitions_opg_h__INCLUDED
