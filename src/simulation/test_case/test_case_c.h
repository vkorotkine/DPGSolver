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

#ifndef DPG__test_case_c_h__INCLUDED
#define DPG__test_case_c_h__INCLUDED
/** \file
 *  \brief Provides additional `complex` functions for test cases.
 */

struct Simulation;

/// \brief Convert a \ref Test_Case_T from `real` to `complex` or vice versa.
void convert_to_Test_Case_rc
	(struct Simulation* sim, ///< \ref Simulation.
	 const char type_rc_o    ///< The output 'r'eal/'c'omplex type.
	);

#endif // DPG__test_case_c_h__INCLUDED
