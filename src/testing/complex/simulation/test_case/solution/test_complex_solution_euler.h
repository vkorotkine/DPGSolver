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

#ifndef DPG__solution_euler_c_h__INCLUDED
#define DPG__solution_euler_c_h__INCLUDED
/** \file
 *  \brief Provides functions relating to the Euler solutions.
 */

struct Test_Case;
struct Simulation;

struct Multiarray_c;

/// \brief `complex` version of \ref convert_variables.
void convert_variables_c
	(struct Multiarray_c* vars, ///< See brief.
	 const char type_i,         ///< See brief.
	 const char type_o          ///< See brief.
	);

#endif // DPG__solution_euler_c_h__INCLUDED
