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

#ifndef DPG__test_unit_nodes_h__INCLUDED
#define DPG__test_unit_nodes_h__INCLUDED
/** \file
 *  \brief Provides functionality for unit testing of the reference coordinates (and associated cubature if relevant).
 */

struct Test_Info;

/// \test Performs unit testing for the nodes.
void test_unit_nodes
	(struct Test_Info*const test_info ///< \ref Test_Info.
	);

#endif // DPG__test_unit_nodes_h__INCLUDED
