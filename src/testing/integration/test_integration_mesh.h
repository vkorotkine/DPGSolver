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

#ifndef DPG__test_integration_mesh_h__INCLUDED
#define DPG__test_integration_mesh_h__INCLUDED
/**	\file
 *	\brief Provides functionality for integration testing of the mesh processing.
 */

struct Test_Info;

/**	\test Performs integration testing for the mesh processing.
 *
 *	Compares members of the following containers with their expected values:
 *	- \ref Mesh_Data;
 *	- \ref Mesh_Connectivity;
 *	- \ref Mesh_Vertices.
 */
void test_integration_mesh
	(struct Test_Info*const test_info, ///< \ref Test_Info.
	 const char*const mesh_name        ///< The test mesh name.
	);

#endif // DPG__test_integration_mesh_h__INCLUDED
