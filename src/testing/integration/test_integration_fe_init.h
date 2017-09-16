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

#ifndef DPG__test_integration_fe_init_h__INCLUDED
#define DPG__test_integration_fe_init_h__INCLUDED
/**	\file
 *	\brief Provides functionality for integration testing of the finite element initialization.
 */

struct Test_Info;

/**	\test Performs integration testing for the finite element initialization.
 *
 *	Compares members of the following containers with their expected values:
 *	- \ref Volume;
 *	- \ref Face.
 */
void test_integration_fe_init
	(struct Test_Info*const test_info, ///< \ref Test_Info.
	 const char*const ctrl_name        ///< The name of the control file.
	);

#endif // DPG__test_integration_fe_init_h__INCLUDED
