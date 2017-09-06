// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_integration_geometry_h__INCLUDED
#define DPG__test_integration_geometry_h__INCLUDED
/**	\file
 *	\brief Provides functionality for integration testing of the finite element initialization.
 */

#include "test_base.h"

/**	\test Performs integration testing for the geometry initialization.
 *
 *	\todo Update the description.
 */
void test_integration_geometry
	(struct Test_Info*const test_info, ///< \ref Test_Info.
	 const char*const ctrl_name        ///< The name of the control file.
	);

#endif // DPG__test_integration_geometry_h__INCLUDED
