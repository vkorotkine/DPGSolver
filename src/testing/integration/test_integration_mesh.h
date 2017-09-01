// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_integration_mesh_h__INCLUDED
#define DPG__test_integration_mesh_h__INCLUDED
/**	\file
 *	Provides functionality for integration testing of the mesh processing.
 */

#include "test_base.h"

/**	\test Performs integration testing for the mesh processing.
 *
 *	\todo Update comments.
 */
void test_integration_mesh
	(struct Test_Info*const test_info, ///< \ref Test_Info.
	 const char*const mesh_name        ///< The test mesh name.
	);

#endif // DPG__test_integration_mesh_h__INCLUDED
