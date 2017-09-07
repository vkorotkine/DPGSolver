// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

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
