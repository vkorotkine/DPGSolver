// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_unit_h__INCLUDED
#define DPG__test_unit_h__INCLUDED
/**	\file
 *	\brief Provides functions/structures for unit testing.
 */

struct Test_Info;

/// \brief Call unit test functions.
void run_tests_unit
	(struct Test_Info*const test_info ///< \ref Test_Info.
	);

#endif // DPG__test_unit_h__INCLUDED
