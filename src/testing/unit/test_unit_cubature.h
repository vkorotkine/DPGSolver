// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_unit_cubature_h__INCLUDED
#define DPG__test_unit_cubature_h__INCLUDED
/** \file
 *  \brief Provides functionality for unit testing of the reference coordinates (and associated cubature if relevant).
 */

struct Test_Info;

/// \test Performs unit testing for the cubature.
void test_unit_cubature
	(struct Test_Info*const test_info ///< \ref Test_Info.
	);

#endif // DPG__test_unit_cubature_h__INCLUDED
