// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_info_h__INCLUDED
#define DPG__test_info_h__INCLUDED
/**	\file
 *	Provides functions/structures for general testing.
 */

#include <stdbool.h>

/// Container for flags marking which integration are to be run.
struct Test_Integration_Run {
	bool equivalence_real_complex, ///< Test equivalence of real and complex functions.
	     equivalence_algorithms,   ///< Test equivalence of the supported algorithms.
	     linearization,            ///< Test the linearizations.
	     conv_order;               ///< Test for expected convergence orders.
};

/// Container for test related information.
struct Test_Info {
	int    nargc; ///< Standard.
	char** argv;  ///< Standard.

	int n_test, ///< The number of tests run.
	    n_pass, ///< The number of tests which passed.
	    n_warn; ///< The number of warnings generated.

	struct Test_Integration_Run t_int; ///< \ref Test_Integration_Run.
};

#endif // DPG__test_info_h__INCLUDED
