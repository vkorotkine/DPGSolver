// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_base_h__INCLUDED
#define DPG__test_base_h__INCLUDED
/**	\file
 *	Provides base functions/structures for general testing.
 */

#include <stdbool.h>
#include <time.h>

// Containers ******************************************************************************************************* //

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

	clock_t ts, ///< Start time.
	        te; ///< End time.

	int n_test, ///< The number of tests run.
	    n_pass, ///< The number of tests which passed.
	    n_warn; ///< The number of warnings generated.

	struct Test_Integration_Run t_int; ///< \ref Test_Integration_Run.
};

// Interface functions ********************************************************************************************** //

/// \brief Increment test counters and print pass/fail related information for the current test.
void test_increment_and_print
	(struct Test_Info*const test_info, ///< \ref Test_Info.
	 const bool pass,                  ///< Flag for whether the test passed or failed.
	 const char*const test_name        ///< The test name string.
	);

/// \brief Print a warning and increment \ref Test_Info::n_warn.
void test_print_warning
	(struct Test_Info*const test_info, ///< \ref Test_Info.
	 const char*const warn_name        ///< The warning message string.
	);

/// \brief Print a failure message.
void test_print_failure
	(struct Test_Info*const test_info, ///< \ref Test_Info.
	 const char*const fail_name        ///< The failure message string.
	);

/// \brief Output test related information.
void output_test_info
	(struct Test_Info*const test_info ///< \ref Test_Info.
	);

#endif // DPG__test_base_h__INCLUDED
