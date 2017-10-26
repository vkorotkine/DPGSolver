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

#ifndef DPG__test_integration_h__INCLUDED
#define DPG__test_integration_h__INCLUDED
/**	\file
 *	\brief Provides functions/structures for integration testing.
 */

#include <stdbool.h>

struct Test_Info;

/// Container for integration test related parameters
struct Integration_Test_Info {
	const char* ctrl_name; ///< The name of the control file used for the current test.

	const bool display_progress; ///< \ref Test_Case::display_progress.

	/** The minimal and maximal reference orders to be used for convergence order testing.
	 *
	 *  This range may differ from \ref Simulation::p_ref, such that convergence order testing can be performed
	 *  without calling the adaptation functions, but it is then expected that the corresponding control files for
	 *  each of the cases are present.
	 */
	const int p_ref[2];

	const int ml[2]; ///< The minimal and maximal mesh levels to be used for convergence order testing.

	const int adapt_type; ///< Analogue of \ref Simulation::adapt_type.
};

/// \brief Call integration test functions.
void run_tests_integration
	(struct Test_Info*const test_info ///< \ref Test_Info.
	);

/** \brief Constructor for the \ref Integration_Test_Info from the data in the control file.
 *  \return Standard. */
struct Integration_Test_Info* constructor_Integration_Test_Info
	(const char*const ctrl_name ///< The control name excluding the full path.
	);

/// \brief Destructor for the \ref Integration_Test_Info.
void destructor_Integration_Test_Info
	(struct Integration_Test_Info* int_test_info ///< Standard.
	);

#endif // DPG__test_integration_h__INCLUDED
