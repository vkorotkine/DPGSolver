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
/** \file
 *  \brief Provides functions/structures for integration testing.
 */

#include <stdbool.h>

struct Test_Info;
struct Simulation;

/// Container for integration test related parameters
struct Integration_Test_Info {
	const char* ctrl_name; ///< The name of the control file used for the current test.

	/// The name of the extension for a specific convergence order test.
	const char* conv_study_extension;

	/** The minimal and maximal reference orders to be used for convergence order testing.
	 *
	 *  This range may differ from \ref Simulation::p_ref, such that convergence order testing can be performed
	 *  without calling the adaptation functions, but it is then expected that the corresponding control files for
	 *  each of the cases are present.
	 */
	const int p_ref[2];

	const int ml[2]; ///< The minimal and maximal mesh levels to be used for convergence order testing.

	const int adapt_type; ///< Analogue of \ref Simulation::adapt_type.

	const double conv_order_discount; ///< The allowed convergence order discount.
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

/// \brief 'Con'/'De'strucor for a \ref Simulation depending on `adapt_type`.
void structor_simulation
	(struct Simulation** sim,    ///< Pointer to the \ref Simulation.
	 const char mode,            ///< Mode of operation. Options: 'c'onstructor, 'd'estructor.
	 const int adapt_type,       ///< \ref Integration_Test_Info::adapt_type.
	 const int p,                ///< The order of the current simulation.
	 const int ml,               ///< The mesh level of the current simulation.
	 const int p_prev,           ///< The order of the previous simulation.
	 const int ml_prev,          ///< The mesh level of the previous simulation.
	 const char*const ctrl_name, ///< The name of the control file (used for \ref constructor_Simulation).
	 const char type_rc,         ///< The 'r'eal/'c'omplex type of the simulation.
	 const bool ignore_static    ///< Flag for whether if conditions based on static variables should be ignored.
	);

/** \brief Set the name of the current file to be used based on the adaptation type, order and mesh level.
 *  \return See brief. */
const char* set_file_name_curr
	(const int adapt_type,      ///< \ref Integration_Test_Info::adapt_type.
	 const int p,               ///< The order of the current simulation.
	 const int ml,              ///< The mesh level of the current simulation.
	 const bool add_missing,    ///< Flag to add missing 'ml' and 'p' parameters to the file name.
	 const char*const file_name ///< The name of the file.
	);

/// \brief Adapt the initial mesh if a adaptation parameters are provided in the input file.
void adapt_initial_mesh_if_required
	(struct Simulation*const sim ///< \ref Simulation.
	);

/// \brief Adapt to the maximum refinement levels in both h and p.
void adapt_to_maximum_refinement
	(struct Simulation*const sim,                           ///< Standard.
	 const struct Integration_Test_Info*const int_test_info ///< Standard.
	);

#endif // DPG__test_integration_h__INCLUDED
