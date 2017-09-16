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

#ifndef DPG__test_support_h__INCLUDED
#define DPG__test_support_h__INCLUDED
/** \file
 *  \brief Provides base support functions for testing.
 */

#include <stdio.h>
#include <stdbool.h>

/// \brief Update the value of pass based on the old/new values.
void update_pass
	(bool*const pass,     ///< The old value for `pass` (Modified in the function if necessary).
	 const bool pass_new  ///< The new value for `pass`.
	);

/** \brief Compare the input `int`s and optionally print them if they differ.
 *  \return `true` if tests passed. */
bool compare_i
	(const int a,              ///< Input 0.
	 const int b,              ///< Input 1.
	 const bool print_enabled, ///< Flag for whether printing is enabled.
	 const char*const var_name ///< The variable name.
	);

/** \brief Compare the input `bool`s and optionally print them if they differ.
 *  \return `true` if tests passed. */
bool compare_b
	(const bool a,             ///< Input 0.
	 const bool b,             ///< Input 1.
	 const bool print_enabled, ///< Flag for whether printing is enabled.
	 const char*const var_name ///< The variable name.
	);

/** \brief Set the variable_name to be printed as part of one of the comparison functions based on a container member.
 *  \return See brief. */
char* set_print_name_container_member
	(const char*const name_container, ///< The name of the container.
	 int ind_container,               ///< The container index.
	 const char*const name_member     ///< The name of the member.
	);

/// \brief Check that the container type is that which is expected.
void check_container_type
	(FILE* data_file,                ///< The file containing the data.
	 const char*const container_type ///< The container type.
	);

/** \brief Check the entries of the input `differences` for `true` values.
 *  \return `true` if an entry in `differences` is `true`; `false` otherwise. */
bool check_diff
	(const int n_entries,          ///< The number of entries.
	 const bool*const differences, ///< The array of difference flags.
	 bool* pass                    ///< Flag for whether the test passed.
	);

#endif // DPG__test_support_h__INCLUDED
