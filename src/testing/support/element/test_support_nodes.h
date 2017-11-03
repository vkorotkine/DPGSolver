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

#ifndef DPG__test_support_nodes_h__INCLUDED
#define DPG__test_support_nodes_h__INCLUDED
/** \file
 *  \brief Provides support functions for testing relating to the \ref Nodes container.
 */

#include <stdbool.h>

// Constructor functions ******************************************************************************************** //

/** \brief `const` version of \ref constructor_file_name_Nodes.
 *  \return Standard. */
const struct const_Nodes* constructor_file_name_const_Nodes
	(const char*const var_name,      ///< Defined for \ref constructor_file_name_Nodes.
	 const char*const file_name_full ///< Defined for \ref constructor_file_name_Nodes.
	);

/** \brief `const` version of \ref constructor_file_name_Plotting_Nodes.
 *  \return See brief. */
const struct const_Plotting_Nodes* constructor_file_name_const_Plotting_Nodes
	(const char*const var_name,      ///< The name of the variable to be read in from the file.
	 const char*const file_name_full ///< The name of the file (including the full path).
	);

// Difference functions ********************************************************************************************* //

/** \brief `const` version of \ref diff_Nodes.
 *  \return See brief. */
bool diff_const_Nodes
	(const struct const_Nodes*const a, ///< Input 0.
	 const struct const_Nodes*const b, ///< Input 1.
	 const double tol                  ///< The tolerance.
	);

/** \brief `const` version of \ref diff_const_Plotting_Nodes.
 *  \return See brief. */
bool diff_const_Plotting_Nodes
	(const struct const_Plotting_Nodes*const a, ///< Input 0.
	 const struct const_Plotting_Nodes*const b, ///< Input 1.
	 const double tol                           ///< The tolerance.
	);

// Printing functions *********************************************************************************************** //

/// \brief `const` version of \ref print_diff_Nodes.
void print_diff_const_Nodes
	(const struct const_Nodes*const a, ///< Input 0.
	 const struct const_Nodes*const b, ///< Input 1.
	 const double tol                  ///< The tolerance.
	);

/// \brief `const` version of \ref print_diff_Plotting_Nodes.
void print_diff_const_Plotting_Nodes
	(const struct const_Plotting_Nodes*const a, ///< Input 0.
	 const struct const_Plotting_Nodes*const b, ///< Input 1.
	 const double tol                           ///< The tolerance.
	);

#endif // DPG__test_support_nodes_h__INCLUDED
